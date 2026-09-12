"""No-fit checkpoint validation and post-run sampling diagnostics on GPU 1."""
import os
os.environ['CUDA_VISIBLE_DEVICES']='GPU-1e502523-2587-176b-91ef-86cb0ea01a82'
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import argparse
from pathlib import Path
import numpy as np
import torch
from ..model import CertifiedModel
from ..rna.ensemble import DeclaredGibbsPrior
from .p4_capability import state_hash, ensemble_forward
from .p4_gpu_math import FastPrior, row_features, pack, adjacencies, batched_forward
from .p4_resources import checked_record, sha, write_record, usage, usage_delta, Budget


def main(path):
    start=usage()
    torch.set_num_threads(1); torch.set_num_interop_threads(1)
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
    torch.use_deterministic_algorithms(True)
    assert torch.cuda.device_count()==1
    torch.cuda.set_per_process_memory_fraction(.15,0)
    device=torch.device('cuda:0')
    cycle=checked_record(path)
    inputs=checked_record(cycle['artifacts']['rna_inputs']['path'])
    rna=checked_record(cycle['artifacts']['rna']['path'])
    checks={}
    with Budget(120,120).enforce():
        for name in ['capability_checkpoint','rna_A_checkpoint','rna_A_plus_B_checkpoint']:
            record=cycle['artifacts'][name]
            assert sha(record['path'])==record['sha256']
            cp=torch.load(record['path'],map_location='cpu',weights_only=True)
            m=CertifiedModel(cp['d_in'],cp['candidates'],cp['n_nodes'],d_hid=cp['hidden'],gamma=cp['gamma']).to(device)
            m.load_state_dict(cp['state_dict'])
            assert state_hash(m)==record['model_sha256']==cp['model_sha256']
            checks[name]={'sha256':record['sha256'],'state_hash_matches':True,'parameters_finite':all(bool(torch.isfinite(p).all()) for p in m.parameters())}
        before=state_hash(m)
        ids=inputs['split']['train'][:16]
        rows=[inputs['rows'][i] for i in ids]
        priors=[FastPrior(r['window']) for r in rows]
        fs=[row_features(r,r['window']) for r in rows]
        packed=pack(fs,priors,device)
        labels=torch.tensor([r['efficacy'] for r in rows],device=device,dtype=torch.float64)
        # Validate the actual 50-nt domain against the original marginal implementation.
        p=priors[0]
        exact=np.asarray(DeclaredGibbsPrior(p.oracle,{e:.5 for e in p.edges}).marginals().mu)
        marginal_error=float(np.max(abs(exact-p.mu)))
        assert marginal_error<1e-10
        As=adjacencies(priors,8,np.random.RandomState(903),device)
        batch=batched_forward(m,packed,As)
        sequential=[]
        for i,p in enumerate(priors):
            local=CertifiedModel(90,p.edges,50,d_hid=12,gamma=m.gamma).to(device)
            local.anchor,local.g=m.anchor,m.g
            sequential.append(ensemble_forward(local,packed[0][i],torch.as_tensor(p.mu,device=device),None,As[i]))
        error=float((batch-torch.stack(sequential)).abs().max().detach().cpu())
        assert error<1e-10
        # Same saved checkpoint and original prediction seed must reproduce saved predictions.
        from .p4_gpu_rna import predict
        selected=inputs['split']['test'][:16]
        rr=[inputs['rows'][i] for i in selected]
        ff=[row_features(r,r['window']) for r in rr]
        pp=[FastPrior(r['window']) for r in rr]
        repro=[]
        for i,f,p in zip(selected,ff,pp):
            value,_=predict(m,[f],[p],64,411+30000+i,Budget(30,30),device)
            repro.append(float(value[0]))
        prediction_error=float(np.max(abs(np.array(repro)-np.asarray(rna['predictions']['actual_A_plus_B'])[selected])))
        assert prediction_error<1e-10
        chosen=m.gamma
        m.set_gamma(1.)  # B training scale, solely a frozen-gradient diagnostic; restored below.
        diagnostics={}
        for count in [8,32]:
            gradients,losses=[],[]
            for seed in [1901,1902,1903,1904]:
                A=adjacencies(priors,count,np.random.RandomState(seed),device)
                pred=batched_forward(m,packed,A)
                loss=(pred-labels).square().mean()
                grad=torch.autograd.grad(loss,tuple(m.g.parameters()))
                gradients.append(torch.cat([g.reshape(-1) for g in grad]).detach().cpu().numpy())
                losses.append(float(loss.detach().cpu()))
            g=np.stack(gradients)
            norm=float(np.linalg.norm(g.mean(axis=0)))
            rms_sd=float(np.sqrt(g.var(axis=0,ddof=1).sum()))
            diagnostics[str(count)]={'draw_replicates':4,'mean_gradient_l2':norm,'gradient_rms_standard_deviation':rms_sd,
               'relative_rms_sd':rms_sd/norm if norm else None,'losses':losses}
        m.set_gamma(chosen)
        assert state_hash(m)==before
        torch.cuda.synchronize()
    resources=usage_delta(start)
    resources.update(gpu_hours=resources['wall_s']/3600,gpu_devices=1,
       physical_gpu_uuid=os.environ['CUDA_VISIBLE_DEVICES'],gpu_charge_basis='conservative whole postcheck wall time',
       cuda_peak_allocated_bytes=torch.cuda.max_memory_allocated(0))
    result={'cycle':path,'checkpoint_checks':checks,'exact_50nt_marginal_max_error':marginal_error,
       'batched_50nt_output_max_error':error,'saved_prediction_reproduction_max_error':prediction_error,
       'gradient_sampling':diagnostics,'gradient_sampling_rows':ids,'gradient_sampling_gamma':1.,
       'diagnostic_status':'post-run sampling diagnostic only; no model/seed/hyperparameter selection',
       'gradient_sampling_scope':'first 16 training rows, frozen final checkpoint, four independent draws per count; conditional sampling variability, not optimization or biological variance',
       'frozen_checkpoint_unchanged':state_hash(m)==before,'optimizer_steps':0,'resources':resources,
       'source_sha256':sha(__file__)}
    dest=path.replace('-cycle.json','-postcheck.json')
    write_record(dest,result)
    print(dest); print(result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('cycle')
    main(parser.parse_args().cycle)
