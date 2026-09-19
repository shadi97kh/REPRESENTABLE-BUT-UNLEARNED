"""Versioned, fully executed P4 cycle under the latest single-GPU instruction."""
from __future__ import annotations
import os
# Pin by UUID before importing any tensor library. The physical GPU is cuda:0 here.
os.environ['CUDA_VISIBLE_DEVICES'] = 'GPU-1e502523-2587-176b-91ef-86cb0ea01a82'
for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[variable] = '1'
os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'

import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import time
import traceback
import yaml
import torch
from .p4_resources import Budget, usage, usage_delta, write_record, checked_record, sha
from .p4_gpu_math import validate
from .p4_gpu_control import run_control
from .p4_gpu_rna import run_rna
from .p4_verifier import run_verifier
from .p4_capability import state_hash


def main():
    start = usage()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    torch.use_deterministic_algorithms(True)
    cfg_path = Path('configs/revision_p4_gpu_v1.yaml')
    cfg = yaml.safe_load(cfg_path.read_text())
    c = cfg['execution']
    assert os.environ['CUDA_VISIBLE_DEVICES'] == c['gpu_uuid']
    assert torch.cuda.is_available() and torch.cuda.device_count() == 1
    torch.cuda.set_device(0)
    torch.cuda.set_per_process_memory_fraction(c['max_gpu_memory_fraction'], 0)
    device = torch.device('cuda:0')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
    prefix = Path('runs/revision') / ('p4_gpu-'+stamp)
    artifacts = {}
    def save_stage(name, result):
        path = str(prefix)+'-'+name+'.json'
        artifacts[name] = {'path':write_record(path,result), 'sha256':sha(path)}
        print(f'saved {name}: {path}',flush=True)
    def save_model(name, model):
        path = Path(str(prefix)+'-'+name+'.pt')
        with path.open('xb') as f:
            torch.save({'state_dict':{k:v.detach().cpu() for k,v in model.state_dict().items()},
                        'gamma':model.gamma,'candidates':model.candidates,'n_nodes':model.n_nodes,
                        'd_in':model.anchor.d_in,'hidden':model.g.lin0.out_features,
                        'model_sha256':state_hash(model)},f)
        Path(str(path)+'.sha256').write_text(f'{sha(path)}  {path.name}\n')
        artifacts[name+'_checkpoint'] = {'path':str(path),'sha256':sha(path),'model_sha256':state_hash(model)}
    gpu_before = subprocess.check_output(['nvidia-smi','--query-gpu=index,uuid,name,memory.total,memory.used,utilization.gpu',
                                           '--format=csv,noheader'],text=True)
    sources = [cfg_path] + sorted(Path('scmp/revision').glob('p4*.py')) + [Path('scmp/model.py'),Path('scmp/rna/ensemble.py'),Path('scmp/oracles/rna_noncrossing.py')]
    source_hashes = {str(p):sha(p) for p in sources}
    historical = {name:{'path':path,'sha256':sha(path)} for name,path in cfg['historical_records'].items()}
    for rec in historical.values():
        checked_record(rec['path'])
    audit = checked_record(c['row_record'])
    pre = {'created_utc':stamp,'config':cfg,'source_hashes':source_hashes,'historical_records':historical,
           'audited_rows_sha256':sha(c['row_record']),'gpu_preflight':gpu_before,
           'physical_gpu_uuid':c['gpu_uuid'],'visible_gpu_count':torch.cuda.device_count(),
           'torch_version':torch.__version__,'cpu_affinity':list(os.sched_getaffinity(0)),
           'outcomes_seen_before_execution':False,'historical_p3_previously_inspected':True,
           'note':'All requested comparisons execute. Capability and utility gates remain necessary for a positive verification claim.'}
    save_stage('prerun',pre)
    result = {'status':'RUNNING','started_utc':stamp,'config_path':str(cfg_path), 'config_sha256':sha(cfg_path),
              'physical_gpu_uuid':c['gpu_uuid'],'visible_gpu_count':1,'stage_wall_s':{}}
    budget = Budget(c['max_cpu_seconds'],c['max_wall_seconds'])
    try:
        with budget.enforce():
            t = time.perf_counter()
            valid = validate(device)
            torch.cuda.synchronize()
            save_stage('validation',valid)
            result['stage_wall_s']['validation'] = time.perf_counter()-t
            print('START capability training on physical GPU 1',flush=True)
            t = time.perf_counter()
            control, model = run_control(cfg['capability'],budget,device)
            torch.cuda.synchronize()
            save_stage('capability',control)
            save_model('capability',model)
            result['stage_wall_s']['capability'] = time.perf_counter()-t
            result['capability_passed'] = control['capability_passed']
            print(f'CAPABILITY GATE: {control["capability_passed"]}; continuing all requested comparisons without changing settings',flush=True)
            frozen = copy.deepcopy(model).cpu()
            print('START RNA actual-model and matched baseline experiment',flush=True)
            t = time.perf_counter()
            rna, rna_model = run_rna(cfg,audit['rows'],budget,device,save_stage,save_model)
            torch.cuda.synchronize()
            save_stage('rna',rna)
            result['stage_wall_s']['rna'] = time.perf_counter()-t
            print('START frozen-predictor verifier comparison (CPU bounds, same trained GPU checkpoint)',flush=True)
            t = time.perf_counter()
            verification = run_verifier(frozen,cfg,budget)
            verification['capability_passed'] = control['capability_passed']
            verification['useful_gain_claim_admissible'] = bool(control['capability_passed'] and verification.get('useful_learned_gain'))
            assert state_hash(frozen) == control['model_sha256']
            save_stage('verifier',verification)
            result['stage_wall_s']['verifier'] = time.perf_counter()-t
            result['status'] = 'COMPLETE_ALL_REQUESTED_COMPARISONS'
            result['optimizer_steps'] = control['optimizer_steps']+rna['optimizer_steps']
            result['useful_verification_gain'] = verification['useful_gain_claim_admissible']
    except Exception as exc:
        result.update(status='INCOMPLETE_EXECUTION',error=repr(exc),traceback=traceback.format_exc())
        print(result['traceback'],flush=True)
    finally:
        torch.cuda.synchronize()
        result['source_hashes_unchanged'] = all(sha(path)==value for path,value in source_hashes.items())
        result['artifacts'] = artifacts
        result['resources'] = usage_delta(start)
        result['resources'].update(gpu_hours=result['resources']['wall_s']/3600,
            gpu_charge_basis='conservative entire execution wall time on exactly one visible GPU',
            gpu_devices=1,physical_gpu_uuid=c['gpu_uuid'],
            cuda_peak_allocated_bytes=torch.cuda.max_memory_allocated(0),cuda_peak_reserved_bytes=torch.cuda.max_memory_reserved(0),
            cuda_memory_fraction_cap=c['max_gpu_memory_fraction'],
            child_processes_spawned=1,child_process_description='read-only nvidia-smi preflight; waited and included in child_cpu_s',
            historical_resource_balance='unreconciled; no claimed balance',cpu_affinity=list(os.sched_getaffinity(0)))
        path = str(prefix)+'-cycle.json'
        write_record(path,result)
        print(json.dumps({'cycle':path,'status':result['status'],'resources':result['resources']}),flush=True)
    return 0 if result['status']=='COMPLETE_ALL_REQUESTED_COMPARISONS' else 2


if __name__ == '__main__':
    raise SystemExit(main())
