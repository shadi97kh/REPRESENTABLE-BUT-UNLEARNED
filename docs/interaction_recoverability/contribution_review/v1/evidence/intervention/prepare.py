"""Small exact fixtures and timing checks, with no optimizer steps or model zoo."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from .accounting import snapshot,elapsed,sha,write_json


def main(out):
    start=snapshot()
    import os
    import resource
    import time
    import numpy as np
    import torch
    from .tasks import family,all_actions,control_predictors
    from .operator import GraphConditionalProcess,response_loss
    from .transforms import FrozenQueryOracle,molecular_input,digest
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    allowed=os.sched_getaffinity(0)
    os.sched_setaffinity(0,{min(allowed)})
    records=[];profiles={}
    before=snapshot()
    for task in ('edge_switch','node_state'):
        fam=family(task,0);base=fam.transform()
        for name,f in control_predictors(task).items():
            model_hash=digest(dict(task=task,coefficients=f.coefficients,intercept=f.intercept))
            oracle=FrozenQueryOracle(f,model_hash,base.construction)
            for a in all_actions():
                state=fam.transform(a)
                effect=oracle.effect(base,state)
                records.append(dict(graph_id=base.identity,task=task,model=name,model_hash=model_hash,
                                    model_kind='independent_analytic_function_NOT_GNN',split='correctness_fixture',
                                    action=json.dumps(a),action_identity=state.identity,
                                    response_kind='model_query_NOT_biological_outcome',
                                    reference_output=oracle.query(base),edited_output=oracle.query(state),effect=effect))
    profiles['80_exact_analytic_queries_plus_graph_construction']=elapsed(before)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    p=out/'responses.csv'
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    Path(str(p)+'.sha256').write_text(f'{sha(p)}  {p.name}\n')
    torch.manual_seed(91226)
    m=GraphConditionalProcess().double()
    fam=family('edge_switch');base=fam.transform()
    targets=[fam.transform(a) for a in all_actions()]
    q=targets[1:4];y=torch.tensor([.2,-.3,.4],dtype=torch.float64)
    state={k:v.clone() for k,v in m.state_dict().items()}
    # Five single-episode forward/backward timings, with NO optimizer creation.
    timings=[]
    for _ in range(5):
        t=snapshot();m.zero_grad(set_to_none=True)
        response_loss(m(base,q,y,targets),torch.zeros(8,dtype=torch.float64)).backward()
        timings.append(elapsed(t))
    assert all(torch.equal(v,state[k]) for k,v in m.state_dict().items())
    profiles['untrained_operator_forward_backward_per_episode']=timings
    t=snapshot()
    for _ in range(5): molecular_input('CCN(CC)CCOC(=O)CCCC')
    profiles['five_small_molecule_feature_rebuilds']=elapsed(t)
    source_rows=list(csv.DictReader(Path('data/intervention_sources/agile.csv').open()))
    t=snapshot()
    # Source-order first five molecules; no model or assay label is evaluated.
    for row in source_rows[:5]: molecular_input(row['combined_mol_SMILES'])
    profiles['five_actual_lipid_feature_rebuilds']=elapsed(t)
    cpus=[t['process_cpu_s'] for t in timings[1:]]
    walls=[t['wall_s'] for t in timings[1:]]
    # Explicit extrapolation, not an authorized allowance or measured training.
    count=13500*8
    estimate=dict(source_predictor_training_cpu_core_hours=0,source_predictor_training_gpu_hours=0,
                  source_generation='analytic coefficients only in prospective CPU screening',
                  cold_start_cpu_s=timings[0]['process_cpu_s'],
                  cold_start_charged_once_per_run=True,
                  forward_backward_episode_count=count,
                  core_hours_from_measured_median=count*float(np.median(cpus))/3600,
                  wall_hours_from_measured_median=count*float(np.median(walls))/3600,
                  planning_core_hour_range=[count*min(cpus)/3600,count*max(cpus)*3/3600],
                  query_cache_cpu_hours_linear=(16384+1024+4096)/80*profiles['80_exact_analytic_queries_plus_graph_construction']['cpu_core_hours'],
                  gpu_hours=0,cpu_threads=1,workers=0,
                  unknowns=['Adam update time and optimizer state peak memory not measured',
                            'auxiliary loss and eight-episode execution overhead only allowed for by 3x upper planning factor',
                            'SPEX external adapter runtime and development baseline selection not measured',
                            'published delivery-predictor inference, full descriptors, 3D construction and fitting not measured',
                            'end-to-end training convergence, runtime variance and generalization unknown'],
                  no_training_authorization=True)
    write_json(out/'profile.json',dict(training_steps=0,predictor_training_steps=0,main_teacher_cache_generated=False,
                                     correctness_fixture_query_count=len(records),profiles=profiles,
                                     parameters=sum(p.numel() for p in m.parameters()),
                                     architecture_equivalence='generic_graph_CNP is same class',
                                     estimates=estimate,resources=elapsed(start)))
    print(json.dumps(estimate,indent=2));print(json.dumps(elapsed(start)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',default='runs/intervention/preparation-v1')
    args=p.parse_args();main(args.out)
