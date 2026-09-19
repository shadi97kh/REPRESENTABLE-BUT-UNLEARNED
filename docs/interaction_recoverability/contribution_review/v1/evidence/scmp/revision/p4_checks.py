"""Measured no-training contract tests for the P4 addition."""
import datetime
import os
from pathlib import Path
for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
if hasattr(os,'sched_getaffinity'):
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
from .p4_resources import usage, usage_delta, sha, write_record, Budget


def main():
    start = usage()
    import torch
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    import pytest
    paths = ['tests/revision/test_p4_decision.py','tests/model/test_certified_model.py']
    with Budget(60,60).enforce():
        code = pytest.main(['-q',*paths])
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d-%H%M%S-%f')
    record = {'kind':'p4_tests','exit_code':int(code),'test_files':paths,
              'source_hashes':{str(p):sha(p) for p in Path('scmp/revision').glob('p4_*.py')},
              'test_hashes':{p:sha(p) for p in paths},'resources':usage_delta(start),
              'training_executed':False,'optimizer_steps':0}
    print(write_record(f'runs/revision/p4_tests-{stamp}.json',record))
    return int(code)


if __name__ == '__main__':
    raise SystemExit(main())
