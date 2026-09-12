"""Preparation-only CLI. No training command exists."""
import argparse
import hashlib
import json
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['diagnose'])
    parser.add_argument('--out',required=True,help='New exclusive output directory')
    args=parser.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    from .diagnostics import run
    result=run();p=out/'diagnostics.json'
    with p.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    (out/'diagnostics.json.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest()+'  diagnostics.json\n')
    for row in result['rows']:print(json.dumps({k:row[k] for k in ('delta','primary_interaction_gap','axis2_w1','joint_w1','axis2_hellinger_squared')}))

if __name__=='__main__':main()
