from pathlib import Path
import json, hashlib, time
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=ROOT/(HERE/'active_run.txt').read_text().strip()
PARENT=ROOT/'runs/sirna_gnn_empirical/v4-20260914T071926Z'
DOC=ROOT/'docs/sirna_gnn_empirical/v5'
PAPER=ROOT/'papers/interaction_recoverability_iclr2027/v10'
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
