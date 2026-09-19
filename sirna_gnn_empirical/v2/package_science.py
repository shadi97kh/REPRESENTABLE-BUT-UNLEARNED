"""Package exact existing scientific artifacts; does not evaluate or train."""
from core import *
import zipfile,subprocess
D=ROOT/'docs/sirna_gnn_empirical/v2';P=ROOT/'papers/interaction_recoverability_iclr2027/v3'
files=[]
for folder in [ROOT/'sirna_gnn_empirical/v2',ROOT/'sirna_gnn_empirical/v1',RUN,OLD,D]:
 for p in folder.rglob('*'):
  if not p.is_file() or p.is_symlink():continue
  rel=p.relative_to(ROOT);parts=rel.parts
  if any(s in parts for s in ['__pycache__','_deps','paper_artifacts']):continue
  if (p.suffix=='.zip' and 'raw' not in parts) or p.suffix=='.pyc' or p.name in ['campaign.lock','scientific_package_manifest.json','package_cost.json','scientific_archive_verification.json','deliverable_hashes.json']:continue
  # Prior archive's scientific raw inputs, saved models, metrics and source identities are all retained.
  files.append(p)
files=sorted(set(files));manifest={p.relative_to(ROOT).as_posix():dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files}
write(D/'scientific_package_manifest.json',manifest)
readme='''This is the completed targeted siRNA follow-up, with code, exact measured inputs, source identities, grouping, every fit and selected/last checkpoint, optimizer histories, held-out predictions, tables, figures, uncertainty and command logs. It also contains the unchanged v1 scientific inputs/checkpoints required for frozen diagnostics and comparison. No TeX manuscript source is in this archive; use manuscript_source.zip for that.

Extract into one root. Completed scientific stages verify using:
  bash sirna_gnn_empirical/v2/run_all.sh audit source_checks plan factorial grouped b2 external external_uncertainty evaluate figures
The exact default interpreter path is recorded in run_all.sh. On another machine, supply the recorded compatible environment and adjust that path explicitly. All local checkpoints were trained from scratch; no RNA feature backend is required. To execute new fits rather than verify completed ones, choose a fresh run location and preserve this run; do not delete completed flags/checkpoints to manufacture a new blind test. APP is retrospective; B2 is within-cohort supervision; Davis S7 is qualified separate-source endpoint assessment. The manuscript and preservation audit require their corresponding source/full original workspace and are separate administrative stages.

No publication, push, wet-lab validation, B3 interaction measurement or GNN architectural novelty is implied. Pair endpoints/control covariance and Davis S1/dose-time gaps remain explicit. See docs/sirna_gnn_empirical/v2/execution_summary.md and the full CSV comparison tables.
'''
start=time.perf_counter();cpu=time.process_time()
with zipfile.ZipFile(RUN/'code_results.zip','w',zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
 z.writestr('README_SCIENTIFIC_PACKAGE.txt',readme);z.writestr('MANIFEST.json',json.dumps(manifest,indent=2))
 for p in files:z.write(p,p.relative_to(ROOT))
write(D/'package_cost.json',dict(wall_s=time.perf_counter()-start,process_cpu_s=time.process_time()-cpu,files=len(files),archive_bytes=(RUN/'code_results.zip').stat().st_size,status='Administrative compression/hashing, separate from experimental stage totals'))
print('scientific archive',len(files),(RUN/'code_results.zip').stat().st_size)
