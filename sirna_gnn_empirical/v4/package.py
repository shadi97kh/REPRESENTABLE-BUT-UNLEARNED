"""Indexed focused evidence archive; reuse old source bytes without duplicating the 6GB archive."""
from common import *
import zipfile,shutil,sys
P=ROOT/'papers/interaction_recoverability_iclr2027/v9';A=P/'artifacts';D=ROOT/'docs/sirna_gnn_empirical/v4';dest=RUN/'delivery';dest.mkdir(exist_ok=True)
assert (A/'independent_build.json').exists() and (RUN/'audit/scientific_audit.json').exists()
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';files={}
def addtree(base,prefix,accept=lambda p:True):
 for p in sorted(base.rglob('*')):
  if p.is_file() and not p.is_symlink() and accept(p.relative_to(base)):files[prefix+'/'+str(p.relative_to(base))]=p
addtree(RUN,'v4',lambda p:p.parts[0]!='delivery' and p.name!='campaign.lock')
addtree(OLD,'v1',lambda p:p.parts[0] not in ['fits','paper_artifacts','superseded'] and p.suffix not in ['.zip','.pt'] and p.name!='campaign.lock')
addtree(PREV,'v2',lambda p:p.parts[0] not in ['fits','superseded'] and p.suffix not in ['.zip','.pt'] and p.name!='campaign.lock')
methods=['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn']
def v3keep(p):
 if p.parts[0] in ['paper_artifacts','legacy','authoring_versions','package_verification']:return False
 if p.suffix=='.zip' or p.name=='campaign.lock':return False
 if p.parts[0]=='fits' and p.suffix in ['.pt','.pkl']:
  return len(p.parts)>4 and p.parts[1:4]==('activity','outer0','final') and p.parts[4] in methods
 return True
addtree(V3,'v3',v3keep)
for v in ['v1','v2','v3','v4']:addtree(ROOT/'sirna_gnn_empirical'/v,'sirna_gnn_empirical/'+v,lambda p:'__pycache__' not in p.parts)
addtree(D,'audits');addtree(P/'typesetting','typesetting')
for n in ['main.pdf','main_text_9_pages.pdf','manuscript_source.zip','typesetting_dependencies.zip','page_boundaries.json','independent_build.json','proof_preservation.json','minimal_source_validation.json','minimal_zip_validation.json','validation_commands.json']:
 files['manuscript/'+n]=A/n
files['instruction/Codex_Targeted_Robustness_and_Minimal_TeX.md']=ROOT/'Codex_Targeted_Robustness_and_Minimal_TeX.md'
# Include exact raw workbook dependencies referred to by adjudication; upstream code/licenses remain as acquired.
manifest=[dict(path=n,bytes=p.stat().st_size,sha256=sha(p)) for n,p in sorted(files.items())]
write(dest/'evidence_manifest.json',manifest)
readme='''# Focused robustness evidence — manuscript v9 / campaign v4

Start with manuscript/main_text_9_pages.pdf, then manuscript/main.pdf and audits/execution_summary.md.
All scientific tables retain full precision. v4/ is the current run; v3/ is the preserved parent evidence,
v1/ and v2/ supply exact historical input/prediction dependencies. The manuscript's evidence index uses
these exact paths. No release label, historical prediction or checkpoint was overwritten.

The archive contains all 792 new fit objects with best/last checkpoints, optimizer states and histories;
40 reused original outer0 final models; all original v3 fit metadata, histories, raw fit/gradient/selection
tables and predictions; raw measured inputs and primary-source/citation evidence; executable campaign
code; source adjudication; all comparison tables; and the four-entry manuscript ZIP with external style files.
Other original v3 checkpoint binaries remain in the preserved full original archive
runs/sirna_gnn_empirical/v3-20260914T013314Z/paper_artifacts/code_results.zip and original workspace.
They are not needed for the focused refits or saved-prediction analyses; this bundle does not duplicate
that approximately 6GB archive. Older v1/v2 fits are likewise preserved locally, not relabeled as new work.

Verify evidence_manifest.json hashes before use. To restore the original relative run layout on a clean
extraction, run `python restore_layout.py`. Then use `bash sirna_gnn_empirical/v4/run_all.sh`.
Completed stage receipts verify immutable outputs before reuse; no checkpoint verification is a new fit.
Presentation/build outputs were assembled after scientific completion. Absolute command paths record
the original execution host; the restore helper supplies equivalent relative data locations.

The manuscript source ZIP contains exactly main.tex, appendix.tex, references.bib and figures/.
Unzip it separately and compile with typesetting/build.sh; do not put style or build products in the source.
No public upload, GitHub push or external communication was performed. These are retrospective measured-data
prediction studies, not new wet-lab or independent external validation. Review files absent from the provided
workspace are identified in audits/review_response.md, rather than claimed read.
'''
restore='''from pathlib import Path
root=Path(__file__).resolve().parent
runs=root/'runs/sirna_gnn_empirical';runs.mkdir(parents=True,exist_ok=True)
versions={'v1':'v1-20260913T185027Z','v2':'v2-20260913T220847Z','v3':'v3-20260914T013314Z','v4':'v4-20260914T071926Z'}
for alias,name in versions.items():
 p=runs/name
 if not p.exists():p.symlink_to(Path('../../')/alias,target_is_directory=True)
import zipfile,shutil
paper=root/'papers/interaction_recoverability_iclr2027/v9';paper.mkdir(parents=True,exist_ok=True)
source=paper/'source';source.mkdir(exist_ok=True)
with zipfile.ZipFile(root/'manuscript/manuscript_source.zip') as z:z.extractall(source)
if not (paper/'typesetting').exists():shutil.copytree(root/'typesetting',paper/'typesetting')
print('Restored run aliases and minimal source; preserved scientific bytes are unchanged.')
'''
archive=dest/'robustness_code_results.zip';t=time.monotonic()
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=2,allowZip64=True) as z:
 z.writestr('README.md',readme);z.writestr('restore_layout.py',restore);z.write(dest/'evidence_manifest.json','evidence_manifest.json')
 for name,p in sorted(files.items()):z.write(p,name)
# Streaming integrity check every member; do not extract another multi-GB copy.
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in manifest:
  h=hashlib.sha256()
  with z.open(r['path']) as f:
   for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
  assert h.hexdigest()==r['sha256'],r['path']
write(dest/'package_verification.json',dict(passed=True,archive=str(archive.relative_to(ROOT)),sha256=sha(archive),archive_bytes=archive.stat().st_size,manifest_members=len(manifest),uncompressed_bytes=sum(x['bytes'] for x in manifest),all_member_hashes_verified=True,wall_s=time.monotonic()-t,scope='Complete focused experiment dependencies and original prediction/fit-record evidence; unrelated old checkpoint binaries remain in the preserved full parent archive.'))
print((dest/'package_verification.json').read_text())
