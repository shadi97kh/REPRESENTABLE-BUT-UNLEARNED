from common import *
import zipfile,shutil,sys,io
D=RUN/'delivery';D.mkdir(exist_ok=True)
parent=PARENT/'delivery/robustness_code_results.zip'
assert sha(parent)=='bec0c6e4282cc6e7d3b7967fc4f7b27a8695263d006df2eec679732c4f9a5ead'
readme='''# Contribution-gated revision — manuscript v10 / investigation v5

Start with manuscript/main_text_9_pages.pdf and audits/execution_summary.md.
G1 FAIL; G2 UNRESOLVED; G3 blocked and unexecuted. No new ML method is established.

The outer archive contains every new code, exact-reference result, proof, primary comparison source, audit and final manuscript asset. parent_evidence.zip is the preserved, self-contained v4 scientific package: original primary data, v3/v4 predictions and tables, all 792 focused fits/checkpoints, forty reused checkpoints, configurations and executed histories. It is intentionally stored intact rather than silently rewriting historical evidence. Its manifest indexes all 16,930 members. v3/v4 scientific CSV paths in the paper are relative to its extracted root. Unrelated historical checkpoint binaries remain in the separate preserved v3 original full archive and are not needed for this revision's comparisons.

Extract this archive into a FRESH directory. Run `python restore_current.py` to extract the parent archive and restore the original relative scientific paths. This writes files only inside that fresh directory. The parent archive's README/manuscript are historical; current paper artifacts are in current_manuscript/ after restoration. Current audits are in current_audits/. No training is launched by restoration.

Clean editable source: manuscript/manuscript_source.zip contains exactly main.tex, appendix.tex, references.bib and figures/. Unzip it into its own empty source directory. Unzip manuscript/typesetting_dependencies.zip separately. Run `bash <typesetting>/build.sh <source> <external-output>`. Official style files are outside source. The source ZIP was independently compiled.

The new runner `bash sirna_gnn_empirical/v5/run_all.sh` enforces contribution gates, verifies completed hashes, and blocks the absent positive pilot. Historical fit code exists to document what ran; do not mistake its presence for authorization to repeat old campaigns. Exact toy examples are synthetic correctness checks, never siRNA measurements. No new fitting, wet lab, push or external contact occurred in v5.

This is a private local review/evidence bundle, not a de-identified conference submission. Original machine logs and paths are intentionally preserved. The manuscript PDF and source themselves are anonymous. Final packaging cost/integrity receipts are supplied alongside because an archive cannot contain its own final hash.
'''
restore='''from pathlib import Path
import zipfile, subprocess, sys, shutil
root=Path(__file__).resolve().parent
shutil.copy2(root/'evidence_manifest.json',root/'current_evidence_manifest.json')
# Preserve current navigation before expanding the intact historical package.
for src,dst in [('manuscript','current_manuscript'),('audits','current_audits')]:
 if not (root/dst).exists():shutil.copytree(root/src,root/dst)
with zipfile.ZipFile(root/'parent_evidence.zip') as z:
 for n in z.namelist():
  p=Path(n)
  if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe archive member')
 z.extractall(root)
subprocess.run([sys.executable,str(root/'restore_layout.py')],cwd=root,check=True)
r=root/'runs/sirna_gnn_empirical/v5-20260914T091218Z';r.parent.mkdir(parents=True,exist_ok=True)
if not r.exists():r.symlink_to(root/'v5',target_is_directory=True)
d=root/'docs/sirna_gnn_empirical/v5';d.parent.mkdir(parents=True,exist_ok=True)
if not d.exists():d.symlink_to(root/'current_audits',target_is_directory=True)
p=root/'papers/interaction_recoverability_iclr2027/v10';p.mkdir(parents=True,exist_ok=True)
for archive,folder in [('manuscript_source.zip','source'),('typesetting_dependencies.zip','typesetting')]:
 (p/folder).mkdir(exist_ok=True)
 with zipfile.ZipFile(root/'current_manuscript'/archive) as z:z.extractall(p/folder)
print('Restored current and historical evidence; no scientific jobs executed.')
'''
files={}
def addtree(base,prefix,skip=None):
 for p in sorted(base.rglob('*')):
  if p.is_file() and '__pycache__' not in p.parts and not (skip and skip(p)):files[prefix+'/'+str(p.relative_to(base))]=p
addtree(HERE,'sirna_gnn_empirical/v5')
addtree(RUN,'v5',lambda p:p.is_relative_to(D) or p.name=='runner.lock')
addtree(DOC,'audits')
for p in (PAPER/'artifacts').glob('*'):
 if p.is_file() and p.suffix in ['.pdf','.zip','.json']:files['manuscript/'+p.name]=p
addtree(PAPER/'artifacts/figures','manuscript/figures')
files['parent_evidence.zip']=parent
manifest=[dict(path=n,bytes=p.stat().st_size,sha256=sha(p)) for n,p in sorted(files.items())]
archive=D/'contribution_code_evidence.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
 z.writestr('README_CURRENT.md',readme);z.writestr('restore_current.py',restore);z.writestr('evidence_manifest.json',json.dumps(manifest,indent=2))
 for n,p in sorted(files.items()):z.write(p,n,compress_type=zipfile.ZIP_STORED if p.suffix=='.zip' else zipfile.ZIP_DEFLATED)
# Verify every written file by streaming decompressed bytes and checking CRC/hash.
import hashlib
with zipfile.ZipFile(archive) as z:
 for row in manifest:
  with z.open(row['path']) as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
  assert actual==row['sha256'],row['path']
write(D/'archive_verification.json',dict(passed=True,archive=str(archive.relative_to(ROOT)),bytes=archive.stat().st_size,sha256=sha(archive),verified_members=len(manifest),parent_archive_sha256=sha(parent),scope='All new artifacts plus intact nested complete focused evidence archive'))
# Compact opinion packet keeps current paper, decisions and interpretable exports.
R=ROOT/'review_packages/sirna_gnn_review_v4';R.mkdir(exist_ok=True)
(R/'manuscript').mkdir(exist_ok=True);(R/'tables').mkdir(exist_ok=True);(R/'research').mkdir(exist_ok=True)
for name in ['main.pdf','main_text_9_pages.pdf','manuscript_source.zip','typesetting_dependencies.zip']:shutil.copy2(PAPER/'artifacts'/name,R/'manuscript'/name)
for p in DOC.glob('*'):
 if p.is_file():shutil.copy2(p,R/'research'/p.name)
for p in (RUN/'tables').glob('*.csv'):shutil.copy2(p,R/'tables'/p.name)
shutil.copy2(RUN/'correctness/results.csv',R/'tables/correctness.csv');shutil.copy2(RUN/'reuse/recomputed_metrics.csv',R/'tables/recomputed_metrics.csv')
shutil.copy2(PARENT/'adjudication/row_reconciliation.csv',R/'tables/source_adjudication.csv')
shutil.copy2(PARENT/'evaluation/ensemble_predictions.csv',R/'tables/ensemble_predictions.csv')
shutil.copytree(PAPER/'artifacts/figures',R/'figures',dirs_exist_ok=True)
(R/'REVIEW_REQUEST.md').write_text('''# Request for an independent opinion — manuscript v10

Draft for the user; no message has been sent. Read manuscript/main_text_9_pages.pdf, then the complete 53-page manuscript/main.pdf. Main pages 1–9, statements10, references11–14, appendix15–53. The workflow is vector artwork on page2.

Please evaluate the empirical contribution without assuming a novel GNN or a positive chemistry result. The new investigation found G1 FAIL, G2 UNRESOLVED and G3 blocked: its finite-response candidate reduces to standard linear optimization and does not justify true-response membership after learned encoding. Nineteen correctness checks passed; zero new fits ran. The existing 792-fit robustness campaign and original negative results are preserved.

Assess (1) whether the exact reduction and scoped theorem are correct; (2) whether any identified primary method changes that conclusion under comparable assumptions; (3) whether the source/label sensitivities and conditional S7 findings support the diagnostic contribution; (4) whether B2/B3 interpretations avoid biological, causal or independent-study overclaims; (5) whether the paper's remaining limitations are clear. Separate confirmed errors from speculative concerns.

Return the strongest defensible contribution, three ranked objections, precise corrections and readiness recommendation. Use research/appendix_audit.md and EVIDENCE_MAP.md. The complete archive is supplied separately; this compact packet does not contain all checkpoint binaries. Do not infer their absence from the original workspace.
''')
(R/'EVIDENCE_MAP.md').write_text('''# Current evidence map

- manuscript/: complete paper, nine-page extract, editable four-entry ZIP and separate official typesetting dependencies.
- research/gate_decision.md and novelty_matrix.md: promotion decisions and primary-source comparison.
- research/method_and_proof.md and theorem_to_code.md: complete conditional proof, counterexamples and actual implementation scope.
- tables/main_comparison.csv, activity_metrics.csv, direct_comparisons.csv and seed_stability.csv: all reused sensitivity scores and separate uncertainty summaries.
- tables/recomputed_metrics.csv: independent saved-prediction rescoring and ensemble/mean-seed identity.
- tables/source_adjudication.csv: every original source row, source cells and unresolved classifications.
- tables/B2_unchanged.csv and B3_visibility.csv: measured chemistry outcomes and negative controls.
- tables/ensemble_predictions.csv: actual held-out sensitivity predictions.
- figures/workflow.pdf/.svg/.png: caption-free current workflow.
- research/appendix_audit.md and claim_audit.md: every main object and scientific claim mapped to complete appendix support or precise full-archive rows.

The separately supplied contribution_code_evidence.zip includes the intact parent_evidence.zip containing primary inputs, 792 focused fit/checkpoint directories, 40 reused checkpoints and full original B2/B3 records. The outer evidence_manifest.json and the nested original manifest identify every file and hash. No new fit/checkpoint is implied by v5 directory names.
''')
compact=ROOT/'review_packages/sirna_gnn_review_v4.zip'
with zipfile.ZipFile(compact,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(R.rglob('*')):
  if p.is_file():z.write(p,str(p.relative_to(R)))
with zipfile.ZipFile(compact) as z:assert z.testzip() is None
write(D/'compact_verification.json',dict(path=str(compact.relative_to(ROOT)),bytes=compact.stat().st_size,sha256=sha(compact),crc_verified=True))
print('Complete archive',archive.stat().st_size,'bytes; compact review',compact.stat().st_size,'bytes')
