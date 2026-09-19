"""Delete only inventoried, untracked build artifacts; never scientific records."""
from pathlib import Path
import os,json,hashlib,subprocess,collections,time,zipfile
R=Path('/home/shadi/iclr2027');D=R/'docs/workspace_cleanup/20260914';P=R/'papers/interaction_recoverability_iclr2027'
start=time.monotonic()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=R).decode().split('\0'))-{''}
tracked_before={n:sha(R/n) for n in tracked if (R/n).is_file()}
status_before=subprocess.check_output(['git','diff','--name-status'],cwd=R,text=True)
candidates={};retained={}
def add(p,reason):
 if not p.is_file() or p.is_symlink():return
 n=str(p.relative_to(R))
 if n in tracked:return
 assert not n.startswith(('runs/','data/','review_packages/','docs/interaction_recoverability/'))
 candidates[n]=dict(reason=reason,bytes=p.stat().st_size,sha256=sha(p))
# Page crops and raster previews have no unique manuscript source content.
for v,sub in [('v1','private/renders'),('v3','private/extract'),('v4','private/extract'),('v6','rendered')]:
 folder=P/v/sub
 if not folder.exists():continue
 pdf=P/v/'main.pdf';assert pdf.is_file();retained[str(pdf.relative_to(R))]=sha(pdf)
 for p in folder.rglob('*'):
  if p.suffix.lower() in ['.png','.jpg','.jpeg','.pdf','.txt']:add(p,'temporary page extract/render; full manuscript PDF retained')
# Raster previews in private authoring areas, without traversing source-build trees.
for v in ['v1','v3','v4']:
 folder=P/v/'private'
 if not folder.exists():continue
 for p in folder.glob('*'):
  if p.suffix.lower() in ['.png','.jpg','.jpeg']:add(p,'temporary page/contact-sheet preview')
for sub in ['sheets','final_sheets']:
 folder=P/'v1/private/visual_revision'/sub
 if folder.exists():
  for p in folder.rglob('*'):
   if p.suffix.lower() in ['.png','.jpg','.jpeg']:add(p,'temporary contact-sheet preview')
for v in ['v7','v8']:
 audit=R/f'docs/sirna_gnn_empirical/manuscript_{v}'
 for p in audit.glob('page-*.png'):add(p,'temporary page-layout inspection image')
 if (audit/'rendered').exists():
  for p in (audit/'rendered').rglob('*'):
   if p.suffix.lower() in ['.png','.jpg','.jpeg']:add(p,'temporary page/contact-sheet preview')
 # Remove extracted source duplicates only after matching every source to ZIP bytes.
 archive=P/v/'manuscript_source.zip';retained[str(archive.relative_to(R))]=sha(archive)
 with zipfile.ZipFile(archive) as z:
  hashes={i.filename:hashlib.sha256(z.read(i)).hexdigest() for i in z.infolist() if not i.is_dir()}
 for name in ['independent_source','independent_source_final']:
  folder=audit/name
  if not folder.exists():continue
  for p in folder.rglob('*'):
   if not p.is_file():continue
   rel=str(p.relative_to(folder))
   if rel in hashes and sha(p)==hashes[rel]:add(p,'byte-identical source already retained in manuscript_source.zip')
   elif p.name in ['main.pdf','main.aux','main.bbl','main.blg','main.log','main.out']:
    # Logs remain useful evidence; preserve them outside the disposable checkout.
    if p.suffix=='.log':continue
    add(p,'generated artifact from completed independent source build; final PDF and audit retained')
# Ignore scientific/library dependency directories; only remove project bytecode caches.
for top in ['certmp','scmp','experiments','tests','interaction_recoverability','intervention','sirna_gnn_empirical']:
 folder=R/top
 if not folder.exists():continue
 for p in folder.rglob('*.pyc'):
  if '__pycache__' in p.parts and '_deps' not in p.parts:add(p,'regenerable Python bytecode cache')
if (R/'.pytest_cache').exists():
 for p in (R/'.pytest_cache').rglob('*'):add(p,'pytest cache')
# TeX intermediates are not final PDFs, bibliography sources, build logs or audit records.
for p in P.rglob('*'):
 if p.suffix in ['.aux','.blg','.out','.toc','.lof','.lot','.fls','.fdb_latexmk'] and '_renderdeps' not in p.parts:
  add(p,'regenerable TeX build intermediate')
# Protect the current deliverables, full historical manuscripts and scientific publication.
for v in P.iterdir():
 if v.is_dir():
  for name in ['main.pdf','manuscript_source.zip','workflow_handoff.zip']:
   p=v/name
   if p.is_file():retained[str(p.relative_to(R))]=sha(p)
for p in P.glob('*/**/*.tex'):
 n=str(p.relative_to(R))
 if n not in candidates:retained[n]=sha(p)
for name in ['files.json']:
 for v in ['v1','v3']:
  file=R/f'docs/sirna_gnn_empirical/publication/{v}/{name}'
  if file.exists():
   for n in json.loads(file.read_text()):
    assert n not in candidates,('published file selected',n)
assert len(candidates)>0
inventory=dict(authorization='User requested removal of unnecessary workspace files',entries=candidates,retained=retained,tracked_before=tracked_before,git_diff_before=status_before)
(D/'deletion_manifest.json').write_text(json.dumps(inventory,indent=2)+'\n')
counts=collections.Counter();sizes=collections.Counter()
# Restrict mutation to reviewed exact files, and detect intervening edits.
for n,x in candidates.items():
 p=R/n;assert sha(p)==x['sha256'];p.unlink();counts[x['reason']]+=1;sizes[x['reason']]+=x['bytes']
# Remove now-empty cache/extraction directories only, deepest first.
parents={parent for n in candidates for parent in (R/n).parents if parent!=R and R in parent.parents}
for folder in sorted(parents,key=lambda p:len(p.parts),reverse=True):
 try:folder.rmdir()
 except OSError:pass
assert all((R/n).is_file() and sha(R/n)==h for n,h in tracked_before.items())
assert all((R/n).is_file() and sha(R/n)==h for n,h in retained.items())
assert subprocess.check_output(['git','diff','--name-status'],cwd=R,text=True)==status_before
report=dict(deleted_files=len(candidates),deleted_bytes=sum(sizes.values()),deleted_MiB=sum(sizes.values())/1024**2,categories={k:dict(files=counts[k],bytes=sizes[k]) for k in counts},tracked_files_unchanged=len(tracked_before),protected_sources_and_deliverables_verified=len(retained),new_scientific_jobs=0,historical_ledgers_changed=False,wall_seconds=time.monotonic()-start)
(D/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
