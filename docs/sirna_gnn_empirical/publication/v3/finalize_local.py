from pathlib import Path
import hashlib,json,shutil,subprocess,datetime
R=Path('/home/shadi/iclr2027');W=Path('/tmp/iclr2027-sirna-v3-publication');D=R/'docs/sirna_gnn_empirical/publication/v3';paths=json.loads((D/'allowlist.json').read_text());before={}
# Mirror only reviewed publication files, backing up any differing local document.
for name in paths:
 src=W/name;dst=R/name
 if dst.exists():
  if src.read_bytes()==dst.read_bytes():continue
  before[name]=hashlib.sha256(dst.read_bytes()).hexdigest();backup=D/'preserved_local_before'/name;backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dst,backup)
 dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
(D/'preserved_local_before.json').write_text(json.dumps(before,indent=2)+'\n')
head=subprocess.check_output(['git','rev-parse','origin/main'],cwd=R,text=True).strip();expected=subprocess.check_output(['git','rev-parse','HEAD'],cwd=W,text=True).strip();assert head==expected
receipt=dict(status='Published; remote main independently fetched and verified',repository='https://github.com/shadi97kh/iclr2027',branch='main',commit=head,commit_url=f'https://github.com/shadi97kh/iclr2027/commit/{head}',changed_files=len(paths),previous_remote_commit='05e600976eecff800960d8d603799d9d8eb89b6c',local_main_preserved=subprocess.check_output(['git','rev-parse','main'],cwd=R,text=True).strip(),publication_checkout=str(W),force_push=False,scientific_fits=0,optimizer_updates=0,published_scope='Scientific code, configuration, aggregate results and support records, supplied workflow PNG/PDF, eleven result figures with PDF/PNG/SVG exports, READMEs and integrity manifest',excluded='Manuscript TeX/PDF, raw inputs, per-observation prediction files, checkpoint binaries, review archives, caches and authoring instructions',administrative_cost='Nonzero inspection, copying, hashing, compilation and network cost; not separately metered. No historical ledger changed.',recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(D/'status.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
