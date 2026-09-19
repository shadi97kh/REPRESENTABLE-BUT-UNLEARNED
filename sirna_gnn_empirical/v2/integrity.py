from core import *
import re
D=ROOT/'docs/sirna_gnn_empirical/v2';D.mkdir(exist_ok=True);before=json.loads((RUN/'preservation_before.json').read_text());changed=[s for s,m in before.items() if not (ROOT/s).is_file() or sha(ROOT/s)!=m['sha256']];assert not changed,changed
fs=[json.loads(p.read_text()) for p in (RUN/'fits').rglob('fit.json')];assert len(fs)==210
for f in fs:assert sha(RUN/f['checkpoint'])==f['checkpoint_sha256']
assert len(fs)+4<=240 and sum(f['optimizer_updates'] for f in fs)+41<400000
# Verify no source/sequence component crosses roles within each outer fit.
for file,group in [('grouped_membership.csv','study_group'),('grouped_membership.csv','sequence_group'),('b2_membership.csv','sequence_group')]:
 df=pd.read_csv(RUN/file);assert (df.groupby(['outer_fold',group]).role.nunique()==1).all()
q=pd.read_csv(RUN/'grouped_predictions.csv');assert q.groupby(['model','seed','record_id']).size().max()==1
b=pd.read_csv(RUN/'b2_predictions.csv');assert b.groupby(['model','seed','pair_id']).size().max()==1
write(D/'integrity_audit.json',dict(preexisting_files_verified=len(before),changed=changed,engine_fits=len(fs),additional_direct_pair_ridge_fits=4,checkpoint_hashes_verified=len(fs),recorded_optimizer_updates=sum(f['optimizer_updates'] for f in fs),profile_updates=20,possible_uncheckpointed_interrupted_updates=[0,21],selected_state_updates=sum(f['selected_updates'] for f in fs),group_role_checks_pass=True,duplicate_outer_prediction_check_pass=True,scientific_code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}))
write(RUN/'integrity_complete.json',dict(preexisting_files_verified=len(before),changed=changed,fit_checkpoint_hashes_verified=len(fs)))
print('Verified',len(before),'preserved files and',len(fs),'fit checkpoints')
