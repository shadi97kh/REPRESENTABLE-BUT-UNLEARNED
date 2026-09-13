"""Real measured-label forward/backward profile, gradient and reload checks."""
import argparse,json,time,resource,os,copy
from pathlib import Path
import numpy as np
import torch
from common import write_json
from models import build,tensors,take
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
torch.set_num_threads(2);torch.manual_seed(701);device='cuda:0' if torch.cuda.is_available() else 'cpu'
if device.startswith('cuda'):torch.cuda.set_per_process_memory_fraction(.30,0);torch.cuda.reset_peak_memory_stats()
m=json.loads((r/'graph_manifest.json').read_text());z=np.load(r/'train_graphs.npz');data={k:torch.tensor(z[k][:128],device=device) for k in ['X','A','mask','C']};y=torch.tensor(np.load(r/'train_targets.npz')['y'][:128],device=device)
model=build('gnn',m,{'hidden':32}).to(device);initial={k:v.detach().clone() for k,v in model.state_dict().items()};opt=torch.optim.AdamW(model.parameters(),lr=.003,weight_decay=1e-4)
if device.startswith('cuda'):torch.cuda.synchronize()
start=time.perf_counter();cpu=time.process_time();history=[];gradients=[]
for step in range(20):
 model.train();opt.zero_grad(set_to_none=True);loss=((model(**data)-y)**2).mean();assert torch.isfinite(loss);loss.backward();gradients.append(sum(float(p.grad.norm()) for n,p in model.named_parameters() if 'message_weights' in n and p.grad is not None));torch.nn.utils.clip_grad_norm_(model.parameters(),5.);opt.step();history.append(float(loss.detach()))
if device.startswith('cuda'):torch.cuda.synchronize()
elapsed=time.perf_counter()-start;cpu_used=time.process_time()-cpu
changed=sum(float((v-initial[k]).abs().sum()) for k,v in model.state_dict().items() if 'message_weights' in k);assert changed>0 and max(gradients)>0
model.eval()
with torch.no_grad():before=model(**data).cpu()
p=r/'profile_checkpoint.pt';torch.save({'model':model.state_dict(),'optimizer':opt.state_dict(),'optimizer_updates':20,'kind':'gnn','config':{'hidden':32},'profile_only':True},p)
reload=build('gnn',m,{'hidden':32}).to(device);reload.load_state_dict(torch.load(p,map_location=device,weights_only=False)['model']);reload.eval()
with torch.no_grad():after=reload(**data).cpu()
err=float((before-after).abs().max());assert err==0
record={'device':device,'torch':torch.__version__,'cuda_runtime':torch.version.cuda,'optimizer_updates':20,'batch_size':len(y),'training_rows_only':True,'wall_s':elapsed,'cpu_s':cpu_used,'seconds_per_update':elapsed/20,'loss_history':history,'message_gradient_norms':gradients,'message_parameter_absolute_change':changed,'checkpoint_reload_max_error':err,'parameters':sum(p.numel() for p in model.parameters()),'cuda_peak_allocated_bytes':torch.cuda.max_memory_allocated() if device.startswith('cuda') else 0,'profile_weights_reused_for_final_training':False}
write_json(r/'profile.json',record);write_json(r/'profile.outputs.json',['profile.json','profile_checkpoint.pt']);print(json.dumps(record,indent=2))
