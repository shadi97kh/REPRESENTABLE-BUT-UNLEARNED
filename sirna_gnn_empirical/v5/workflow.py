"""Caption-free, native vector rendering of the verified existing architecture."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(16,9.45));fig.subplots_adjust(0,0,1,1);ax.set(xlim=(0,16),ylim=(0,9.45));ax.axis('off')
ink='#24394A';muted='#5E7080';green='#E3F2EC';purple='#ECE5F5';blue='#E4F0F9';peach='#FAEBDC'
def box(x,y,w,h,text,color,size=18,bold=False):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.05,rounding_size=0.12',lw=1.4,edgecolor='#BBC8D3',facecolor=color))
 ax.text(x+w/2,y+h/2,text,ha='center',va='center',color=ink,fontsize=size,fontweight='bold' if bold else 'normal',linespacing=1.35)
def arrow(x,y,u,v,color=muted,style='-'):
 ax.add_patch(FancyArrowPatch((x,y),(u,v),arrowstyle='-|>',mutation_scale=18,lw=1.8,color=color,linestyle=style))
# Four scientific panels; every node describes an actual forward or training operation.
for x,w,title,col in [(0.1,3.15,'a  Encode duplex',green),(3.5,5.0,'b  Typed message passing',purple),(8.75,3.15,'c  Ordered readout',blue),(12.15,3.75,'d  Evaluate',peach)]:
 box(x,2.8,w,6.5,'','#FBFCFE');box(x,8.68,w,.62,title,col,20,True)
box(.28,7.32,2.79,1.03,'Sequence + chemistry\nGuide / passenger',green,19)
box(.28,5.67,2.79,1.25,'X: 64 × 93\nA: 8 × 64 × 64\nPadding mask: 64',green,19)
arrow(1.68,7.26,1.68,6.98)
box(.28,4.27,2.79,1.03,'Training-only\nfeature/context masks',green,18)
arrow(1.68,5.59,1.68,5.36)
box(.28,3.03,2.79,.83,'Eight context fields',green,18)
arrow(3.28,9.0,3.45,9.0)
box(3.72,7.32,4.55,1.02,'Linear 93 → 32; ReLU; mask\nH⁰: 64 × 32',purple,19)
box(3.72,5.95,4.55,.98,'mᵢ = Σᵣⱼ Aᵣᵢⱼ hⱼ W̃ᵣ / dᵢ\ndᵢ = max(1, total degree)',purple,19)
arrow(6,7.25,6,7.04)
box(3.72,4.48,4.55,1.1,'Self transform + message\nReLU → dropout 0.1\nResidual → LayerNorm → mask',purple,18)
arrow(6,5.88,6,5.67)
ax.text(6,4.18,'Repeat twice · hidden width 32',ha='center',va='center',fontsize=18,color=ink,fontweight='bold')
box(3.72,2.98,4.55,.88,'R1: shared directional bases\n+ train-supported residuals',purple,18)
arrow(8.55,9.0,8.7,9.0)
box(8.98,7.24,2.69,1.15,'Flatten H² in order\n2,048 coordinates\n+ 8 context fields',blue,18)
box(8.98,5.68,2.69,1.15,'Dense 2,056 → 64\nReLU + dropout 0.1',blue,18)
arrow(10.33,7.17,10.33,6.93)
box(8.98,4.35,2.69,.89,'Linear 64 → 1: z',blue,19)
arrow(10.33,5.61,10.33,5.35)
box(8.98,3.05,2.69,.89,'ŷ = μT + sT z',blue,24)
arrow(10.33,4.28,10.33,4.04)
arrow(11.98,9.0,12.1,9.0)
box(12.38,7.32,3.29,1.02,'Freeze weights\nDisable dropout',peach,19)
box(12.38,5.93,3.29,1.03,'B1: activity / ranking\nOne observed condition',peach,18)
arrow(14.03,7.25,14.03,7.06)
box(12.38,4.52,3.29,1.03,'B2: ŷ₁ − ŷ₀\nSame sequence + context',peach,18)

box(12.38,3.02,3.29,1.12,'B3: ŷ₁₁ − ŷ₁₀ − ŷ₀₁ + ŷ₀₀\nFour measured conditions\nShared endpoint identities',peach,17)

# Distinct training paths; both end at updates of the same model weights.
box(.1,.15,7.65,2.25,'',blue);box(8,.15,7.9,2.25,'',peach)
ax.text(.35,2.03,'ACTIVITY TRAINING',fontsize=21,fontweight='bold',color=ink)
box(.35,.67,2.28,.94,'Training activity\nỹ = (y − μT)/sT',blue,17)
box(3.0,.67,2.07,.94,'Weighted MSE\n(z − ỹ)²',blue,18)
box(5.42,.67,2.05,.94,'Update shared\nmodel weights',blue,17)
arrow(2.7,1.14,2.91,1.14);arrow(5.12,1.14,5.34,1.14)
ax.text(8.25,2.03,'PAIR-SUPERVISED FOLLOW-UP',fontsize=21,fontweight='bold',color=ink)
box(8.25,.67,2.1,.94,'Training pairs only\nΔy = y₁ − y₀',peach,17)
box(10.67,.67,2.37,.94,'Activity loss + λ\nweighted pair loss',peach,17)
box(13.37,.67,2.27,.94,'Update the same\nshared weights',peach,17)
arrow(10.42,1.14,10.61,1.14);arrow(13.11,1.14,13.3,1.14)
ax.text(4,.38,'Separate protocol · measured training labels',ha='center',fontsize=15,color=muted)
ax.text(12,.38,'Held-out measured labels enter evaluation only',ha='center',fontsize=15,color=muted)
D=PAPER/'artifacts/figures';D.mkdir(parents=True,exist_ok=True)
for ext in ['pdf','svg','png']:fig.savefig(D/('workflow.'+ext),dpi=220,facecolor='white')
plt.close(fig)
print('Created native vector PDF/SVG and raster preview')
