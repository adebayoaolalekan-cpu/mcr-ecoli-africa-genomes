import sys; sys.path.insert(0,'.')
from base import load
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

d=load()
cats=['IncI2','IncX4','IncHI2','Other']
def catmap(x): return x if x in ('IncI2','IncX4','IncHI2') else 'Other'
d['rc']=d['rep'].map(catmap)
colors={'IncI2':'#2c6fbb','IncX4':'#e08214','IncHI2':'#5aae61','Other':'#9e9e9e'}

def stack(ax, groups, glabel):
    labels=list(groups.keys())
    x=np.arange(len(labels)); w=0.62
    bottoms=np.zeros(len(labels))
    for c in cats:
        solid=np.array([groups[g].get((c,'Same contig'),0) for g in labels])
        mob  =np.array([groups[g].get((c,'MOB-Recon'),0) for g in labels])
        ax.bar(x,solid,w,bottom=bottoms,color=colors[c],edgecolor='white',linewidth=0.6)
        bottoms=bottoms+solid
        ax.bar(x,mob,w,bottom=bottoms,color=colors[c],edgecolor='white',linewidth=0.6,hatch='////')
        bottoms=bottoms+mob
    ax.set_xticks(x); ax.set_xticklabels(labels,rotation=35,ha='right',fontsize=9)
    ax.set_ylabel('Genomes',fontsize=10)
    ax.set_title(glabel,fontsize=11,loc='left')
    for s in ['top','right']: ax.spines[s].set_visible(False)

def build_groups(col, order):
    g={}
    for key in order:
        sub=d[d[col]==key]
        gg={}
        for _,r in sub.iterrows():
            k=(r['rc'],r['evid']); gg[k]=gg.get(k,0)+1
        g[key]=gg
    return g

secorder=['Poultry','Human','Food (non-poultry)','Other animal','Environment','Unspecified']
seclabels={'Food (non-poultry)':'Food','Other animal':'Animal','Unspecified':'Unspec.'}
gs=build_groups('Sector (derived)',secorder)
gs={seclabels.get(k,k):v for k,v in gs.items()}

cc=d['Country'].value_counts()
countryorder=[c for c in cc.index if cc[c]>=4]
gc=build_groups('Country',countryorder)

fig,axes=plt.subplots(1,2,figsize=(10.2,4.4),gridspec_kw={'width_ratios':[1.05,1]})
stack(axes[0],gs,'A. By source sector')
stack(axes[1],gc,'B. By country (n ≥ 4)')

legend=[Patch(facecolor=colors[c],label=c) for c in cats]
legend+=[Patch(facecolor='white',edgecolor='#444',label='on mcr contig'),
         Patch(facecolor='white',edgecolor='#444',hatch='////',label='MOB-Recon')]
axes[1].legend(handles=legend,fontsize=8,loc='upper right',frameon=False,ncol=1)
plt.tight_layout()
plt.savefig('/home/claude/ms/figure1.png',dpi=200,bbox_inches='tight')
print('saved')
