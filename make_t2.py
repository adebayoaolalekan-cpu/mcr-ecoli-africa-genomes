import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from base import load
import pandas as pd, json
d=load()
clen={}
for l in open('galaxy_typing.tsv'):
    p=l.rstrip().split('\t')
    if p[0]=='B': clen[p[1]]=int(p[8])
d['mcrlen']=d['Galaxy input ID'].map(clen)
def tier(r):
    s=str(r['Plasmid replicon carrying mcr - to fill'])
    if 'MOB' not in s:
        return 'T1' if str(r['Galaxy input ID'])=='GCA_007993675.2' else 'T2'
    L=r['mcrlen']
    return 'T3' if (pd.isna(L) or L>=5000) else 'T4'
d['tier']=d.apply(tier,axis=1)
def variant(r):
    v=str(r['mcr (Galaxy AMRFinderPlus)']).strip()
    if v.startswith('mcr-9'): return 'mcr-9 group'
    if v.startswith('mcr-10'): return 'mcr-10.1'
    if v=='mcr-1.5': return 'mcr-1.5'
    if v=='mcr-1.1': return 'mcr-1.1'
    return 'mcr-1 (no allele)'
d['var']=d.apply(variant,axis=1)
order=['IncI2','IncX4','IncHI2','Mixed MOB-Recon bin','Other plasmid','Chromosomal or unresolved']
vorder=['mcr-1.1','mcr-1 (no allele)','mcr-1.5','mcr-9 group','mcr-10.1']
hdr=['Replicon linked to mcr','n (%)','Closed','Same contig','MOB ≥ 5 kb','MOB < 5 kb']+vorder
rows=[hdr]
for r in order:
    sub=d[d.rep==r]; n=len(sub)
    rows.append([r,f"{n} ({round(100*n/118,1)})"]+
        [str((sub.tier==t).sum()) for t in ['T1','T2','T3','T4']]+
        [str((sub['var']==v).sum()) for v in vorder])
n=len(d)
rows.append(['Total',f"{n} (100)"]+[str((d.tier==t).sum()) for t in ['T1','T2','T3','T4']]+
            [str((d['var']==v).sum()) for v in vorder])
T=json.load(open('tables.json'))
T['T2']=rows
json.dump(T,open('tables.json','w'),indent=0)
for x in rows: print(x)
