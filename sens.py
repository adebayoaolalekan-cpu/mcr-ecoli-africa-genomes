from base import *
d=load()
cats=['IncI2','IncX4','IncHI2','Mixed MOB-Recon bin','Other plasmid','Chromosomal or unresolved']
out={}
def dist(x,w=None):
    if w is None: w=pd.Series(1.0,index=x.index)
    t=w.groupby(x).sum(); tot=w.sum()
    return {c:(round(float(t.get(c,0)),1),round(100*float(t.get(c,0))/tot,1)) for c in cats}, round(float(tot),1)
out['all']=dist(d.rep)
e=d[d.BioProject!='PRJNA1273464']; out['excl_PRJNA1273464']=dist(e.rep)
w=1/d.groupby('BioProject').BioSample.transform('count'); out['bioproject_weighted']=dist(d.rep,w)
m=d.groupby('BioProject').rep.agg(lambda s:s.value_counts().index[0]); out['one_per_bioproject_modal']=dist(m)
for name,sub in [('poultry_all',d[d['Sector (derived)']=='Poultry']),('poultry_excl',e[e['Sector (derived)']=='Poultry']),('human_all',d[d['Sector (derived)']=='Human'])]:
    out[name]=dist(sub.rep)
pw=d[d['Sector (derived)']=='Poultry']; out['poultry_bpweighted']=dist(pw.rep,1/pw.groupby('BioProject').BioSample.transform('count'))
out['n_ST_all']=int(d['ST (Achtman) - to fill'].nunique()); out['n_ST_excl']=int(e['ST (Achtman) - to fill'].nunique())
out['n_countries_excl']=int(e.Country.nunique()); out['n_excl']=int(len(e))
out['sector_excl']=e['Sector (derived)'].value_counts().to_dict()
cl=json.load(open('../snp_clusters.json')); memb={b:c['erd'] for c in cl for b in c['ours']}
d['erd']=d.BioSample.map(memb).fillna(d.BioSample)
dd=d.drop_duplicates('erd'); out['one_per_snp_cluster']=dist(dd.rep); out['n_snpdedup']=int(len(dd))
out['n_bioprojects']=int(d.BioProject.nunique())
bp=d.groupby('BioProject').agg(n=('BioSample','count'),country=('Country',lambda s:'/'.join(sorted(set(str(x) for x in s)))),years=('Year',lambda s:(lambda v:(f"{int(v.min())}-{int(v.max())}" if len(v) and v.min()!=v.max() else (str(int(v.min())) if len(v) else '')))(s.dropna()))).sort_values('n',ascending=False)
bp.to_csv('bioproject_summary.csv')
json.dump(out,open('sens.json','w'),indent=1,default=str)
for k,v in out.items(): print(k, v if not isinstance(v,tuple) else '')
