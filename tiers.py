from base import *
d=load()
# mcr contig length for the 39 MOB-Recon assignments
clen={}
for l in open('galaxy_typing.tsv'):
    p=l.rstrip().split('\t')
    if p[0]=='B': clen[p[1]]=int(p[8])
# map Galaxy input id -> contig length
d['gid']=d['Galaxy input ID']
d['mcrlen']=d['gid'].map(clen)
# assembly level (closed) from PD asm_level
lvl=dict(zip(d['BioSample'],d.get('asm_level','')))
def tier(r):
    s=str(r['Plasmid replicon carrying mcr - to fill'])
    closed = str(r['gid'])=='GCA_007993675.1' or str(r['gid'])=='GCA_007993675.2'
    if 'MOB-Recon' not in s and 'MOB-recon' not in s and 'MOB' not in s:
        return 'T1 closed/complete' if closed else 'T2 same contig (draft)'
    # MOB-Recon assignment
    L=r['mcrlen']
    if pd.isna(L): return 'T3 MOB-Recon contig >=5 kb'
    return 'T3 MOB-Recon contig >=5 kb' if L>=5000 else 'T4 MOB-Recon contig <5 kb'
d['tier']=d.apply(tier,axis=1)
print(d['tier'].value_counts().to_string())
print('--- tier x replicon ---')
print(pd.crosstab(d['rep'],d['tier']).to_string())
print('--- T4 genomes (tentative) ---')
t4=d[d.tier.str.startswith('T4')]
print(t4[['BioSample','gid','gene','rep','mcrlen']].to_string())
# save per-genome tier for supp
d[['BioSample','gid','gene','rep','evid','mcrlen','tier']].to_csv('evidence_tiers.csv',index=False)
