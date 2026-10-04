import pandas as pd, collections
rows=[l.rstrip('\n').split('\t') for l in open(__import__('sys').argv[1]) if l.strip()]
A=[r for r in rows if r[0]=='A']; P=[r for r in rows if r[0]=='P']; R=[r for r in rows if r[0]=='R']; M=[r for r in rows if r[0]=='M']
mcr=collections.defaultdict(list)
for r in A:
    if r[3].startswith('mcr'): mcr[r[1]].append((r[3],r[2],r[6],r[7],r[8]))
pf=collections.defaultdict(list)
for r in P: pf[r[1]].append((r[2],r[3],r[4],r[5]))
st={}
for r in M:
    st[r[1]]=(r[3],' '.join(r[4:]))
genomes=sorted(set(st))
out=[]
for g in genomes:
    calls=mcr.get(g,[])
    links=[]
    for sym,contig,meth,cov,idt in calls:
        reps=[p[1] for p in pf.get(g,[]) if p[0]==contig]
        links.append(f"{sym}@{contig}:{'/'.join(reps) if reps else 'NOMATCH'}")
    out.append(dict(genome=g,ST=st[g][0],alleles=st[g][1],mcr=';'.join(c[0]+'('+c[2]+','+c[3]+'/'+c[4]+')' for c in calls),link=' | '.join(links),reps=','.join(sorted(set(p[1] for p in pf.get(g,[]))))))
df=pd.DataFrame(out); df.to_csv(__import__('sys').argv[2],sep='\t',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_colwidth',90)
print(df[['genome','ST','mcr','link']].to_string())
