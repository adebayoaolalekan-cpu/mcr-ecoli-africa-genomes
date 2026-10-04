import sys; sys.path.insert(0,'.')
from base import load
import pandas as pd, json
d=load()
chk=[]
def c(name,cond,val): chk.append((('OK ' if cond else 'XX '),name,val))
# totals
c('n genomes =118', len(d)==118, len(d))
c('countries =13', d.Country.nunique()==13, d.Country.nunique())
c('poultry =62', (d['Sector (derived)']=='Poultry').sum()==62, int((d['Sector (derived)']=='Poultry').sum()))
c('human =30', (d['Sector (derived)']=='Human').sum()==30, int((d['Sector (derived)']=='Human').sum()))
# mcr
c('mcr-1 total =109', (d.gene=='mcr-1').sum()==109, int((d.gene=='mcr-1').sum()))
c('mcr-9 =8', (d.gene=='mcr-9').sum()==8, int((d.gene=='mcr-9').sum()))
c('mcr-10 =1', (d.gene=='mcr-10').sum()==1, int((d.gene=='mcr-10').sum()))
c('all poultry carry mcr-1', ((d['Sector (derived)']=='Poultry')&(d.gene=='mcr-1')).sum()==62, int(((d['Sector (derived)']=='Poultry')&(d.gene=='mcr-1')).sum()))
# replicon
c('IncI2 =59', (d.rep=='IncI2').sum()==59, int((d.rep=='IncI2').sum()))
c('IncX4 =28', (d.rep=='IncX4').sum()==28, int((d.rep=='IncX4').sum()))
c('IncHI2 =20', (d.rep=='IncHI2').sum()==20, int((d.rep=='IncHI2').sum()))
c('same-contig =79', (d.evid=='Same contig').sum()==79, int((d.evid=='Same contig').sum()))
c('MOB =39', (d.evid=='MOB-Recon').sum()==39, int((d.evid=='MOB-Recon').sum()))
# poultry IncI2 46
c('poultry IncI2 =46', ((d['Sector (derived)']=='Poultry')&(d.rep=='IncI2')).sum()==46, int(((d['Sector (derived)']=='Poultry')&(d.rep=='IncI2')).sum()))
# human IncX4 14
c('human IncX4 =14', ((d['Sector (derived)']=='Human')&(d.rep=='IncX4')).sum()==14, int(((d['Sector (derived)']=='Human')&(d.rep=='IncX4')).sum()))
# ST
st=d['ST (Achtman) - to fill'].astype(str)
nassigned=d[~st.isin(['Untypeable','nan','None','-'])]['ST (Achtman) - to fill'].nunique()
c('distinct assigned ST =48', nassigned==48, nassigned)
# CTX-M count from data would need AMR scan; trust 26
for pre,name,val in chk: print(pre,name,'->',val)
bad=[x for x in chk if x[0]=='XX ']
print('\nFAILS:',len(bad))
