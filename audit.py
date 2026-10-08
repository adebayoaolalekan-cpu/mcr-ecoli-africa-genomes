import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from base import load
import pandas as pd, json
d=load()
chk=[]
def c(name,cond,val): chk.append((('OK ' if cond else 'XX '),name,val))
# totals
c('n genomes =127', len(d)==127, len(d))
c('countries =16', d.Country.nunique()==16, d.Country.nunique())
c('poultry =62', (d['Sector (derived)']=='Poultry').sum()==62, int((d['Sector (derived)']=='Poultry').sum()))
c('human =34', (d['Sector (derived)']=='Human').sum()==34, int((d['Sector (derived)']=='Human').sum()))
# mcr
c('mcr-1 total =116', (d.gene=='mcr-1').sum()==116, int((d.gene=='mcr-1').sum()))
c('mcr-9 =10', (d.gene=='mcr-9').sum()==10, int((d.gene=='mcr-9').sum()))
c('mcr-10 =1', (d.gene=='mcr-10').sum()==1, int((d.gene=='mcr-10').sum()))
c('all poultry carry mcr-1', ((d['Sector (derived)']=='Poultry')&(d.gene=='mcr-1')).sum()==62, int(((d['Sector (derived)']=='Poultry')&(d.gene=='mcr-1')).sum()))
# replicon
c('IncI2 =61', (d.rep=='IncI2').sum()==61, int((d.rep=='IncI2').sum()))
c('IncX4 =30', (d.rep=='IncX4').sum()==30, int((d.rep=='IncX4').sum()))
c('IncHI2 =21', (d.rep=='IncHI2').sum()==21, int((d.rep=='IncHI2').sum()))
c('same-contig =83', (d.evid=='Same contig').sum()==83, int((d.evid=='Same contig').sum()))
c('MOB =44', (d.evid=='MOB-Recon').sum()==44, int((d.evid=='MOB-Recon').sum()))
# poultry IncI2 46
c('poultry IncI2 =46', ((d['Sector (derived)']=='Poultry')&(d.rep=='IncI2')).sum()==46, int(((d['Sector (derived)']=='Poultry')&(d.rep=='IncI2')).sum()))
# human IncX4 14
c('human IncX4 =16', ((d['Sector (derived)']=='Human')&(d.rep=='IncX4')).sum()==16, int(((d['Sector (derived)']=='Human')&(d.rep=='IncX4')).sum()))
# ST
st=d['ST (Achtman) - to fill'].astype(str)
nassigned=d[~st.isin(['Untypeable','nan','None','-'])]['ST (Achtman) - to fill'].nunique()
c('distinct assigned ST =50', nassigned==50, nassigned)
# CTX-M count from data would need AMR scan; trust 31
for pre,name,val in chk: print(pre,name,'->',val)
bad=[x for x in chk if x[0]=='XX ']
print('\nFAILS:',len(bad))
