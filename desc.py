from base import *
d=load()
o={}
# year bins
def yb(y):
    if pd.isna(y): return 'Not recorded'
    y=int(y)
    return '2005-2014' if y<=2014 else ('2015-2018' if y<=2018 else ('2019-2021' if y<=2021 else '2022-2025'))
d['yb']=d.Year.map(yb)
ct=pd.crosstab(d.yb,d.rep); ct['n']=ct.sum(1); print(ct)
ctg=pd.crosstab(d.yb,d.gene); print(ctg)
ctc=pd.crosstab(d.yb,d.Country); print(ctc)
print(d[d.Year.isna()][['BioSample','BioProject','Country','Collection date (as submitted)']])
print(d[d.Year<=2012][['BioSample','BioProject','Country','Year','Host (as submitted)','Isolation source (as submitted)','mcr variant(s)','ST (Achtman) - to fill','Plasmid replicon carrying mcr - to fill','Galaxy input ID']].to_string())
# sector mapping
sm=d.groupby(['Host (as submitted)','Isolation source (as submitted)','Sector (derived)','Poultry niche (derived)'],dropna=False).size().reset_index(name='n')
sm.to_csv('sector_mapping.csv',index=False); print(sm.to_string())
