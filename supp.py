import json, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
import sys; sys.path.insert(0,'.')
from base import load
d=load()

AU=["Algeria","Angola","Benin","Botswana","Burkina Faso","Burundi","Cameroon","Cape Verde","Central African Republic","Chad","Comoros","Congo (Republic)","Cote d'Ivoire","Democratic Republic of the Congo","Djibouti","Egypt","Equatorial Guinea","Eritrea","Eswatini","Ethiopia","Gabon","Gambia","Ghana","Guinea","Guinea-Bissau","Kenya","Lesotho","Liberia","Libya","Madagascar","Malawi","Mali","Mauritania","Mauritius","Morocco","Mozambique","Namibia","Niger","Nigeria","Rwanda","Sao Tome and Principe","Senegal","Seychelles","Sierra Leone","Somalia","South Africa","South Sudan","Sudan","Tanzania","Togo","Tunisia","Uganda","Zambia","Zimbabwe"]
variants={"Cape Verde":"Cabo Verde","Cote d'Ivoire":"Ivory Coast","Eswatini":"Swaziland","Democratic Republic of the Congo":"DR Congo / Zaire","Congo (Republic)":"Republic of the Congo","Tanzania":"United Republic of Tanzania"}
present=set(d.Country.dropna().unique())
wb=Workbook()
def sheet(name,header,rows,widths=None):
    ws=wb.create_sheet(name)
    ws.append(header)
    for c in ws[1]:
        c.font=Font(bold=True); c.alignment=Alignment(wrap_text=True,vertical='top')
    for r in rows: ws.append(r)
    if widths:
        for i,w in enumerate(widths,1): ws.column_dimensions[chr(64+i)].width=w
    return ws

# S2 country terms
rows=[]
for c in AU:
    rows.append([c, variants.get(c,''), 'yes' if c in present else 'no'])
sheet('S2_country_terms',['African Union state','Spelling variant or alias also searched','Genomes found'],rows,[30,34,16])

# S3 sector mapping
sm=pd.read_csv('sector_mapping.csv').fillna('(not given)')
sheet('S3_sector_mapping',['Host (as submitted)','Isolation source (as submitted)','Sector (derived)','Poultry niche (derived)','n'],sm.values.tolist(),[26,34,20,30,6])

# S4 sensitivity
s=json.load(open('sens.json'))
cats=['IncI2','IncX4','IncHI2','Mixed MOB-Recon bin','Other plasmid','Chromosomal or unresolved']
def srow(label,key):
    dd=s[key][0]; n=s[key][1]
    return [label,round(n,1)]+[f"{round(dd[c][0],1)} ({dd[c][1]}%)" for c in cats]
rows=[srow('All genomes','all'),
      srow('Excluding PRJNA1273464 (39 Algerian genomes)','excl_PRJNA1273464'),
      srow('Weighted by inverse BioProject size','bioproject_weighted'),
      srow('One genome per BioProject (modal)','one_per_bioproject_modal'),
      srow('One genome per SNP cluster','one_per_snp_cluster'),
      srow('Poultry sector, all','poultry_all'),
      srow('Poultry sector, excl. PRJNA1273464','poultry_excl')]
sheet('S4_sensitivity',['Analysis','Genomes (or weighted n)']+cats,rows,[38,18]+[16]*6)

# S5 evidence tiers
et=pd.read_csv('evidence_tiers.csv').fillna('')
sheet('S5_evidence_levels',['BioSample','Galaxy input','mcr gene','Replicon linked to mcr','Evidence source','mcr contig length (bp)','Evidence level'],et.values.tolist(),[18,18,10,30,14,18,26])

# S6 snp clusters
cl=json.load(open('../snp_clusters.json'))
rows=[]
for c in sorted(cl,key=lambda x:-x['n_ours']):
    rows.append([c['erd'],c['n_ours'],'; '.join(f"{k}:{v}" for k,v in c['countries'].items()),
                 '; '.join(f"{k}:{v}" for k,v in c['sectors'].items()),
                 '; '.join(f"ST{k}:{v}" for k,v in c['sts'].items()),
                 c['members'],'; '.join(f"{k}:{v}" for k,v in c['geos'].items())])
sheet('S6_snp_clusters',['PD SNP cluster','Our genomes in cluster','Our countries','Our sectors','Sequence types','Total PD members','Geography of all members'],rows,[18,16,22,28,28,16,34])

del wb['Sheet']
out='/mnt/user-data/outputs/Supplementary_Tables_S2-S6.xlsx'
wb.save(out)
print('saved',out)
print('countries matched:',len(present),sorted(present))
