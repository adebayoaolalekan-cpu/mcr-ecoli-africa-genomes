import pandas as pd, re, json, collections
import os
F=os.environ.get('MCR_WORKBOOK','Supplementary_Table_S1_per_genome_results.xlsx')
def load():
    d=pd.read_excel(F,sheet_name='Isolates')
    d=d[d['BioSample'].notna()].copy()
    g=pd.read_excel(F,sheet_name='Galaxy_results')
    d=d.merge(g,on='BioSample',how='left',suffixes=('','_g'))
    def cat(s):
        s=str(s)
        if 'mixed' in s: return 'Mixed MOB-Recon bin'
        if s.startswith('Chromosomal'): return 'Chromosomal or unresolved'
        if 'HI2' in s: return 'IncHI2'
        if s.startswith('IncI2'): return 'IncI2'
        if s.startswith('IncX4'): return 'IncX4'
        return 'Other plasmid'
    # QC exclusion: genomes failing the CheckM2 contamination ceiling stated in the
    # Methods are retained in the workbook but excluded from every analysis.
    inc=d.get('Included in analysis')
    if inc is not None:
        d=d[inc.astype(str).str.strip().str.lower()!='no'].copy()
    d['rep']=d['Plasmid replicon carrying mcr - to fill'].map(cat)
    d['evid']=d['Plasmid replicon carrying mcr - to fill'].map(lambda s:'MOB-Recon' if 'MOB-Recon' in str(s) else 'Same contig')
    d['gene']=d['mcr variant(s)'].map(lambda s:'mcr-9' if 'mcr-9' in s else ('mcr-10' if 'mcr-10' in s else 'mcr-1'))
    return d
