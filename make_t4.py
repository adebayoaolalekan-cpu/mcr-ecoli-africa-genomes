import json
s=json.load(open('sens.json'))
def row(label,key):
    d=s[key][0]; n=s[key][1]
    g=lambda c: f"{int(round(d[c][0]))} ({d[c][1]})"
    return [label, str(int(round(n))), g('IncI2'), g('IncX4'), g('IncHI2'), g('Other/mixed/unresolved') if 'Other/mixed/unresolved' in d else str(int(round(d['Mixed MOB-Recon bin'][0]+d['Other plasmid'][0]+d['Chromosomal or unresolved'][0])))]
T4=[['Analysis','Genomes','IncI2 n (%)','IncX4 n (%)','IncHI2 n (%)','Other n']]
T4.append(row('All genomes','all'))
T4.append(row('Excluding the 39-genome Algerian submission (PRJNA1273464)','excl_PRJNA1273464'))
T4.append(row('Weighted by inverse BioProject size','bioproject_weighted'))
T4.append(row('One genome per BioProject (modal replicon)','one_per_bioproject_modal'))
T4.append(row('One genome per SNP cluster','one_per_snp_cluster'))
T=json.load(open('tables.json'))
T['T4']=T4
json.dump(T,open('tables.json','w'),indent=0)
for r in T4: print(r)
