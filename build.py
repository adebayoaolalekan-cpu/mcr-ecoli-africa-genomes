import pandas as pd, re, collections
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

def load_compact(fn):
    rows=[l.rstrip('\n').split('\t') for l in open(fn) if l.strip()]
    return rows
rows=load_compact('compactB.tsv')+[r for r in load_compact('compactA.tsv') if r[0] in 'APRM']
amr=collections.defaultdict(list); pf=collections.defaultdict(list); ab=collections.defaultdict(lambda: collections.defaultdict(list)); mlst={}
for r in rows:
    if r[0]=='A' and r[3].startswith('mcr'): amr[r[1]].append(dict(sym=r[3],contig=r[2],meth=r[6],cov=float(r[7]),idt=float(r[8])))
    elif r[0]=='P': pf[r[1]].append((r[2],r[3],r[4],r[5]))
    elif r[0]=='R': ab[r[2]][r[1]].append(r[4])
    elif r[0]=='M': mlst[r[1]]=(r[2],r[3],' '.join(r[4:]))
mob={}
for fn in ['mob_B.tsv','mob_A.tsv']:
    for l in open(fn):
        f=l.rstrip('\n').split('\t'); mob[f[0]]=dict(mol=f[1],cl=f[2],reps=f[3],relax=f[4],binsize=int(f[5]),ncont=int(f[6]),contig_len=int(f[7]) if len(f)>7 else None)
# QC
qc={}
q=[l.rstrip('\n').split() for l in open('quast_B.txt')]
d={r[0]:r[1:] for r in q}
for i,n in enumerate(d['Assembly']):
    g='_'.join(n.split('_')[:2]); g=g[:-2]+'.'+g[-1] if False else n.split('_')[0]+'_'+n.split('_')[1]+'.'+n.split('_')[2]
    qc[g]=dict(contigs=int(d['Contigs'][i]),length=int(d['TotalLength'][i]),n50=int(d['N50'][i]),gc=float(d['GC'][i]))
q=[l.rstrip('\n').split('\t') for l in open('quast_A.tsv') if l.strip()]
d={r[0]:r[1:] for r in q}
for i,n in enumerate(d['Assembly']):
    qc[n]=dict(contigs=int(d['# contigs'][i]),length=int(d['Total length'][i]),n50=int(d['N50'][i]),gc=float(d['GC (%)'][i]))
fastp={l.split('\t')[0]:int(l.split('\t')[4]) for l in open('fastp_A.tsv')}

def norm(name):
    n=re.sub(r'_\d+(_Delta)?$','',name)
    return n
def mob_reps(s):
    reps=[x for x in s.split(',') if x and x!='-']
    inc=[x for x in dict.fromkeys(reps) if not x.startswith('rep_cluster')]
    return inc if inc else list(dict.fromkeys(reps))

orig=pd.read_excel('orig.xlsx',sheet_name='Isolates')
wb=load_workbook('orig.xlsx')
ws=wb['Isolates']
hdr={c.value:c.column for c in ws[1]}
det=[]
for ri in range(2,ws.max_row+1):
    bios=ws.cell(ri,hdr['BioSample']).value
    gca=ws.cell(ri,hdr['Assembly (GCA)']).value
    run=ws.cell(ri,hdr['SRA Illumina PE run(s)']).value
    cat=ws.cell(ri,hdr['Recommended Galaxy input']).value
    pdmcr=ws.cell(ri,hdr['mcr variant(s)']).value
    gid = run if cat=='SRA Illumina PE reads' else gca
    arm = 'A (reads, Shovill)' if cat=='SRA Illumina PE reads' else ('B (ENA assembly)' if gid=='GCA_988286365.1' else 'B (GenBank assembly)')
    notes=[]; st=''; rep=''
    calls=amr.get(gid,[])
    # de-duplicate identical calls on same contig
    seen=set(); uc=[]
    for c in calls:
        k=(c['sym'],c['contig'])
        if k in seen: continue
        seen.add(k); uc.append(c)
    calls=uc
    exclude = (gid=='SRR9099604')
    s=mlst.get(gid)
    if s:
        if s[1]=='-' or s[0]!='ecoli_achtman_4':
            if s[0]!='ecoli_achtman_4':
                notes.append(f'MLST returned the {s[0]} scheme despite ecoli_achtman_4 being selected; rerun needed')
                st='Untypeable'
            elif '~' in s[2] or '?' in s[2]:
                notes.append('Novel or partial MLST allele ('+s[2]+')'); st='Untypeable'
            else:
                notes.append('All seven alleles known but no ST assigned; likely novel ST ('+s[2]+')'); st='Untypeable'
        else:
            st=int(s[1])
    links=[]
    for c in calls:
        reps=list(dict.fromkeys(norm(p[1]) for p in pf.get(gid,[]) if p[0]==c['contig']))
        links.append((c,reps))
    contig_reps=' / '.join(sorted(set(r for _,rr in links for r in rr)))
    mobinfo=mob.get(gid)
    if gid=='GCA_007993675.2':
        rep='IncI2 (second copy on chromosome-size contig, see notes)'
        notes.append('Two mcr-1.1 copies in a closed genome; one on the IncI2 plasmid, one on the 4.96 Mb contig CP042470.1, which also carries IncHI2, IncFIB(AP001918), IncI2 and IncX1 replicons. Check whether the second copy is chromosomal or a co-integrate')
    elif contig_reps:
        rep=contig_reps
    elif mobinfo:
        r=mob_reps(mobinfo['reps'])
        if mobinfo['mol']=='chromosome':
            rep='Chromosomal or unresolved (MOB-Recon)'
        elif not r:
            rep='Plasmid, no replicon detected (MOB-Recon)'
        else:
            mixed = len(set(x.replace('IncHI2A','IncHI2') for x in r if x!='IncQ1'))>1 and not all(x.startswith('IncF') for x in r) and not (any(x.startswith('IncHI1') for x in r) and all(x.startswith('IncHI1') or x.startswith('IncF') for x in r))
            rep='/'.join(r)+' (MOB-Recon bin'+(', mixed replicons' if mixed else '')+')'
            if mixed: notes.append('MOB-Recon bin holds more than one plasmid family; mcr plasmid type not resolved with short reads')
        cl=mobinfo.get('contig_len')
        if cl and cl<5000: notes.append(f'mcr contig only {cl:,} bp; bin assignment is tentative')
    elif calls:
        rep='No replicon on mcr contig'
    for c,_ in links:
        if not c['meth'].startswith('ALLELE'):
            notes.append(f"{c['sym']} call is {c['meth']} ({c['cov']:.1f}% coverage, {c['idt']:.2f}% identity); check before counting")
    galmcr=';'.join(dict.fromkeys(c['sym'] for c in calls))
    pdset=set(str(pdmcr).replace(' ','').split(';')) if pdmcr else set()
    galset=set(c['sym'] for c in calls)
    if galset!=pdset and not exclude:
        notes.append(f'Galaxy AMRFinderPlus call ({galmcr or "none"}) differs from PD ({pdmcr})')
    if gid=='GCA_988286365.1': notes.append('Assembly not in NCBI Datasets; downloaded from ENA (10 contigs, 5.36 Mb)')
    if exclude:
        notes=[f'Excluded at read QC: {fastp[gid]/1e6:.1f} Mb clean bases (~{fastp[gid]/5.1e6:.1f}x, below 30x). Shovill assembly gave ST{st}; GenBank assembly GCA_006348965.1 could be typed instead']
        st=''; rep=''
    if st!='': ws.cell(ri,hdr['ST (Achtman) - to fill']).value=st
    if rep: ws.cell(ri,hdr['Plasmid replicon carrying mcr - to fill']).value=rep
    if notes: ws.cell(ri,hdr['Notes - to fill']).value='. '.join(notes)+'.'
    q=qc.get(gid,{})
    m=mobinfo or {}
    det.append([bios,gid,arm,fastp.get(gid),q.get('contigs'),q.get('length'),q.get('n50'),q.get('gc'),
                'Excluded (depth)' if exclude else ('Pass' if q else 'Not run (ENA upload)'),
                pdmcr, galmcr, '; '.join(c['meth'] for c in calls), '; '.join(f"{c['cov']:.2f}" for c in calls), '; '.join(f"{c['idt']:.2f}" for c in calls),
                '; '.join(dict.fromkeys(c['contig'] for c in calls)),
                ';'.join(dict.fromkeys(ab[gid].get('ncbi',[]))) or 'none', ';'.join(dict.fromkeys(ab[gid].get('resfinder',[]))) or 'none',
                ' / '.join(dict.fromkeys(p[1] for c,_ in links for p in pf.get(gid,[]) if p[0]==c['contig'])) or 'none',
                m.get('mol',''), m.get('reps',''), m.get('relax',''), m.get('binsize') or None,
                (mlst[gid][1] if gid in mlst else ''), (mlst[gid][2] if gid in mlst else ''),
                ', '.join(sorted(set(p[1] for p in pf.get(gid,[])))) ])
cols=['BioSample','Galaxy input ID','Galaxy arm','Clean bases after fastp','Est. depth (x, 5.1 Mb genome)','Contigs (QUAST)','Total length (bp)','N50 (bp)','GC (%)','QC status',
      'mcr (PD)','mcr (Galaxy AMRFinderPlus)','AMRFinderPlus method','% coverage','% identity','mcr contig',
      'ABRicate NCBI mcr','ABRicate ResFinder mcr','PlasmidFinder replicon(s) on mcr contig','MOB-Recon molecule','MOB-Recon bin replicon(s)','MOB-Recon relaxase(s)','MOB-Recon bin size (bp)',
      'MLST ST (raw)','MLST alleles','All PlasmidFinder replicons in genome']
ns=wb.create_sheet('Galaxy_results')
hfont=Font(name='Arial',bold=True,color='FFFFFFFF'); hfill=PatternFill('solid',fgColor='FF1F4E78'); bfont=Font(name='Arial',size=10)
for j,c in enumerate(cols,1):
    x=ns.cell(1,j,c); x.font=hfont; x.fill=hfill; x.alignment=Alignment(wrap_text=True,vertical='top')
for i,r in enumerate(det,2):
    for j,v in enumerate(r,1):
        if j==5: v=f'=IF(D{i}="","",D{i}/5100000)' 
        x=ns.cell(i,j,v); x.font=bfont
    ns.cell(i,5).number_format='0.0'
    for j in (4,6,7,8,23): ns.cell(i,j).number_format='#,##0'
for j in range(1,len(cols)+1): ns.column_dimensions[get_column_letter(j)].width=16
for j in (12,19,20,25,26): ns.column_dimensions[get_column_letter(j)].width=30
ns.freeze_panes='C2'; ns.auto_filter.ref=f'A1:{get_column_letter(len(cols))}{len(det)+1}'
ns.row_dimensions[1].height=45
for c in ('Z','AA','AB'):
    for ri in range(2,ws.max_row+1): ws[f'{c}{ri}'].alignment=Alignment(wrap_text=False)
ws.column_dimensions['AA'].width=40; ws.column_dimensions['AB'].width=60

# Summary additions
s=wb['Summary']
sectors=['Poultry','Human','Other animal','Food (non-poultry)','Environment','Unspecified']
start=39
s.cell(start,1,'Genomes by sequence type (Achtman) and sector, STs seen in two or more genomes').font=Font(name='Arial',bold=True)
s.cell(start+1,1,'ST').font=Font(name='Arial',bold=True)
for j,sec in enumerate(sectors,2): s.cell(start+1,j,sec).font=Font(name='Arial',bold=True)
s.cell(start+1,8,'Total').font=Font(name='Arial',bold=True)
stv=[ws.cell(ri,hdr['ST (Achtman) - to fill']).value for ri in range(2,ws.max_row+1)]
cnt=collections.Counter(v for v in stv if isinstance(v,int))
sts=[k for k,v in sorted(cnt.items(),key=lambda x:(-x[1],x[0])) if v>=2]
r=start+2
for stn in sts+['Untypeable']:
    s.cell(r,1,stn)
    for j in range(2,8):
        col=get_column_letter(j)
        s.cell(r,j,f'=COUNTIFS(Isolates!$Z$2:$Z$119,$A{r},Isolates!$O$2:$O$119,{col}${start+1})')
    s.cell(r,8,f'=SUM(B{r}:G{r})'); r+=1
s.cell(r,1,'Other STs (single genome)')
for j in range(2,8):
    col=get_column_letter(j)
    s.cell(r,j,f'=COUNTIFS(Isolates!$O$2:$O$119,{col}${start+1},Isolates!$Z$2:$Z$119,"<>")-SUM({col}{start+2}:{col}{r-1})')
s.cell(r,8,f'=SUM(B{r}:G{r})'); r+=1
s.cell(r,1,'Total typed or untypeable')
for j in range(2,9):
    col=get_column_letter(j); s.cell(r,j,f'=SUM({col}{start+2}:{col}{r-1})')
r+=2
s.cell(r,1,'Genomes by plasmid replicon linked to mcr and sector (wildcard match on the Isolates replicon column; mixed MOB-Recon bins count under each family)').font=Font(name='Arial',bold=True)
r+=1; h=r
s.cell(h,1,'Replicon').font=Font(name='Arial',bold=True)
for j,sec in enumerate(sectors,2): s.cell(h,j,sec).font=Font(name='Arial',bold=True)
s.cell(h,8,'Total').font=Font(name='Arial',bold=True)
r+=1
for lab,pat in [('IncI2','*IncI2*'),('IncX4','*IncX4*'),('IncHI2','*IncHI2*'),('IncHI1','*IncHI1*'),('IncF (FII/FIB/FIA)','*IncF*'),('Other rep cluster (no Inc type)','rep_cluster*'),('Plasmid, no replicon detected','Plasmid, no replicon*'),('Chromosomal or unresolved','Chromosomal*')]:
    s.cell(r,1,lab)
    for j in range(2,8):
        col=get_column_letter(j)
        s.cell(r,j,f'=COUNTIFS(Isolates!$AA$2:$AA$119,"{pat}",Isolates!$O$2:$O$119,{col}${h})')
    s.cell(r,8,f'=SUM(B{r}:G{r})'); r+=1
s.cell(r+1,1,'Note: ST and replicon columns were filled from the Galaxy run (usegalaxy.eu, histories "mcr Africa A reads", "mcr Africa B assemblies", "mcr Africa C typing"). Per-genome detail is in the Galaxy_results sheet.')
for row in s.iter_rows(min_row=start,max_row=r+1):
    for c in row:
        if c.font is None or not c.font.bold: c.font=Font(name='Arial',bold=bool(c.font and c.font.bold))
# README note
rd=wb['README']
rr=rd.max_row+2
for t in ['Galaxy run (October 2026)',
          'Tools on usegalaxy.eu: fasterq-dump 3.1.1, fastp 1.3.7 (Q20, min length 50), Shovill 1.4.2 (SPAdes, 5.1M, depth 100, min contig 200), QUAST 5.3.0, NCBI Datasets 18.33.1, AMRFinderPlus 4.2.7 (database V4.2-2026-05-15.1, organism Escherichia, plus genes on), ABRicate 1.4.0 (NCBI and ResFinder 80/80, PlasmidFinder 95/60), MLST 2.22.0 (ecoli_achtman_4), MOB-Recon 3.1.9.',
          'Galaxy_results: QC, mcr call, mcr contig, replicon linkage and MLST for each of the 118 genomes. Yellow columns in Isolates were filled from this sheet.']:
    rd.cell(rr,1,t).font=Font(name='Arial',bold=(t.startswith('Galaxy run'))); rr+=1
wb.save('African_mcr_Ecoli_NCBI_dataset_Galaxy_results.xlsx')
print('done', len(det), sts)
