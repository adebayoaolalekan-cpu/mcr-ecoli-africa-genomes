# mcr-positive *Escherichia coli* from Africa in NCBI Pathogen Detection

Analysis code and intermediate data for the study of plasmid replicon types and
sequence types of *mcr*-positive *E. coli* genomes from Africa held in NCBI
Pathogen Detection (118 genomes, 13 countries, 2005–2025).

This repository lets a reader reproduce the summary tables, the figure, the
sensitivity analyses, the evidence-level grading and the SNP-cluster analysis
from the per-genome typing output. The primary typing was run on the public
Galaxy server at https://usegalaxy.eu; the two workflows are described in
`WORKFLOWS.md`.

## Layout

```
scripts/   Python scripts that build the tables, figure and analyses
data/      per-genome typing output and intermediate tables
```

### scripts
- `base.py` – loads the per-genome results workbook and derives replicon
  categories, evidence source and *mcr* gene group.
- `recon.py`, `rebuild.py` – rebuild the results sheet from the Galaxy outputs
  and reconcile every filled cell against its source (0 discrepancies).
- `sens.py` – sensitivity of the replicon distribution to uneven sampling
  (exclude the large Algerian submission, inverse-BioProject weighting, one
  genome per BioProject, one genome per SNP cluster).
- `tiers.py` – four-level evidence grading of each replicon assignment.
- `desc.py` – distribution by year and the host/source to sector mapping.
- `figure.py` – Figure 1 (replicon by sector and by country, MOB-derived
  segments hatched).
- `make_t4.py`, `supp.py` – build Table 4 and Supplementary Tables S2–S6.
- `analyze.py`, `build.py`, `pd_parse.py` – helper scripts used during the
  original analysis and the Pathogen Detection metadata parse.

### data
- `compactA.tsv`, `compactB.tsv` – AMRFinderPlus, PlasmidFinder, ABRicate and
  MLST calls per genome (arm A reads, arm B deposited assemblies).
- `mob_A.tsv`, `mob_B.tsv` – MOB-Recon plasmid bins (molecule, primary cluster,
  replicons, relaxase, bin size, *mcr* contig length).
- `galaxy_typing.tsv`, `galaxy_qc.tsv` – reconciliation inputs.
- `quast_A.tsv`, `quast_B.txt`, `fastp_A.tsv` – assembly QC.
- `snp_clusters.json` – Pathogen Detection single-linkage SNP-cluster membership.
- `sens.json`, `evidence_tiers.csv`, `sector_mapping.csv`,
  `bioproject_summary.csv` – computed outputs.
- `species_vfdb.txt` – VFDB species confirmation of the three genomes submitted
  under a non-*E.-coli* name.

## Reproducing

```bash
pip install pandas openpyxl matplotlib
cd scripts
python sens.py      # sensitivity analyses -> sens.json
python tiers.py     # evidence levels -> evidence_tiers.csv
python figure.py    # Figure 1
python supp.py      # Supplementary Tables S2-S6
```

`base.py` expects the per-genome results workbook
(`African_mcr_Ecoli_NCBI_dataset_Galaxy_results.xlsx`, Supplementary Table S1)
in the path set at the top of the file.

## Tool and database versions

SRA Toolkit 3.1.1, fastp 1.3.7, Shovill 1.4.2 (SPAdes 3.15.5), QUAST 5.3.0,
NCBI Datasets 18.33.1, CheckM2 1.1.0 (database 1.0.2),
AMRFinderPlus 4.2.7 (reference gene database 2026-05-15.1; Pathogen Detection
genotypes used 2026-03-24.1), ABRicate 1.4.0 (bundled NCBI, ResFinder and
PlasmidFinder databases), mlst 2.22.0 (PubMLST *E. coli* Achtman scheme),
MOB-suite 3.1.9, ISEScan 1.7.3, BLAST+ 2.16.0.

## Licence

Code released under the MIT Licence. Genome data are public and remain under the
terms of their original submitters.
