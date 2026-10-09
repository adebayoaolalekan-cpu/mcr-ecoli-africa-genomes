# mcr-positive *Escherichia coli* from Africa in NCBI Pathogen Detection

Analysis code and data for the study of plasmid replicon types and sequence
types of *mcr*-positive *E. coli* genomes from Africa held in NCBI Pathogen
Detection (127 genomes retrieved from 16 countries, 2005–2025; 126 analysed
after quality control).

> **v2.1.0 applies the stated quality-control criterion.** One retrieved genome,
> GCA_015208565.1 (SAMN14219481, South African wastewater, *mcr-9*), returns 6.00% estimated
> CheckM2 contamination against the 5% ceiling set in the Methods. v2.0.0 retained it with a
> justification; v2.1.0 excludes it, as the criterion requires. The genome stays in
> `Supplementary_Table_S1_per_genome_results.xlsx` with its CheckM2 values and the reason for
> exclusion recorded, and the retrieval files still hold all 127, so both sets can be
> reconstructed. `base.py` filters on the `Included in analysis` column, so every script in this
> directory now reports the 126-genome analysed set.

> **v2.0.0 superseded v1.0.0.** The earlier release held 118 genomes from 13 countries. It was
> retrieved by filtering Pathogen Detection on a 45-name African country list at query time, and
> that list omitted nine of the 54 African Union member states, including South Africa, the Central
> African Republic and the Democratic Republic of the Congo. Nine genomes were missed. v2.0.0
> retrieves the global mcr-positive set for one pinned build and assigns country locally instead.
> All 118 earlier genomes are retained and none was dropped. See `CHANGELOG.md`.

Everything needed to regenerate the manuscript's tables, figure, sensitivity
analyses and evidence grading from the per-genome typing output is in this one
directory. The primary typing was run on the public Galaxy server at
https://usegalaxy.eu; the two workflows are included as `.ga` files and
described in `WORKFLOWS.md`.

## Quick start

```bash
pip install pandas openpyxl matplotlib
python sens.py      # sensitivity of the replicon distribution to uneven sampling
python tiers.py     # four-level evidence grading of each replicon assignment
python make_t2.py   # Table 2 (replicon by evidence level and mcr variant)
python make_t4.py   # Table 4 (replicon distribution under alternative sampling)
python desc.py      # distribution by year, host/source to sector mapping
python figure.py    # Figure 1 -> figure1.png
python supp.py      # Supplementary Tables S2-S6 -> .xlsx
python audit.py     # re-checks the manuscript's numbers against the data
```

All eight run from this directory with no arguments and were verified to do so
in a clean checkout. `audit.py` is the useful one for a reader who wants to
confirm the paper: it re-derives the headline counts and prints a pass/fail line
for each.

## Files

**Per-genome results**
- `Supplementary_Table_S1_per_genome_results.xlsx` — the master table. One row
  per genome with accessions, metadata, sector, *mcr* variant, sequence type,
  replicon, evidence level, CheckM2 completeness and contamination, and the
  IS30-family flag. Every script reads this. Override its location with the
  `MCR_WORKBOOK` environment variable if you move it.

**Typing output from Galaxy**
- `compactA.tsv`, `compactB.tsv` — AMRFinderPlus, PlasmidFinder, ABRicate and
  MLST calls per genome (arm A reads, arm B deposited assemblies).
- `mob_A.tsv`, `mob_B.tsv` — MOB-Recon plasmid bins (molecule, primary cluster,
  replicons, relaxase, bin size, *mcr* contig length).
- `galaxy_typing.tsv`, `galaxy_qc.tsv` — reconciliation inputs.
- `quast_A.tsv`, `quast_B.txt`, `fastp_A.tsv` — assembly QC.
- `checkm2_isescan_raw.txt` — raw CheckM2 completeness and contamination, and
  the ISEScan IS-family counts, per genome.
- `species_vfdb.txt` — VFDB species confirmation of the three genomes submitted
  under a non-*E.-coli* name. All three are *E. coli* and all three are retained.

**Retrieval** (added in v2)
- `retrieval/mcr_africa_rebuild.py` — the cohort retrieval pipeline. Downloads the full
  per-isolate metadata table for a pinned Pathogen Detection build, selects mcr-positive
  isolates by tokenised regex with qualified calls counted separately, and assigns country
  locally against a 54-state alias table. 70 boundary cases are asserted on every run, covering
  Niger against Nigeria, Guinea against Equatorial Guinea and Guinea-Bissau, South Africa against
  the phrase "southern Africa", the historical names, and the ISO code collisions that put
  Viet Nam in Namibia and `USA:GA` in Gabon.
- `retrieval/01_final_cohort.tsv` — the 127 retained isolates, with the assigned country, the
  alias that matched and a fallback flag.
- `retrieval/02_global_mcr_positive.tsv.gz` — all 12,560 mcr-positive isolates in the build
  before any geographic filtering.
- `retrieval/03_unresolved_location_residual.tsv.gz` — every mcr-positive isolate worldwide whose
  country could not be resolved, for manual inspection.
- `retrieval/04_country_alias_lookup.tsv` — the alias table as used, with the number of isolates
  each term matched, including the terms that matched none.
- `retrieval/06_reconciliation_added.tsv`, `retrieval/08_reconciliation_retained.txt` — the
  comparison against the v1 cohort.
- `retrieval/09_final_counts.tsv`, `retrieval/10_log.txt` — final counts, and the provenance log
  with a count at every filtering step.

**Typing of the nine genomes added in v2**
- `galaxy_D/` — AMRFinderPlus, MOB-Recon and MLST output for the nine genomes added by the
  corrected retrieval, run in the Galaxy history "mcr Africa D typing" with the same two
  workflows and the MLST scheme set explicitly, plus `nine_genomes_typed.tsv` summarising ST,
  mcr contig length, replicon and evidence level per genome.

**Computed outputs** (regenerated by the scripts above)
- `sens.json`, `evidence_tiers.csv`, `sector_mapping.csv`,
  `bioproject_summary.csv`, `snp_clusters.json`, `tables.json`.

**Galaxy workflows**
- `mcr_Africa_read_to_assembly.ga`, `mcr_Africa_typing.ga` — import into any
  Galaxy server.

**Provenance scripts** (kept for transparency, not expected to run standalone)
- `recon.py`, `rebuild.py` — rebuilt the results sheet from the Galaxy outputs
  and reconciled every filled cell against its source, finding 0 discrepancies.
- `build.py`, `analyze.py`, `add_checkm.py` — original build and annotation
  steps. These ran against live Galaxy outputs during the analysis.

## Tool and database versions

SRA Toolkit 3.1.1, fastp 1.3.7, Shovill 1.4.2 (SPAdes 3.15.5), QUAST 5.3.0,
NCBI Datasets 18.33.1, CheckM2 1.1.0 (database 1.0.2),
AMRFinderPlus 4.2.7 (reference gene database 2026-05-15.1; the genotypes
imported from Pathogen Detection used 2026-03-24.1), ABRicate 1.4.0 (bundled
databases all built 2026-04-03: ncbi 8,232 sequences, resfinder 3,206,
plasmidfinder 488, vfdb 4,592), mlst 2.22.0 (PubMLST *E. coli* Achtman scheme),
MOB-suite 3.1.9, ISEScan 1.7.3, BLAST+ 2.16.0.

## Licence

Code released under the MIT Licence. Genome data are public and remain under the
terms of their original submitters.
