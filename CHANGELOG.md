# Changelog

## v2.1.0 — 2026-10-08

One genome is excluded, and the analysed collection becomes **126 of the 127 retrieved**.

### Why

GCA_015208565.1 (SAMN14219481), the South African wastewater isolate carrying mcr-9, returns 6.00%
estimated CheckM2 contamination. The Methods set a 5% ceiling. v2.0.0 retained the genome and
argued the case in the Limitations; v2.1.0 applies the criterion as written and leaves it out. It
is the only genome in the collection to exceed the ceiling, and the highest contamination among the
126 analysed genomes is 3.64%.

The genome is not deleted. It stays in `Supplementary_Table_S1_per_genome_results.xlsx` with its
CheckM2 values, and two new columns record the decision:

```
Included in analysis    Yes / No
Reason for exclusion    CheckM2 contamination 6.00%, above the 5% ceiling stated in the Methods
```

`base.py` filters on that column, so every script in this directory reports the analysed set while
the retrieval files under `retrieval/` continue to hold all 127. Either set can be reconstructed
from what is here.

### What moves

| | v2.0.0 | v2.1.0 |
|---|---|---|
| Genomes analysed | 127 | 126 |
| Countries | 16 | 16 |
| mcr-1 | 116 | 116 |
| mcr-9 | 10 | 9 |
| mcr-10.1 | 1 | 1 |
| Environment sector | 4 | 3 |
| IncI2 | 61 (48.0%) | 61 (48.4%) |
| IncX4 | 30 (23.6%) | 30 (23.8%) |
| IncHI2 | 21 (16.5%) | 21 (16.7%) |
| Same-contig assignments | 83 | 83 |
| MOB-Recon assignments | 44 | 43 |
| Evidence levels 1/2/3/4 | 1/82/26/18 | 1/82/26/17 |
| Genomes with a numbered ST | 122 | 121 |
| Distinct numbered STs | 50 | 50 |
| BioProjects | 43 | 42 |
| Poultry share | 48.8% | 49.2% |

The substantive consequence is for mcr-9. The excluded genome was the only environmental mcr-9 in
the collection, so the claim that mcr-9 had been found in the environment is withdrawn. What
remains is one food genome, a South African cucumber carrying mcr-9.1. The association between
mcr-9 and IncHI2 still rests on six MOB-Recon assignments, now out of nine mcr-9 genomes rather
than ten.

### Files changed

- `base.py` — filters on `Included in analysis`.
- `audit.py` — assertions updated, and three new checks: 127 retrieved rows, exactly one
  QC-excluded row, and no analysed genome above 5% contamination. All 19 checks pass.
- `Supplementary_Table_S1_per_genome_results.xlsx` — two columns added to `Isolates`; the `Summary`
  sheet now counts the analysed set and carries a retrieval-and-quality-control block.
- `Supplementary_Tables_S2-S6.xlsx`, `tables.json`, `sens.json`, `figure1.png`,
  `evidence_tiers.csv`, `bioproject_summary.csv` — regenerated.
- `README.md` — cohort description updated.

### Corrections carried in the same release

Three numbers in the v2.0.0 manuscript draft were wrong independently of the exclusion and are
corrected here and in the paper: the count of genomes sharing a SNP cluster with another African
genome (47, not 53), the CTX-M count among human isolates (14 of 34, not 15), and the statement
that every unresolved mcr contig was shorter than 3.5 kb, which does not hold for a Togolese
mcr-9.1 on a 7.9 kb contig.

## v2.0.0 — 2026-10-08

The cohort grows from 118 genomes in 13 countries to **127 genomes in 16 countries**. All 118 genomes
from v1.0.0 are retained. None was dropped, so no record in the earlier set had been revised or
suppressed by NCBI in the interval.

### Why the cohort changed

v1.0.0 was built by querying NCBI Pathogen Detection on 26 September 2026 and filtering on a list of
African country names **at query time**. That list carried 45 names for the 54 African Union member
states. Nine member states were missing from it, among them South Africa, the Central African Republic
and the Democratic Republic of the Congo, and nine genomes were lost with them.

The missed genomes were not awkward edge cases. Every one carried a plain canonical location string,
several reading simply `South Africa` or `Central African Republic`, and two had been deposited in 2018
and 2019, years before the original query. Nothing in the v1.0.0 output signalled that anything was
absent.

v2.0.0 was rebuilt the other way round: the complete global set of mcr-positive isolates was downloaded
for one pinned Pathogen Detection build and country was assigned **locally**, against an auditable alias
table that reports how many isolates each term matched, including the terms that matched none.

### Source build

```
organism group   E.coli and Shigella
build            PDG000000004.6345   (published 2026-10-06 12:12 US Eastern)
isolates         587,057
metadata SHA-256 30b72681a8ae923026c867f53214360d286dcc823a6474fbb6d22608cba56c35
```

12,560 isolates in that build are mcr-positive worldwide before any geographic filtering. 127 resolve to
an African Union member state.

### Two matching defects found and fixed during the rebuild

Both would have corrupted a country assignment made by naive string matching, and both are now asserted
in the test set that runs on every execution of `retrieval/mcr_africa_rebuild.py`.

- **ISO alpha-3 matched as a token.** `Viet Nam` contains the whole token `nam`, Namibia's alpha-3 code,
  which assigned 178 Vietnamese isolates to Namibia. Alpha-3 codes are now accepted only when a candidate
  string is exactly the code.
- **ISO alpha-2 against US state abbreviations.** `USA:GA` resolved to Gabon, `USA:SC` to the Seychelles,
  `USA:SD` to Sudan, `USA:MA` to Morocco, `USA:NE` to Niger and `USA:TN` to Tunisia. Two-letter codes are
  no longer matched at all.

A first pass reported 330 isolates across 20 countries; 203 of those were these two defects.

### The nine added genomes

| BioSample | Assembly | Country | Source | mcr | ST | Replicon | Level |
|---|---|---|---|---|---|---|---|
| SAMEA120683848 | GCA_977857845.1 | Central African Republic | patient | mcr-1.1 | 226 | IncX4 | 2 |
| SAMEA120683850 | GCA_977859995.1 | Central African Republic | patient | mcr-1.1 | 226 | IncX4 | 2 |
| SAMEA5988742 | GCA_902703155.1 | South Africa | patient | mcr-1.1 | 602 | IncI2 | 2 |
| SAMN09303237 | GCA_003795445.1 | South Africa | pig, nares | mcr-1.1 | 9440 | IncHI2 (bin) | 4 |
| SAMN14219481 | GCA_015208565.1 | South Africa | wastewater | mcr-9 | 69 | unresolved | 4 |
| SAMN28613555 | GCA_025336805.1 | DR Congo | *Papio cynocephalus*, faeces | mcr-1.1 | 69 | unresolved | 4 |
| SAMN28613556 | GCA_025336825.1 | DR Congo | *Papio cynocephalus*, faeces | mcr-1.1 | 101 | IncI2 | 2 |
| SAMN36943049 | GCA_034465125.1 | South Africa | patient | mcr-1.1 ×2 | 226 | unresolved | 4 |
| SAMN49772998 | GCA_053709775.1 | South Africa | cucumber | mcr-9.1 | 297 | unresolved | 4 |

They were typed in the Galaxy history `mcr Africa D typing` using the same two workflows, with the MLST
scheme `ecoli_achtman_4` set explicitly.

### Headline numbers

| | v1.0.0 | v2.0.0 |
|---|---|---|
| Genomes | 118 | 127 |
| Countries | 13 | 16 |
| mcr-1 | 109 | 116 |
| mcr-9 | 8 | 10 |
| mcr-10.1 | 1 | 1 |
| IncI2 | 59 (50.0%) | 61 (48.0%) |
| IncX4 | 28 (23.7%) | 30 (23.6%) |
| IncHI2 | 20 (16.9%) | 21 (16.5%) |
| Same-contig assignments | 79 | 83 |
| MOB-Recon assignments | 39 | 44 |
| Evidence levels 1/2/3/4 | 1/78/26/13 | 1/82/26/18 |
| Numbered sequence types | 48 | 50 |
| Poultry share | 52.5% | 48.8% |

Two results that did **not** move are worth stating, because the opposite might have been expected.
Neither added mcr-9 genome could be placed on a plasmid: both sit on short contigs in chromosome bins.
The association between mcr-9 and IncHI2 therefore still rests on the same six MOB-Recon assignments,
now out of 10 mcr-9 genomes rather than 8, so adding data weakened the proportion supporting it. And
poultry, while still the largest single sector at 62 genomes, is no longer a majority of the collection.

### Files added

```
retrieval/mcr_africa_rebuild.py     the retrieval pipeline, with 70 asserted boundary cases
retrieval/01_final_cohort.tsv       the 127 retained isolates
retrieval/02_global_mcr_positive.tsv.gz   all 12,560 mcr-positive isolates worldwide, pre-geography
retrieval/03_unresolved_location_residual.tsv.gz   every isolate whose country could not be resolved
retrieval/04_country_alias_lookup.tsv     the alias table, with per-term match counts
retrieval/06_reconciliation_added.tsv     the nine, with why the earlier search missed each
retrieval/08_reconciliation_retained.txt  the 118 retained accessions
retrieval/09_final_counts.tsv       counts by country, mcr variant and sector
retrieval/10_log.txt                provenance, every filtering step with its count
galaxy_D/                           typing output for the nine
```

### Files changed

- `Supplementary_Table_S1_per_genome_results.xlsx` — nine rows added to `Isolates` and `Galaxy_results`,
  one to `SRA_runs_for_Galaxy`, and the `Summary` sheet rebuilt for 16 countries with ranges extended to
  row 128.
- `galaxy_typing.tsv` — nine arm-B rows added, so `tiers.py` and `make_t2.py` grade the new genomes.
- `audit.py` — assertions updated to the v2 numbers. All 16 checks pass.
- `README.md` — cohort description updated.

### Not yet done

ISEScan and SNP cluster membership have not been extracted for the nine added genomes, so the
IS30-family count describes 109 mcr-1 genomes and the SNP cluster analysis describes 118. The
one-genome-per-SNP-cluster row of the sensitivity table treats each of the nine as its own cluster,
which is the conservative assumption. QUAST and CheckM2 were not rerun for the nine, so the assembly
quality statistics also describe the first 118.

## v1.0.0 — 2026-10-04

First release. 118 genomes, 13 countries, 2005 to 2025.
