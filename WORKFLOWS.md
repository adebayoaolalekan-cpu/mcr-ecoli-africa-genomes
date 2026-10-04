# Galaxy workflows

Both workflows were built and run on the public Galaxy server
https://usegalaxy.eu. The exported workflow files are included here as `.ga`
files in `workflows/`, ready to import into any Galaxy server. The workflow
identifiers on usegalaxy.eu are recorded below so the authors can share public
links at publication.

## 1. mcr Africa read to assembly
Workflow id: `fa9b3146518ac2e1` · file: `workflows/mcr_Africa_read_to_assembly.ga`

Steps: fasterq-dump (SRA Toolkit 3.1.1) → fastp 1.3.7 (Q20, min length 50,
minimum 153 Mb clean sequence retained) → Shovill 1.4.2 with SPAdes 3.15.5
(genome size 5.1 Mb, subsample to 100x, min contig 200 bp) → QUAST 5.3.0.

## 2. mcr Africa typing
Workflow id: `89e375fc444fb0b0` · file: `workflows/mcr_Africa_typing.ga`

Steps: AMRFinderPlus 4.2.7 (nucleotide mode, organism *Escherichia*, plus genes)
→ ABRicate 1.4.0 (NCBI, ResFinder, PlasmidFinder) → mlst 2.22.0 (scheme
ecoli_achtman_4 set explicitly under advanced settings) → MOB-Recon 3.1.9.

## Additional tools run outside the two workflows
- NCBI Datasets 18.33.1 (deposited GenBank assemblies, arm B)
- CheckM2 1.1.0, database 1.0.2 (completeness and contamination)
- ABRicate 1.4.0 against VFDB (species confirmation)
- ISEScan 1.7.3 (IS family context of *mcr-1*)
- tblastn, BLAST+ 2.16.0 (QseC/QseB downstream of *mcr-9*)

## Database build dates
All ABRicate bundled databases carried the build date 2026-04-03:
ncbi (8,232 sequences), resfinder (3,206), plasmidfinder (488), vfdb (4,592).
AMRFinderPlus used reference gene database 2026-05-15.1; the genotypes imported
from NCBI Pathogen Detection had been computed with 2026-03-24.1.

## Analysis histories
Three Galaxy histories hold the read-based arm, the deposited-assembly arm and
the combined typing. Public share links will be inserted at publication.
