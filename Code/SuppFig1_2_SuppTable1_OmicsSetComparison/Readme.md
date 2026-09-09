# Omics Set Comparison (SuppFig1_2_SuppTable1_OmicsSetComparison)

This directory contains scripts for comparing compounds/sequences in this study against curated databases and large single-omics studies.

## Input data

### 1) Curated database comparison

- **MCDB (small molecules):**
  - `Data/DatabaseComparison/MCDB/milk_metabolites/milk_metabolites_summary_manual.csv`
- **MBPDB (peptidomics):**
  - concatenated from files in `Data/DatabaseComparison/MBPDB/`
- **BoMiProt (proteomics):**
  - protein sequences are stored in `Data/DatabaseComparison/BoMiProt/bomiprot_fasta.zip`

### 2) Largest single-omics study comparison

- **MCDB (small molecules):**
  - `Data/DatabaseComparison/MCDB/milk_metabolites/milk_metabolites_summary_manual.csv`
  - this file records source studies; the main study (PMID: 30994344) is treated as the largest single-omics experiment
- **Peptidomics:**
  - `Data/DatabaseComparison/SingleStudies/Peptidomics/PeptideSequence_.Ningetalcsv.csv`
  - corresponding downloaded supplementary source file is in the same directory
- **Proteomics:**
  - `Data/DatabaseComparison/SingleStudies/Proteomics/ProteinSequence_Chopraetal.csv`
  - corresponding downloaded supplementary source files are in the same directory

## Source code files

### 1) Curated database comparison

- `CompoundSetComparison_MCDB.R`
- `Peptidomics_Comparison_MBPDB.R`
- `Proteomics_Comparison_BoMiProt.R`

### 2) Largest single-omics study comparison

- `SmallMolecules_Comparison_MCDB_MainStudy.R`
- Peptidomics and proteomics comparisons are performed directly in R using simple `intersect()` and `setdiff()` operations.

## Output

- Main output: intersection counts between this study and published databases/studies.
- Venn diagrams and final summary tables are prepared manually from these outputs.
