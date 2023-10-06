# Feather mite genomics experiment
Alix Matthews


Separate READMEs for each step, but as a general overview, this project aims to compare the sequencing efficiency and population genetics of an experimental set of feather mite samples (1 mite vs. 5 mites vs. 20 mites in a single extraction - each from the same feather - of *Amerodectes protonotaria*).

## Files needed at various steps:
- These may or may not be in the subdirectories
`20210816_exp_filenames.txt`: all files, for the initial PP steps; used in `00_PP_20210903` steps
`20210816_exp_filenames_PROW953_only.txt`: only PROW953R1, which had a lot of contamination; used in `01_aTRAM_COI_20220901` step
`20210816_exp_filenames_SNPpipeline.txt`: all samples, just added an extra useful column, used in `04_SNP_20221013` steps
`20210816_exp_filenames_SNPpipeline_PROW953_only.txt`: only PROW953R1, which had a lot of contamination; used in `01_RemoveBacteria_20220901` steps

--

## Directories below and order of operations

## 00_PP_20210903

Initial pre-processing steps for the experimental samples (n=12). These are high in GC content, so will probably have to do several passes of pre-processing. This is a "standard" initial pre-processing.

## 01_aTRAM_COI_20220826

Making sure all samples are *A. protonotaria* ... but had issue with one (PROW953R1), so had to do a few more extra steps

## 01_RemoveBacteria_20220901

Removing bacteria reads from PROW953R1 and compare mapping rates

## 01_aTRAM_COI_20220901

Redoing aTRAM (COI) on all samples

## 03_IndexRef

Indexing the reference (NCBI deposited genome)

## 04_SNP_20221013

SNP pipeline for all samples

## 05_final-analysis-files

All files needed for the FINAL analyses. Lots of intermediate/unused files exist (in my other private repo, `alixmatthews/20210816_projects/20210816_exp`), but not in this cleaned up repo :)


