#!/bin/bash

#SBATCH --job-name=03_SNPg_filterminac_NCBIref_exp
#SBATCH --partition=comp01
#SBATCH --output=03_SNPg_filterminac_NCBIref_exp_%j.txt
#SBATCH --error=03_SNPg_filterminac_NCBIref_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:20:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate BCFTools # version 1.15.1

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW_EXP 

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_exp

## REF DIR, REF
REF_DIR=/scrfs/storage/amatthews/20210816_projects/20221012_finalgenome/IndexRef/PROW981R1_NCBI
REF=scaffolds_reduced_contigs_kept_NCBI-200000000

## SNP DIR 
SNP_DIR=${PROJECT_DIR}/04_SNP_20221013

## RESULTS DIR
RESULTS_DIR=${SNP_DIR}/RESULTS_${SPP}



##----------------------------------------------------------------------------------------
## Filter out SNPs that different from the reference but are equal to each other across samples

bcftools view --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1.vcf.gz --min-ac 1:minor ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30.vcf.gz








