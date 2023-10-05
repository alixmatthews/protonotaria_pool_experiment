#!/bin/bash

#SBATCH --job-name=03_SNPe_callvars_NCBIref_exp
#SBATCH --partition=comp01
#SBATCH --output=03_SNPe_callvars_NCBIref_exp_%j.txt
#SBATCH --error=03_SNPe_callvars_NCBIref_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:30:00

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
## Call variants on the mpileup file, joint-genotyping

bcftools call --threads 32 --skip-variants indels --multiallelic-caller --variants-only ${RESULTS_DIR}/bam/${SPP}_ALL_align_sort_dm_mq20_bam.mpileup > ${RESULTS_DIR}/vcf/${SPP}_ALL.vcf


##----------------------------------------------------------------------------------------
## Rename samples

## First, check the order using the following once this slurm is complete:
## bcftools query -l ${SPP}_ALL.vcf 

## Then create a file with the reheaders of the samples in the correct order (e.g., ${SPP}_reheader.txt).

## Then run this:
## bcftools reheader --samples ${SNP_DIR}/${SPP}_reheader.txt ${RESULTS_DIR}/vcf/${SPP}_ALL.vcf --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed.vcf

## Last, check the of the reheader'd samples and make sure it is the same as the original:
## bcftools query -l ${SPP}_ALL_renamed.vcf 


