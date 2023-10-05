#!/bin/bash

#SBATCH --job-name=03_SNPj_filtermissing_exp
#SBATCH --partition=pcon06
#SBATCH --output=03_SNPj_filtermissing_exp_%j.txt
#SBATCH --error=03_SNPj_filtermissing_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load vcftools/0.1.15

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW_EXP 

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_exp

## REF DIR, REF
REF_DIR=$PROJECT_DIR/03_IndexRef/PROW981R1_ref_full_no_microbes
REF=scaffolds_reduced_contigs_kept

## SNP DIR 
SNP_DIR=$PROJECT_DIR/04_SNP_20220723

## RESULTS DIR
RESULTS_DIR=$SNP_DIR/RESULTS_$SPP



##----------------------------------------------------------------------------------------
## Filter out SNP sites that have any missing data (genotyping was successful in 100% of individuals...) using maf05 and maf10


# DP threshold of 5

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp05.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp05.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100




# DP threshold of 10

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp10.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp10.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100




# DP threshold of 15

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp15.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp15_maxmiss100

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp15.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp15_maxmiss100




# DP threshold of 20

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp20.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp20_maxmiss100

vcftools --gzvcf ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp20.vcf.gz --max-missing 1.0 --recode --recode-INFO-all --out ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp20_maxmiss100



