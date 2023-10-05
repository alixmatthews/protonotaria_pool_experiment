#!/bin/bash

#SBATCH --job-name=03_SNPi_filterdp_exp
#SBATCH --partition=comp01
#SBATCH --output=03_SNPi_filterdp_exp_%j.txt
#SBATCH --error=03_SNPi_filterdp_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

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
REF_DIR=$PROJECT_DIR/03_IndexRef/PROW981R1_ref_full_no_microbes
REF=scaffolds_reduced_contigs_kept

## SNP DIR 
SNP_DIR=$PROJECT_DIR/04_SNP_20220723

## RESULTS DIR
RESULTS_DIR=$SNP_DIR/RESULTS_$SPP



##----------------------------------------------------------------------------------------
## Filter by seq depth, several different thresholds (5, 10, 15, 20) with maf05 and maf10

## So apparently there is a bug for the +setGT command: https://github.com/samtools/bcftools/issues/1607. It says it was fixed in the newest version (the one I'm using), but maybe I am just thinking about the --include vs. --exclude backwards.

## The solution is use the filtering expression which kind of the opposite of what you would think. So here I am filtering by DP. I want to change the GT for a sample to ./. if DP is less than 5/10/15/20. Intuitively, I would put {--exclude 'FMT/DP<5}; however, that is not the case (by checking the resulting .vcf files, which did not do what I wanted at all - it was changing everything backwards or in some weird pattern). Instead, here, we are using {--include 'FORMAT/DP<5}. The resulting .vcf files are behaving as expected.

## maf05
bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp05.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<5' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp10.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<10' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp15.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<15' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf05_dp20.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<20' 

## maf10
bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp05.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<5' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp10.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<10' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp15.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<15' 

bcftools +setGT ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10.vcf.gz --output-type z --output ${RESULTS_DIR}/vcf/${SPP}_ALL_renamed_q30_minac1_maf10_dp20.vcf.gz -- --target-gt q --new-gt . --include 'FORMAT/DP<20' 






