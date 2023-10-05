#!/bin/bash

#SBATCH --job-name=04_plink_LDvari_PCA_99_NCBIref_exp
#SBATCH --partition=comp01
#SBATCH --output=04_plink_LDvari_PCA_99_NCBIref_exp_%j.txt
#SBATCH --error=04_plink_LDvari_PCA_99_NCBIref_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:10:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load plink

# version is PLINK v1.90b5.2

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

## PLINK DIR PER SPECIES
mkdir ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014
PLINK_DIR=${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014




##----------------------------------------------------------------------------------------
## Use PLINK to identify variants that are potentially in linkage disequilibrium and separate them. Output files will be: the prune.in (variants to keep) and prune.out (variants that are to be removed because potentially linked)

## Also make a .bim and .bed file for ADMIXTURE and .eigenval and .eigenvec files for plink-based PCA


## First identify variants that are potentially linked. Consider a window size of 100 SNPs, calculate LD between each pair of SNPs in the window, remove one of a pair of SNPs if the r2 value is greater than 0.001, 0.01, 0.1, 0.3, 0.4, then shift the window down 50 SNPs forward and repeat the procedure. 

## Do for all dp datasets and all maf datasets (0% missing) so can compare PCAs


cd ${PLINK_DIR}



## maf01, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99


## maf01, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99


## maf05, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99


## maf05, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99




## maf10, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99


## maf10, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 100 50 0.99 --out ${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99










## Then, make the .bim and .bed files for ADMIXTURE and files to create a PLINK-based PCA


## maf01, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf01_dp05_maxmiss100_ld99



## maf01, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf01_dp10_maxmiss100_ld99





## maf05, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf05_dp05_maxmiss100_ld99



## maf05, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99







## maf10, dp5

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf10_dp05_maxmiss100_ld99



## maf10, dp10

plink --vcf ../${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100.recode.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99.prune.in --make-bed --pca --out ${SPP}_ALL_renamed_q30_minac1_maf10_dp10_maxmiss100_ld99




