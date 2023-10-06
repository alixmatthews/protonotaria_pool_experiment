#!/bin/bash

#SBATCH --job-name=03_SNPk_calculateSNPdepth-sambambaonly_NCBIref_exp
#SBATCH --partition=comp01
#SBATCH --output=03_SNPk_calculateSNPdepth-sambambaonly_NCBIref_exp_%j.txt
#SBATCH --error=03_SNPk_calculateSNPdepth-sambambaonly_NCBIref_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load samtools/1.15.1
# module load bedtools2/2.25.0
module load sambamba/0.8.2


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
## Compare samtools depth and sambamba coverage (removed bedtools)




################################### samtools #############################################

#### no missing data (maxmiss100) - ld001

# samtools depth -b ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.prune.in.samtools.bed -f ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --min-BQ 30 -H -o ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001_samtools_depth.txt


#### no missing data (maxmiss100) - ld99

# samtools depth -b ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99.prune.in.samtools.bed -f ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --min-BQ 30 -H -o ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99_samtools_depth.txt


#### includes missing data - ld001

# samtools depth -b ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld001.prune.in.samtools.bed -f ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --min-BQ 30 -H -o ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld001_samtools_depth.txt


#### includes missing data - ld99

# samtools depth -b ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld99.prune.in.samtools.bed -f ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --min-BQ 30 -H -o ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld99_samtools_depth.txt








################################### sambamba #############################################

BAMFILE_01=${RESULTS_DIR}/bam/PROW_953_R1_TGACAACC-CTGTTGAC_align_sort_dm_mq20.bam
BAMFILE_02=${RESULTS_DIR}/bam/PROW_953_R5_CAGTGCTT-TGAGGTGT_align_sort_dm_mq20.bam
BAMFILE_03=${RESULTS_DIR}/bam/PROW_954_R20_ACAGCAAG-TTGCAGAC_align_sort_dm_mq20.bam
BAMFILE_04=${RESULTS_DIR}/bam/PROW_981_R1_TGCTGTGA-GAATCGTG_align_sort_dm_mq20.bam
BAMFILE_05=${RESULTS_DIR}/bam/PROW_984_R1_TAGTGCCA-TCTTGACG_align_sort_dm_mq20.bam
BAMFILE_06=${RESULTS_DIR}/bam/PROW_953_R20_TCACTCGA-GTGCCATA_align_sort_dm_mq20.bam
BAMFILE_07=${RESULTS_DIR}/bam/PROW_954_R1_CTGACTAC-ACATTGCG_align_sort_dm_mq20.bam
BAMFILE_08=${RESULTS_DIR}/bam/PROW_984_R5_TCGTGCAT-TGGTACAG_align_sort_dm_mq20.bam
BAMFILE_09=${RESULTS_DIR}/bam/PROW_981_R5_CAACACAG-AGGCTTCT_align_sort_dm_mq20.bam
BAMFILE_10=${RESULTS_DIR}/bam/PROW_954_R5_GTGATCCA-CGTTGAGT_align_sort_dm_mq20.bam
BAMFILE_11=${RESULTS_DIR}/bam/PROW_984_R20_CTACATCC-CTAGGTGA_align_sort_dm_mq20.bam
BAMFILE_12=${RESULTS_DIR}/bam/PROW_981_R20_CCACATTG-CTGGAGTA_align_sort_dm_mq20.bam



#### no missing data (maxmiss100) - ld001

sambamba depth region --regions=${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.prune.in.sambamba.bed ${BAMFILE_01} ${BAMFILE_02} ${BAMFILE_03} ${BAMFILE_04} ${BAMFILE_05} ${BAMFILE_06} ${BAMFILE_07} ${BAMFILE_08} ${BAMFILE_09} ${BAMFILE_10} ${BAMFILE_11} ${BAMFILE_12} -o ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001_sambamba_depth.txt



#### no missing data (maxmiss100) - ld99

sambamba depth region --regions=${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99.prune.in.sambamba.bed ${BAMFILE_01} ${BAMFILE_02} ${BAMFILE_03} ${BAMFILE_04} ${BAMFILE_05} ${BAMFILE_06} ${BAMFILE_07} ${BAMFILE_08} ${BAMFILE_09} ${BAMFILE_10} ${BAMFILE_11} ${BAMFILE_12} -o ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld99_sambamba_depth.txt



#### includes missing data - ld001

sambamba depth region --regions=${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld001.prune.in.sambamba.bed ${BAMFILE_01} ${BAMFILE_02} ${BAMFILE_03} ${BAMFILE_04} ${BAMFILE_05} ${BAMFILE_06} ${BAMFILE_07} ${BAMFILE_08} ${BAMFILE_09} ${BAMFILE_10} ${BAMFILE_11} ${BAMFILE_12} -o ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld001_sambamba_depth.txt



#### includes missing data - ld99

sambamba depth region --regions=${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld99.prune.in.sambamba.bed ${BAMFILE_01} ${BAMFILE_02} ${BAMFILE_03} ${BAMFILE_04} ${BAMFILE_05} ${BAMFILE_06} ${BAMFILE_07} ${BAMFILE_08} ${BAMFILE_09} ${BAMFILE_10} ${BAMFILE_11} ${BAMFILE_12} -o ${RESULTS_DIR}/vcf/${SPP}_plink_includemissing_LD_PCA_20230123/${SPP}_ALL_renamed_q30_minac1_maf05_dp10_ld99_sambamba_depth.txt









############################## bedtools - ignore #########################################



# this bedtools one is mostly working (not sure what the output means yet)...

# bedtools coverage -a ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/PROW_EXP_ALL_renamed_q30_minac1_maf05_dp10_maxmiss100_ld001.prune.in.bed -b ${RESULTS_DIR}/bam/PROW_953_R1_TGACAACC-CTGTTGAC_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_953_R5_CAGTGCTT-TGAGGTGT_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_954_R20_ACAGCAAG-TTGCAGAC_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_981_R1_TGCTGTGA-GAATCGTG_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_984_R1_TAGTGCCA-TCTTGACG_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_953_R20_TCACTCGA-GTGCCATA_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_954_R1_CTGACTAC-ACATTGCG_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_984_R5_TCGTGCAT-TGGTACAG_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_981_R5_CAACACAG-AGGCTTCT_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_954_R5_GTGATCCA-CGTTGAGT_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_984_R20_CTACATCC-CTAGGTGA_align_sort_dm_mq20.bam ${RESULTS_DIR}/bam/PROW_981_R20_CCACATTG-CTGGAGTA_align_sort_dm_mq20.bam > ${RESULTS_DIR}/vcf/${SPP}_plink_LD_PCA_20221014/PROW_EXP_bedtools_coverage.txt





