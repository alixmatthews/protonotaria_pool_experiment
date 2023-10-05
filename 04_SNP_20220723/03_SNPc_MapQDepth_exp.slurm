#!/bin/bash

#SBATCH --job-name=03_SNPc_MapQDepth_exp
#SBATCH --partition=pcon06
#SBATCH --output=03_SNPc_MapQDepth_exp_%j.txt
#SBATCH --error=03_SNPc_MapQDepth_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load samtools/1.15.1

cd /local_scratch/$SLURM_JOB_ID/


##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW_EXP

## FILENAMES
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_exp/20210816_exp_filenames_SNPpipeline.txt" 

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
## Loop through samples
## Skip alignments with MAPQ smaller than 20, calculate average depth, and index *align_sort_dm_mq20.bams

while IFS=" " read -r value1 value2
do {

rsync -av $RESULTS_DIR/bam/${value1}_align_sort_dm.bam .

# skips alignments with mapping quality less than 20
samtools view --threads 30 --min-MQ 20 $RESULTS_DIR/bam/${value1}_align_sort_dm.bam --bam --output $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam

# calculate average depth and output to a .txt file
# -a means output all positions
samtools depth --threads 30 -a $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average =",sum/NR; print "Stdev =",sqrt(sumsq/NR - (sum/NR)**2)}' >> $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20_depth.txt;

# index *_align_sort_dm_mq20.bam files (i.e., aligned, sorted, duplMarked, mapping quality of 20)
samtools index $RESULTS_DIR/bam/${value1}_align_sort_dm_mq20.bam



} done <"$filenames"






