#!/bin/bash

#SBATCH --job-name=01_ConvertUnmappedBam2Fastq_exp_PROW_953_R1
#SBATCH --partition=comp01
#SBATCH --output=01_ConvertUnmappedBam2Fastq_exp_PROW_953_R1_%j.txt
#SBATCH --error=01_ConvertUnmappedBam2Fastq_exp_PROW_953_R1_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00

##----------------------------------------------------------------------------------------
## LOAD MODULES

module load bwa/0.7.17
module load samtools/1.15.1
module load java/sunjdk_1.8.0 # for gatk and picard
module load gatk/4.2.6.1
module load fastqc/0.11.5 # workaround for picard-tools
module load picard-tools/2.17.10


cd /local_scratch/$SLURM_JOB_ID/

##----------------------------------------------------------------------------------------
## "ADJUSTABLES"

## SPECIES OF INTEREST
SPP=PROW_EXP_PROW_953_R1 

## FILENAMES
filenames="/scrfs/storage/amatthews/20210816_projects/20210816_exp/20210816_exp_filenames_SNPpipeline_PROW953_only.txt"

## MAIN DIR
PROJECT_DIR=/scrfs/storage/amatthews/20210816_projects/20210816_exp


## FASTQ LOCATIONS
FASTQ=$PROJECT_DIR/00_PP_20210903/Adapter_Removed_bb

## REF DIR, REF, AND REFTYPE
# it's a version of this reference, best assembly from spades on 20211014 of PROW 981 R1: /scrfs/storage/amatthews/20210816_projects/20210816_exp/01_ASSEMBLY/spades_20211014/PROW_981_R1_TGCTGTGA-GAATCGTG/scaffolds.fasta.... but I have removed the microbe seqs and all contigs <5000bp (scaffolds_reduced_contigs_kept.fasta)


REF_DIR=$PROJECT_DIR/01_RemoveBacteria_20220901
REF=Burkholderia_Paraburkholderia_reference


## SNP DIR (mkdir beforehand)
SNP_DIR=$REF_DIR/01_SNP_20220901

## RESULTS DIR PER SPECIES
RESULTS_DIR=$SNP_DIR/RESULTS_$SPP




##----------------------------------------------------------------------------------------
## convert .bam to .fastq

cd $RESULTS_DIR/bam/


while IFS=" " read -r value1 value2
do {

# sort unmapped (sort by name -n)
samtools sort --threads 30 -n ${value1}_align.unmapped.bam -o ${value1}_align_sort.unmapped.bam

# save fastq reads into separate R1 and R2 files
samtools fastq --threads 30 ${value1}_align_sort.unmapped.bam -1 ${FASTQ}/${value1}/${value1}_dd_bb_trim_unmapped_R1.fastq -2 ${FASTQ}/${value1}/${value1}_dd_bb_trim_unmapped_R2.fastq -0 /dev/null -s /dev/null -n


} done <"$filenames"


