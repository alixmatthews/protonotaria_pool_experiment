#!/bin/bash

#SBATCH --job-name=01_IndexRef_MapBacteria_Fastq2Bam_exp_PROW_953_R1
#SBATCH --partition=comp01
#SBATCH --output=01_IndexRef_MapBacteria_Fastq2Bam_exp_PROW_953_R1_%j.txt
#SBATCH --error=01_IndexRef_MapBacteria_Fastq2Bam_exp_PROW_953_R1_%j.err
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

## MAKE SUBDIRS IN RESULTS_DIR
mkdir -p $RESULTS_DIR # necessary to make multiple dirs within
mkdir $RESULTS_DIR/bcf
mkdir $RESULTS_DIR/bam
mkdir $RESULTS_DIR/vcf
mkdir $RESULTS_DIR/sam
mkdir $RESULTS_DIR/flagstat




##----------------------------------------------------------------------------------------


## INDEX, .FAI, AND .DICT THE REFERENCES
## only have to do this once per reference
## REF_FULL
bwa index ${REF_DIR}/${REF}.fasta

## ensure .fai file is created
samtools faidx ${REF_DIR}/${REF}.fasta

## create sequence dictionary for reference
java -jar /share/apps/bioinformatics/picard/picard-tools-2.17.10/picard.jar CreateSequenceDictionary REFERENCE=${REF_DIR}/${REF}.fasta OUTPUT=${REF_DIR}/${REF}.dict




##----------------------------------------------------------------------------------------


## LOOP THROUGH SAMPLES
## fastq to indexed, align_sort.bam

while IFS=" " read -r value1 value2
do {

# align fastq reads to reference genome using BWA, output is a .sam file
# -M: tells bwa to consider split reads as secondary
# -R: read group info
bwa mem -t 32 -M -R "@RG\tID:${value1}\tSM:${value1}\tLB:${value2}\tPL:ILLUMINA" $REF_DIR/$REF.fasta $FASTQ/${value1}/${value1}_dd_bb_trim_R1.fastq $FASTQ/${value1}/${value1}_dd_bb_trim_R2.fastq > ${value1}_align.sam

# convert .sam to .bam file (bam is a compressed binary version of sam) using samtools
samtools view --threads 30 -S -b ${value1}_align.sam > ${value1}_align.bam

# sort the bam files using samtools
samtools sort --threads 30 ${value1}_align.bam -o ${value1}_align_sort.bam 

# index the align_sort.bam files
samtools index ${value1}_align_sort.bam

# mapping stats
samtools flagstat --threads 30 -O tsv ${value1}_align_sort.bam > ${value1}_align_sort_flagstat_out.txt


# rsync results back over to home directory
rsync -av --remove-source-files ${value1}_align.sam $RESULTS_DIR/sam/
rsync -av --remove-source-files ${value1}_align.bam $RESULTS_DIR/bam/
rsync -av --remove-source-files ${value1}_align_sort.bam $RESULTS_DIR/bam/
rsync -av --remove-source-files ${value1}_align_sort_flagstat_out.txt $RESULTS_DIR/flagstat/


} done <"$filenames"


