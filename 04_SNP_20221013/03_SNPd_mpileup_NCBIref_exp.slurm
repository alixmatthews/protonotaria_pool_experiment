#!/bin/bash

#SBATCH --job-name=03_SNPd_mpileup_NCBIref_exp
#SBATCH --partition=comp06
#SBATCH --output=03_SNPd_mpileup_NCBIref_exp_%j.txt
#SBATCH --error=03_SNPd_mpileup_NCBIref_exp_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alix.matthews@smail.astate.edu  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00

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
## Do mpileup with supplying a list that has the absolute paths of the files
## -Ou flag indicates uncompressed bcf output to speed up performance by removing unnecessary compression/decompression 
## max-depth flag is related to: "at a position, read maximally INT reads per input file... Max raw per-file depth; avoids excessive memory usage". The default is 250, but they suggest the user has responsibility to change it based on their data. Our avg depth is around 12 and suggest to 2-3x the average depth, so maybe 35 is good. You can also filter based on depth later (using GATK hard filters). Point is to avoid really high coverages (above the average) because that can suggest misalignment and higher changes of erroneous variant calls.

bcftools mpileup -Ou --fasta-ref ${REF_DIR}/${REF}.fsa --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --bam-list ${RESULTS_DIR}/bam/${SPP}_align_sort_dm_mq20_bam.list --max-depth 35 > ${RESULTS_DIR}/bam/${SPP}_ALL_align_sort_dm_mq20_bam.mpileup






