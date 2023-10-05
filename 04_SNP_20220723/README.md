# SNP calling pipeline

**Note:** for all of these slurms, we are using the PROW981_R1 "reduced full reference with microbes removed" 
(how we got to that point with that reference assembly is outlined in 20210816_projects/20210816_snp/20220513_snp/02_IndexRef/)

---

## 1. Convert .fastq to indexed, align_sort.bam files 

- slurm: `03_SNPa_Fastq2Bam_exp.slurm` 
- Get flagstat summary (output: `03_SNPa_Fastq2Bam_exp_FLAGSTAT_SUMMARY.xlsx`)

```
grep "" *.txt > PROW_EXP_flagstat_out_all_files.txt
# import to Excel
# copy and paste first couple of lines and do flash fill to separate #s following : from the sample name
# use filters to filter our what category want to look at and summarize the output
```



## 2. Mark duplicates and index the align_sort_dm.bams
- slurm: `03_SNPb_MrkDup_exp.slurm`



## 3. Skip alignments with MAPQ smaller than 20, calculate average depth, and index align_sort_dm_mq20.bams 
- slurm: `03_SNPc_MapQDepth_exp.slurm`
- summary is here (produced using code below): `03_SNPc_MapQDepth_exp_output.xlsx`
- get a list of bam files for mpileup (next step):

```
cd /scrfs/storage/amatthews/20210816_projects/20210816_exp/04_SNP_20220723/RESULTS_PROW_EXP/bam
find "$PWD" -type f -name '*mq20.bam' > PROW_EXP_align_sort_dm_mq20_bam.list


# grep the avg and stdev: 
grep "Average" *mq20_depth.txt
grep "Stdev" *mq20_depth.txt
```
    


## 4. Create mpileup files, max read depth per file is set at 35
- slurm: `03_SNPd_mpileup_exp.slurm`



## 5. Call variants on the mpileup file, joint-genotyping
- slurm: `03_SNPe_callvars_exp.slurm`
    + need `PROW_EXP_reheader.txt` with this slurm
- Check the order of the names in the .vcf and then will want to rename to make sure they go in the same order as the original

 ```
## Rename samples
cd 
/scrfs/storage/amatthews/20210816_projects/20210816_exp/04_SNP_20220723/RESULTS_PROW_EXP/vcf

module load python/anaconda-3.9
source /share/apps/bin/conda-3.9.sh
conda activate BCFTools # version 1.15.1

## First, check the order using the following once this slurm is complete:
bcftools query -l PROW_EXP_ALL.vcf 

## Then create a file with the reheaders of the samples in the correct order (e.g., ${SPP}_reheader.txt).

## Then run this to conduct the reheader'ing:
bcftools reheader --samples /scrfs/storage/amatthews/20210816_projects/20210816_exp/04_SNP_20220723/PROW_EXP_reheader.txt /scrfs/storage/amatthews/20210816_projects/20210816_exp/04_SNP_20220723/RESULTS_PROW_EXP/vcf/PROW_EXP_ALL.vcf --output /scrfs/storage/amatthews/20210816_projects/20210816_exp/04_SNP_20220723/RESULTS_PROW_EXP/vcf/PROW_EXP_ALL_renamed.vcf

## Last, check the of the reheader'd samples and make sure it is the same as the original:
bcftools query -l PROW_EXP_ALL_renamed.vcf 
 ```



## 6. Filter by phred QUAL (threshold is 30)
- slurm: `03_SNPf_filterqual_exp.slurm`



## 7. Filter out SNPs that differ from the reference but are equal to each other across samples
- slurm: `03_SNPg_filterminac_exp.slurm`



## 8. Filter on minor allele freq (MAF), several thresholds
- slurm: `03_SNPh_filterminaf_exp.slurm`



## 9. Filter by sampling seq depth (FORMAT/DP), several thresholds
- slurm: `03_SNPi_filterdp_exp.slurm`



## 10. Filter out SNP sites with missing data, 100% only 
(other thresholds are possible, but not really probably a good idea to do dif ones with this particular project)
- slurm: `03_SNPj_filtermissing_exp.slurm`



### File with overview of each filtering step: `PROW_EXP_NumberOfSNPs_filtering.xlsx`


---

## 13. LD-pruning; move to next directory (`05_Analyses`)
- these slurms are in `05_Analyses/05_plink_LD_PCA_$DATE`


### Directory with the number of SNPs ~ linkage threshold `PROW_EXP_NumSNPs_LD` contains summary of summary file, R script, and figures

