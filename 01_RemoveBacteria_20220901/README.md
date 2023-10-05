## Address problematic sample (PROW 953 R1) (...and all other samples)

Obviously something is wrong with this sample, but because there is no clear pattern of problems with the other samples (R5 or R20) from this individual or from the other R1 samples, it is likely an isolated problem that is unique to this sample. Becuase we know that the reads will map to the PROW reference (I did this mapping step prior to recognizing that this sample does not have a 'correct' COI sequence), it is obviously still a MITE, and 99.9999999% likely it is a *Amerodectes protonotaria* mite (the other samples from this feather are *A. protonotaria* in the COI seqs and they all map to the reference. It may have an overabundance of bacteria in the reads that aTRAM was only able to recover bacterial reads in the COI gene. So the goal here is to remove the bacterial reads by mapping them to a bacteria reference, then taking the reads that DON'T map (unmapped reads, which are the non-bacterial reads) and putting those through aTRAM again.

---

**1.** Blast and download genomes

**1.1.** BLAST weird sample

- BLASTn of weird sample shows that it is most closely related to bacteria (Burkholderia and Paraburkholderia)
- top hit: Paraburkholderia hospita strain mHSR1 plasmid pmHSR1_P, complete sequence Sequence ID: CP024940.1 Length: 1235162 ... Identities: 1334/1428(93%)

**1.2.** Download genomes of these top hits bacteria

- *Concept:* Concept here is to use these bacteria as a reference to map the weird sample's reads to. I'll then extract the unmapped (i.e., non-bacterial) reads and re-run aTRAM

- Downloaded recent full genome references of the bacteria listed in **1.1.** and concatenate them together to use a reference (`Burkholderia_Paraburkholderia_reference.fasta`)

---

**2.** Deal with the bacteria reads

**2.1.** Map reads to bacteria reference
- Name file with only the werid sample: `20210816_exp_filenames_SNPpipeline_PROW953_only.txt`
- Slurm: `01_IndexRef_MapBacteria_Fastq2Bam_exp_PROW_953_R1.slurm`
- Flagstat summary of this sample against the bacteria reference: `PROW_953_R1_TGACAACC-CTGTTGAC_align_sort_flagstat_out.txt`
- - Shows that there is 51.63% of the reads mapping to these bacteria genomes! So no wonder only ~43% mapped to the feather mite genome! Still some that must be to something other than these bacteria and feather mite, but maybe removing these 51% of reads will be enough to extract the COI and get the species identity of this mite confirmed!

**2.2.** Extract unmapped (non-bacterial) reads
- Slurm: `01_ExtractUnmappedNonBacReads_exp_PROW_953_R1.slurm`

**2.3.** Convert unmapped .bam to R1 and R2 .fastqs
- Slurm: `01_ConvertUnmappedBam2Fastq_exp_PROW_953_R1.slurm`

---

**3.** I went ahead and mapped all samples to these bacteria references, just to compare the mapping rate
- Slurm: `01_IndexRef_MapBacteria_Fastq2Bam_exp_allsamples.slurm`
- Flagstat summary of all samples against the bacteria reference: `PROW_exp_allsamples_bacteria_flagstat_out.txt`
---


### Next, move to 01_aTRAM_COI_20220901 directory




