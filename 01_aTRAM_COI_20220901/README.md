## Now, let's figure out what this sample is

Now that we have the unmapped reads (i.e., those reads that did not map to the bacterial genome reference we gave it), we can run aTRAM and check the identity of this sample (and hopefully confirm it is *A. protonotaria*).

**1.** Pull COI from these unmapped reads and see what it is!

**1.1.** Run aTRAM
- Files needed: `COI_ref.fasta` and `20210816_exp_filenames_PROW953_only.txt`
- Slurm: `01_aTRAM_COI_unmappedPROW953R1_20220901.slurm`

**1.2.** BLASTn
- BLASTn of the `stitched_exons.fasta` output file from aTRAM
- It is now BLASTing to *A. protonotaria*!
  - Amerodectes sp. 3 AEM-2018 voucher STAR AM153 PROW330A cytochrome oxidase subunit I (COX1) gene, partial cds; mitochondrial Sequence ID: KY491605.1 Length: 1190... Identities: 1061/1070(99%)

---

**2.** Estimate phylogeny

**2.1.** Align sequences
- Using MAFFT v7.503 on Linux machine

```
mafft --auto Experimental_PROW_plus_all_others_COI.fasta > ./iqtree_20220901/Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta
```

- Output: `Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta`

**2.2.** Run IQTree
- Using IQTree v. 1.6.12 on Linux machine

```
iqtree -s ./iqtree_20220901/Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta -m MFP -merit AICc -bb 1000 -pre COI_AICc -nt AUTO
```
- Ouputs: `COI_AICc.treefile` and `COI_AICc.iqtree`
- Looking at output (`COI_AICc.treefile_bootstraps.pdf`)...... the sample is in the correct spot! We're good to go! 

---

### So, this confirms we've got all *A. protonotaria* samples. Can confidently move forward with the SNP results
