## Confirm species identity of 12 experimental samples 

**1.1.** aTRAM (round 1)
- filenames: `20210816_exp_filenames.txt`
- `01_aTRAM_COI_expmites_20220826.slurm`


**1.2.** Prep .fasta files
- Download all `stitched_exons.fasta` files
- cat all 12 .fastas together
- File with the 12 experimental mites: `Experimental_PROW_COI.fasta`

---

**2.1.** Align sequences
- Using MAFFT v7.503 on Linux machine

```
mafft --auto Experimental_PROW_plus_all_others_COI.fasta > ./iqtree_20220826/Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta
```

- Output: `Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta`
- One sample (PROW 953 R1) looks very weird. Let's BLAST it in **3.1.**

**2.2.** Run IQTree
- Using IQTree v. 1.6.12 on Linux machine

```
iqtree -s ./iqtree_20220826/Experimental_PROW_plus_all_others_COI.fasta.mafftaligned.fasta -m MFP -merit AICc -bb 1000 -pre COI_AICc -nt AUTO
```
- Ouput: `COI_AICc.treefile`
- Looking at output here, the weird sample (PROW 953 R1) is obviously very different!

---

### Next, move to 01_RemoveBacteria_20220901 directory

