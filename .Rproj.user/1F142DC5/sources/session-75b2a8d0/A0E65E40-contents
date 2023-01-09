# MissR_geneflow
 
# How to cite:
Yue Shi, Kristen L. Bouska, Garrett J. McKinney, William Dokai, Andrew Bartels, Megan V. McPhee, Wesley A. Larson. Gene flow influences the genomic architecture of local adaptation in six riverine fish species. *Molecular Ecology* (2021). [https://www.biorxiv.org/content/10.1101/2021.05.18.444736v1](https://www.biorxiv.org/content/10.1101/2021.05.18.444736v1) (It just got accepted on Dec 1st, 2021. Link will be available soon!)

# Data

Input and intermediate data for various analyses for the six study species. We use the following code for species names.

  - **bhmw** -> Bullhead Minnow
  - **blgl** -> Bluegill
  - **fwdm** -> Freshwater Drum
  - **cncf** -> Channel Catfish
  - **gzsd** -> Gizzard Shad
  - **ersn** -> Emerald Shiner

`./data/popmap/`: popmap (after filtering) for each species.

`./data/env/`: 
  - `missR_env_bypop_reduced.csv`: standardized environmental data (subtract the mean and then devided by the standard deviation of the variable across populations via R `scale()` function, `center=TRUE, scale=TRUE`).
  
`./data/rowID_tag_tagpos/`: list of genotyped SNPs with row index, tag, and tag position. 

`./data/outlier_neutral/`: list of outlier SNPs and neutral SNPs with tag ($chrom) and tag position ($pos).
  - `*fstol.list`: *F<sub>ST</sub>* outliers.
  - `*geaol.list`: GEA outliers.
  - `*fst_gea_ol.list`: outlier SNPs (union of *F<sub>ST</sub>* outliers and GEA outliers).
  - `*neutral.list`: neutral SNPs (before thinning).
  - `*neutral_thin.list`: neutral SNPs (after thinning).
  - `*outlier_summary.txt`: identified outliers using each of the 7 methods, including Bayescan (bs), Arlequin (arl), OutFLANK (ofk), pcadapt, RDA(rda), LFMM2(lfmm), and Bayenv2 (bf).  

`./data/alignment/`: 
  - `*bwa.uni.map20.aln.txt`: list of aligned SNPs (after filtering) for each species.
  - `*bwa.aln.stats.txt`: chromosomal position, Fstp, Ho, maf for each aligned SNP (after filtering) for each species.
  - `*aln.hmm.islands.stats.txt`: alignment info plus HMM islands info.

`./data/hmm_inversions/`: 
 - `*aln.Fstp_3states_HMMstates.txt`: HMM states for each aligned SNP.
 - `*aln.fstp.txt`: Fstp value for each aligned SNP. 
 - `ersn.danRer11.ws20.*.ol.txt`: list of SNPs within each putative inversion.
 - `ersn.danRer11.ws20.*.ol.cluster.genotype.txt`: cluster membership (0, 1, 2) of individuals for each putative inversion. 
 - `ersn.chr*.inv.recode.vcf`: vcf file for each putative inversion.
 
`./data/geno_input/`:
 - `*rda.geno`: genotype file with the number of the non-reference allele (0, 1, 2) seen at each locus in each sample. 
 - `ersn*.gen`: genepop file for each putative inversion.

 
`./data/dxy/`:
 - `*.dxy.txt`: per-SNP dxy value. 

`./data/ld/`:
 - `*.ld`: R2 value for pairs of SNPs (maf>0.01) on each of the 5 chromosomes with over-clustered outlier SNPs. 


**Note**:

Demultiplexed RAD sequencing data used in this study are archived in the NCBI with BioProject ID PRJNA674918.

vcf files (all SNPs post filtering) and genepop files (neutral SNPs after thinning) are archived on DRYAD [(https://doi.org/10.5061/dryad.fn2z34tvx)](https://doi.org/10.5061/dryad.fn2z34tvx).

# Scripts

All analyses were performed in parallel within each species and the results were compared among species. The following scripts are demonstrated using one species as an example and for reference only unless indicated otherwise. Please adjust accordingly for your computing environment and working directory. 

`./scripts/0_read_processing/`: 
 - 00_process_radtags.slurm
 - 01_clone_filter.slurm
 - 02_ustacks.slurm
 - 03_cstacks.slurm
 - 04_sstacks.slurm
 - 05_stacks_TGP.slurm
 
`./scripts/1_snp_filtering/`:
 - 10_hdplot.r
 - 11_vcf_keep_highest_MAF.py
 - 12_countHetsMissing_genepop_sample-ncode3.pl
 
 
`./scripts/2_envVars_outlier_neutral/`:
 - 20_fig1C_figS1A_figS1B.r: Figure 1C and Figure S1.
 - 21_rda.r
 - 22_outlier_summary.r: Table S5.
 - 23_neutral_thinning.r
 - 24_neutral_pwfst.r: Table S6.
 - 25_neutral_PCA.r: Figure 3.
 
 
`./scripts/3_hmm_inversion/`:
 - 30_hmm_islands.r: filter hmm islands.
 - 31_lostruct.r: run the core lostruct algorithm.
 - 32_lostruct_ol.plot.r: export lostruct windows with distinct 3-cluster PCA patterns. 
 - 33_ld_byChr.r: LD decay plot function for Figure S2. 
 - 34_ld_heatmap.r: LD heatmap plot function for Figure 4D&E, Figure S3, Figure S4D&E, and Figure S5D&E.
 - 35_fig4ABC_figS4ABC_figS5ABC.r: Figure 4A&B&C, Figure S4A&B&C, and Figure S5A&B&C.
 - 36_inversion_genotype_heatmap.r: Figure S6. 
 - 37_fig5.r: Figure 5. 
 - 38_fig2.r: Figure 2. 
 
`./scripts/4_cluster/`:
 - 40_NND_permutation.r: Table S8.
 - 41_fst_ho_dxy_ld_permutation.r: Table S8. 
 - 42_fig6.r: Figure 6. 
  
`./scripts/permutation_functions/`: shared functions to conduct permutation tests. 
  


 


