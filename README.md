# DNAmixture

# How to cite:
Yue Shi, Cory M. Dick, Kirby Karpan, Diana Baetscher, Mark J. Henderson, Suresh A. Sethi, Megan V. McPhee, Wesley A. Larson. Towards absolute abundance for conservation applications: estimating the number of contributors via microhaplotype genotyping of mixed-DNA samples. *Molecular Ecology Resources* (2023) *in review*.

# Data

`./data/alignment/`: alignment summary

`./data/sam/`: list of sam files. See `./scripts/slurm/` for alignment pipeline. 

`./data/meta/`: sample meta info

`./data/panel/`: panel info

`./data/microhaplot_input/`: input files to run MICROHAPLOT. 

`./data/haplotypes/`: output files from MICROHAPLOT. 

`./data/NOC/`: NOC estimates. 

`./data/mock_subsamping_output/`: results of DNA mixture subsampling experiments
  
**Note**:

Demultiplexed GT-seq data used in this study are archived in the NCBI Sequence Read Archive with a BioProject ID, PRJNA917209. Consensus sequences (.fasta) and SNP info (.vcf ) of the microhaplotype panel with 125 markers used in the study to run MICROHAPLOT are archived on DRYAD (TBD). 

# Scripts

`./scripts/slurm/`: slurms files for alignment to generate the list of sam files in `./data/sam/`.  

`./scripts/`:
  - `1_panel_optimization.r`: Fig S1
  - `2_alignment_summary.r`: Fig 1 and Fig S3
  - `3.*_*.r`: Fig S2
    - run MICROHAPLOT to extract haplotypes from SAM files
    - conduct haplotype filtering
    - calculate haplotype frequency
    - check cross-contamination.
  - `4.*_chinook_mockmix*`: Fig 2, Fig 4, and Fig S4. 
  - `5.*_*.r`: Fig 3, Fig 5, and Fig S5. 
  - `SeithSA_DNAMixture.r`: likelihood-based model to estimate NOC in DNA mixtures. 
  - `chinook_mockmix_subsampling_core.r`: DNA mixture subsampling experiments.
  


 



