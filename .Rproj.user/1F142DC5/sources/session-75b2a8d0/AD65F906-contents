library(foreach)
library(doParallel)
library(doRNG) # needed for 'set seed' across parallelization capabilities
library(RcppAlgos) # faster combinations?
library(forensim)
library(gtools)
library(tidyverse)
library(patchwork)

rm(list=ls())
source("SeithSA_DNAMixture.r")
source("chinook_mockmix_subsampling_core.r")

##################
###Input files ###
##################

haplo_noNX_full <- read.table("./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", header=TRUE)
mock <- haplo_noNX_full %>% 
  filter(group=="mock")

mock.meta <- read.table("./data/meta/mock_meta.txt", header=T)
mock.meta$Sample_ID <- str_replace_all(mock.meta$Sample_ID, "_", "-")
mock <- left_join(mock, mock.meta, by=c("id"="Sample_ID"))

chnk_haplo.af <- read.table("./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt",header=TRUE)
chnk.hap <- read.table("chinook.haplot.pF.20.0.2.txt", header=TRUE)
loci.whitelist <- unique(chnk.hap$locus)
haplo_key <- read.table("chinook_haplot_key.txt", header=TRUE)

#filter based on the list of loci/haplotypes in the chinook source samples###
mock.filter <- mock %>% 
  mutate(locus_haplo=paste0(locus, "_", haplo)) %>% 
  filter(locus %in% loci.whitelist) %>% 
  filter(locus_haplo %in% haplo_key$haplo_info)

loc.list <- unique(mock.filter$locus)

#################################
###Filtering parameter set up ###
#################################
d <- 0 #coverage of a haplotype at a locus within an individual
r <- 0.002 #read depth ratio of a haplotype within an individual
c <- 0 #coverage per loc per sample


#################################################
###run subsampling function and repeat N times###
#################################################
n=50 #repeat n times for each subsampling
prop.list <- seq(0.1, 0.9, 0.1) #the number of loci


for (prop in prop.list){
  print(paste0("Start proportion ", prop))
  
  mle.full.sim <- NULL
  
  for (i in 1:n){
    print(paste0("Simulation", i))
    mle.full <- NOC_subsample(mock.filter, chnk_haplo.af, mock.meta, loc.list, d, r, c, prop)
    mle.full$sim=i
    mle.full.sim <- rbind(mle.full.sim, mle.full)
    }
  write.table(mle.full.sim, paste0("./mock_subsampling_output/r_", r, "_prop_", prop, "_sim_", n, ".txt"), quote=F, row.names = F)

}




