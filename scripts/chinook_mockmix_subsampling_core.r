library(foreach)
library(doParallel)
library(doRNG) # needed for 'set seed' across parallelization capabilities
library(RcppAlgos) # faster combinations?
library(forensim)
library(gtools)
library(tidyverse)
library(patchwork)

rm(list=ls())
source("./scripts/SeithSA_DNAMixture.r")

NOC_subsample <- function(mock.filter, chnk_haplo.af, mock.meta, loc.list, d, r, c, prop){
  
  #subsampling loci
  loc.sub <- sample(loc.list, round(prop*length(loc.list)), replace=FALSE)
  mock.filter.sub <- mock.filter %>% 
    filter(locus %in% loc.sub)
  
  #filter by quality, d,r, c
  mock.filter.sub <- mock.filter.sub %>% 
    group_by(locus, id) %>% 
    mutate(read.ratio=depth/depth[1]) %>% 
    filter(depth>=d & read.ratio>=r)
  
  mock.filter.sub <- mock.filter.sub %>%
    group_by(locus,id) %>% 
    dplyr::mutate(n.hap.perLocInd=n_distinct(haplo),
                  read.ratio.new=depth/depth[1],
                  rank=row_number(),
                  coverage_perLocInd=sum(depth)) %>% 
    ungroup() %>% 
    select(id, locus, haplo, locus_haplo, depth, coverage_perLocInd, n.hap.perLocInd, read.ratio, rank)
  
  mock.filter.sub <- mock.filter.sub %>% 
    filter(coverage_perLocInd>=c)
  
  #MLE estimate
  mock.filter.sub.af <- left_join(mock.filter.sub,chnk_haplo.af,by= c("locus_haplo"="haplo_info", "locus", "haplo")) %>% 
    select(id, locus, haplo_id, af)
  
  samples <- mock.filter.sub.af %>% 
    distinct(id) %>% 
    pull()
  
  mle <- data.frame(id=NA,seith_ori=NA)
  max.n <- 20
  j <- 1
  
  for (i in samples) {
    pv <- mock.filter.sub.af %>% 
      filter(id==i) %>% 
      spread(locus,af) %>%
      select(-c(id, haplo_id)) %>% 
      as.list()
    pv <- lapply(pv,function(x) x[!is.na(x)])
    mle[j,]=c(i,
              MixtureLikelihood_multi(max.n, pv)$Contributors)
    j=j+1
  }
  
  ###MLE magic ends here####
  mle.full <- left_join(mle,mock.meta,by=c("id"="Sample_ID"))
  
  mle.full <- mle.full %>% 
    select(id, Sample_Well, n_ind=inds, mle=seith_ori, rep) %>% 
    mutate(mle=as.numeric(mle)) %>% 
    mutate(mle_bias=mle-n_ind)
  
  return(mle.full)
}




