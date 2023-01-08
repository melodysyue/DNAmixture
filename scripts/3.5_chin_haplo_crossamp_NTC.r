library(tidyverse)

rm(list=ls())
haplo_noNX_full <- read.table("./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", header=TRUE)

unique(haplo_noNX_full$group)

d <- 20 #coverage of a haplotype at a locus wtihin an individual
r <- 0.2 #read depth ratio of a haplotype within an individual

##################################################################
###Step 1: cross amplification in  LMB and CCF fin clip samples###
##################################################################

other <- haplo_noNX_full %>% 
  filter(group %in% c("lmb", "ccf"))


other.filter <- other %>% 
  filter(depth>=d & read.ratio>=r) %>% 
  mutate(locus_haplo=paste0(locus, "-", haplo)) %>% 
  select(id, locus, haplo, locus_haplo, depth) 

write.table(other.filter, paste0("./data/haplotypes/lmb.ccf.haplot.pF.", d, ".", r, ".txt"), quote=F, col.names = TRUE, row.names = FALSE)


#########################################
###Step 2: check contamination in NTCs###
#########################################

ntc <- haplo_noNX_full %>% 
  filter(group == "ntc")

length(unique(ntc$id))

summary(ntc)

ntc.filter <- ntc %>% 
  filter(depth>=d & read.ratio>=r) %>% 
  mutate(locus_haplo=paste0(locus, "-", haplo)) %>% 
  select(id, locus, haplo, locus_haplo, depth) 

write.table(ntc.filter, paste0("./data/haplotypes/ntc.haplot.pF.", d, ".", r, ".txt"), quote=F, col.names = TRUE, row.names = FALSE)

###################################
###compare with the chinook data###
###################################

other.filter <- read.table("./data/haplotypes/lmb.ccf.haplot.pF.20.0.2.txt", header=T)

other.list <- unique(other.filter$locus_haplo)


chnk <- read.table("./data/haplotypes/chinook.haplot.pF.20.0.2.txt", header=T)

chnk.list <- chnk %>% 
  mutate(locus_haplo=paste0(locus, "-", haplo)) %>% 
  pull(locus_haplo) %>% 
  unique()

cross.amp <- other.filter %>% 
  filter(locus_haplo %in% chnk.list) 


length(unique(cross.amp$id))
length(unique(cross.amp$locus))
length(unique(cross.amp$locus_haplo))

length(unique(other.filter$id))
other.filter %>% 
  #filter(locus_haplo %in% chnk.list) %>% 
  group_by(id) %>% 
  mutate(total_cov=sum(depth)) %>% 
  select(id, total_cov) %>% 
  unique() %>% 
  arrange(desc(total_cov))





