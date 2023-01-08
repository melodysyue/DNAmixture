library(tidyverse)

rm(list=ls())
haplo_noNX_full <- read.table("./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", header=TRUE)
keep_N114 <- read.table("./data/panel/DNAMixture_loci2keep_N114.list", header=T)

##################################################
###Step 1: select only chinook fin clip samples###
##################################################


unique(haplo_noNX_full$group)

chnk <- haplo_noNX_full %>% 
  filter(group=="chin")

length(unique(chnk$id))
length(unique(chnk$locus))

unique(chnk$locus)[!unique(chnk$locus) %in% keep_N114$x]

####################
###Step 2: filter###
####################

d <- 20 #coverage of a haplotype at a locus within an individual
r <- 0.2 #read depth ratio of a haplotype within an individual

chnk %>% 
  filter(depth>=d & read.ratio>=r) %>% 
  pull(id) %>% 
  unique() %>% 
  length()

chnk.filter <- chnk %>% 
  filter(depth>=d & read.ratio>=r) %>% 
  group_by(locus) %>% 
  mutate(n.hap.perloc=n_distinct(haplo)) %>% 
  filter(n.hap.perloc>1) #remove monomorphic loci

length(unique(chnk.filter$locus))
length(unique(chnk.filter$id))


#remove loci with more than 2 haplotypes in any individual
dup.loci <- chnk.filter %>% 
  group_by(id, locus) %>% 
  mutate(n.hap.perLocInd=n_distinct(haplo)) %>% 
  filter(n.hap.perLocInd>2) %>% 
  pull(locus) %>% 
  unique()

chnk.filter2 <- chnk.filter %>% 
  filter(!locus %in% dup.loci) %>% 
  group_by(locus) %>% 
  mutate(n.hap.perloc=n_distinct(haplo)) %>% 
  ungroup() %>% 
  group_by(id, locus) %>% 
  mutate(n.hap.perLocInd=n_distinct(haplo)) %>% 
  select(id, locus, haplo, n.hap.perloc, n.hap.perLocInd, depth, read.ratio, rank_new) 

length(dup.loci)
length(unique(chnk.filter2$locus))
length(unique(chnk.filter2$id))


#iterative filtering by genotyping rate

#round 1, remove ind<30% and loci<10%

num.loci <- length(unique(chnk.filter2$locus)) #96
badInd <- chnk.filter2 %>% 
  group_by(id) %>% 
  mutate(genoInd=n_distinct(locus)) %>% 
  filter(genoInd<ceiling(num.loci*0.3)) %>%  
  pull(id) %>% 
  unique()

chnk.i30 <- chnk.filter2 %>% 
  filter(!id%in%badInd)

num.id <- length(unique(chnk.i30$id)) #906
badLoci <- chnk.i30 %>% 
  group_by(locus) %>% 
  mutate(genoLoc=n_distinct(id)) %>% 
  filter(genoLoc<ceiling(num.id*0.1)) %>% 
  pull(locus) %>% 
  unique()

chnk.i30.l10 <- chnk.i30 %>% 
  filter(!locus %in% badLoci) 


#round 2, remove ind<50% and loci<30%

num.loci <- length(unique(chnk.i30.l10$locus)) #95
badInd <- chnk.i30.l10%>% 
  group_by(id) %>% 
  mutate(genoInd=n_distinct(locus)) %>% 
  filter(genoInd<ceiling(num.loci*0.5)) %>%  
  pull(id) %>% 
  unique()

chnk.i30.l10.i50 <- chnk.i30.l10 %>% 
  filter(!id%in%badInd)


num.id <- length(unique(chnk.i30.l10.i50$id))
badLoci <- chnk.i30.l10.i50 %>% 
  group_by(locus) %>% 
  mutate(genoLoc=n_distinct(id)) %>% 
  filter(genoLoc<ceiling(num.id*0.3)) %>% 
  pull(locus) %>% 
  unique()

chnk.i30.l10.i50.l30 <- chnk.i30.l10.i50 %>% 
  filter(!locus %in% badLoci) 


#round 3, remove ind<70% and loci<50%

num.loci <- length(unique(chnk.i30.l10.i50.l30$locus)) #88
badInd <- chnk.i30.l10.i50.l30%>% 
  group_by(id) %>% 
  mutate(genoInd=n_distinct(locus)) %>% 
  filter(genoInd<ceiling(num.loci*0.7)) %>%  
  pull(id) %>% 
  unique()


chnk.i30.l10.i50.l30.i70 <- chnk.i30.l10.i50.l30 %>% 
  filter(!id%in%badInd)

num.id <- length(unique(chnk.i30.l10.i50.l30.i70$id)) #696
badLoci <- chnk.i30.l10.i50.l30.i70 %>% 
  group_by(locus) %>% 
  mutate(genoLoc=n_distinct(id)) %>% 
  filter(genoLoc<ceiling(num.id*0.5)) %>% 
  pull(locus) %>% 
  unique()

chnk.i30.l10.i50.l30.i70.l50 <- chnk.i30.l10.i50.l30.i70 %>% 
  filter(!locus %in% badLoci) 


#round 4, remove ind<80% and loci<70%

num.loci <- length(unique(chnk.i30.l10.i50.l30.i70.l50$locus)) #86

badInd <- chnk.i30.l10.i50.l30.i70.l50%>% 
  group_by(id) %>% 
  mutate(genoInd=n_distinct(locus)) %>% 
  filter(genoInd<ceiling(num.loci*0.8)) %>%  
  pull(id) %>% 
  unique()

chnk.i30.l10.i50.l30.i70.l50.i80 <- chnk.i30.l10.i50.l30.i70.l50 %>% 
  filter(!id%in%badInd)

num.id <- length(unique(chnk.i30.l10.i50.l30.i70.l50.i80$id)) #565

badLoci <- chnk.i30.l10.i50.l30.i70.l50.i80 %>% 
  group_by(locus) %>% 
  mutate(genoLoc=n_distinct(id)) %>% 
  filter(genoLoc<ceiling(num.id*0.7)) %>% 
  pull(locus) %>% 
  unique()

chnk.i30.l10.i50.l30.i70.l50.i80.l70 <- chnk.i30.l10.i50.l30.i70.l50.i80 %>% 
  filter(!locus %in% badLoci) 

#summary after filtering
length(unique(chnk.i30.l10.i50.l30.i70.l50.i80.l70$id)) #565
length(unique(chnk.i30.l10.i50.l30.i70.l50.i80.l70$locus)) #74
summary(chnk.i30.l10.i50.l30.i70.l50.i80.l70)
length(dup.loci)

chnk.i30.l10.i50.l30.i70.l50.i80.l70 %>% 
  ungroup() %>% 
  select(locus, haplo) %>% 
  unique() %>% 
  dim()

write.table(chnk.i30.l10.i50.l30.i70.l50.i80.l70, paste0("./data/haplotypes/chinook.haplot.pF.", d, ".", r, ".txt"), quote=F, col.names = TRUE, row.names = FALSE)

