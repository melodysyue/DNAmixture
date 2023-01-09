library(tidyverse)

rm(list=ls())

chnk.hap <- read.table("./data/haplotypes/chinook.haplot.pF.20.0.2.txt", header=T)

length(unique(chnk.hap$id))
length(unique(chnk.hap$locus))
chnk.hap %>% 
  select(locus,haplo) %>% 
  unique() %>% 
  dim()

head(chnk.hap)
summary(chnk.hap)

chnk.hap %>% 
  select(locus, haplo, n.hap.perloc) %>% 
  unique() %>% 
  summary()


######################################
### Step 1: Convert to genepop file###
######################################

# Need to identify the homozygotes so we can duplicate the haplotype twice
homo2dup <- chnk.hap %>% 
  filter(n.hap.perLocInd==1) %>% 
  mutate(tmp_id = paste0(id,locus)) %>% 
  pull(tmp_id)

homo2add <- chnk.hap %>% 
  mutate(tmp_id = paste0(id,locus)) %>%
  filter(tmp_id %in% homo2dup) %>%
  mutate(rank.new = 2, read.ratio = NA) %>% 
  select(id, locus, haplo, n.hap.perloc, n.hap.perLocInd, depth, read.ratio, rank_new)

# bind rows so that each fish has 2 haplos per locus
colnames(chnk.hap)
colnames(homo2add) 
chnk_haplo_final <- rbind(chnk.hap,homo2add) %>% 
  mutate(haplo_info=paste0(locus,"_",haplo))

# double check there are two entries for each locus+ind combo
chnk_haplo_final %>% 
  group_by(locus, id) %>%
  count() %>% 
  ungroup() %>%
  distinct(n)

write.table(chnk_haplo_final, "./data/haplotypes/chinook.haplot_addHOMO.pF.txt", quote=F, row.names = F)


# before we go to wide two-column format we need to convert our haplotypes into numbers

haplo_key <- chnk_haplo_final %>% 
  distinct(haplo_info) %>% 
  mutate(haplo_id = 1:nrow(.))

hap_recode <- left_join(chnk_haplo_final, haplo_key, by = "haplo_info")

# now go wide to a 1 column format to calculate missing data
##hap2col <- hap_recode %>%
#  select(group,id,locus,rank.pF,haplo_id) %>%
#  unite(tmp, locus:rank.pF,sep = ".") %>%
#  spread(tmp,haplo_id)

hap1col <- hap_recode %>% 
  select(id, locus, haplo_id) %>% 
  group_by(id,locus) %>% 
  arrange(haplo_id) %>% 
  summarise(genotype=paste(haplo_id,collapse="/")) %>% 
  spread(locus, genotype) %>% 
  ungroup()

# calculate missing loci per individual
hap1col %>%
  mutate(n_miss = rowSums(is.na(.))) %>% 
  select(id, n_miss) %>%
  arrange(desc(n_miss))


# Let's print it out and do some analyses.
write.table(hap1col, file = "./data/haplotypes/chinook_haplot1col_genotype.txt", sep = "\t", row.names = F, quote = F)
write.table(haplo_key, file="./data/haplotypes/chinook_haplot_key.txt", sep="\t", row.names = F, quote=F)

######################################
###Step 2: Get haplotype frequency ###
######################################
#rm(list=ls())
#chnk_haplo_final <- read.table("./data/haplotypes/chinook.haplot_addHOMO.pF.txt", header=T)


chnk_haplo.count <- chnk_haplo_final %>% 
  group_by(locus) %>% 
  add_count(name="total_perLocus") %>% #total copies per locus;
  ungroup() %>% 
  group_by(locus, haplo_info) %>% 
  add_count(name="haplo.copies") %>% #how many copies per haplotype per locus;
  ungroup() %>% 
  select(id, locus, haplo, haplo_info, total_perLocus,n.hap.perloc, haplo.copies)


chnk_haplo.af <- chnk_haplo.count %>% 
  mutate(af=haplo.copies/total_perLocus) %>% 
  select(locus, haplo, haplo_info, total_perLocus, n.hap.perloc, haplo.copies, af) %>% 
  distinct()

chnk_haplo.af <- left_join(chnk_haplo.af, haplo_key, by="haplo_info") %>% 
  select(locus, haplo, haplo_info, haplo_id, total_perLocus, n.hap.perloc, haplo.copies, af)


write.table(chnk_haplo.af, "./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt", quote=F, row.names = F)
