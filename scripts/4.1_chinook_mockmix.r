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

##################
###MockMixture ###
##################

haplo_noNX_full <- read.table("./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", header=TRUE)
unique(haplo_noNX_full$group)

mock <- haplo_noNX_full %>% 
  filter(group=="mock")

mock.meta <- read.table("./data/meta/mock_meta.txt", header=T)

mock.meta$Sample_ID <- str_replace_all(mock.meta$Sample_ID, "_", "-")

mock <- left_join(mock, mock.meta, by=c("id"="Sample_ID"))

length(unique(mock$id))
length(unique(mock$locus))

#################################
###Filtering parameter set up ###
#################################
d <- 0 #coverage of a haplotype at a locus within an individual
r <- 0.02 #read depth ratio of a haplotype within an individual
c <- 0

##############################################################################################
### Step 1: filter based on the list of loci/haplotypes in the chinook source samples###
##############################################################################################

chnk.hap <- read.table("./data/haplotypes/chinook.haplot.pF.20.0.2.txt", header=TRUE)
loci.whitelist <- unique(chnk.hap$locus)
haplo_key <- read.table("./data/haplotypes/chinook_haplot_key.txt", header=TRUE)

mock.filter <- mock %>% 
  mutate(locus_haplo=paste0(locus, "_", haplo)) %>% 
  filter(locus %in% loci.whitelist) %>% 
  filter(locus_haplo %in% haplo_key$haplo_info)

table(mock.filter$locus_haplo %in% haplo_key$haplo_info)
table(mock.filter$locus %in% chnk.hap$locus)

length(unique(mock.filter$id))
length(unique(mock.filter$locus))


###############################
###Step 2: filter by quality###
###############################

mock.filter <- mock.filter %>% 
  group_by(locus, id) %>% 
  mutate(read.ratio=depth/depth[1]) %>% 
  filter(depth>=d & read.ratio>=r)

#update information
mock.filter <- mock.filter %>%
  group_by(locus,id) %>% 
  dplyr::mutate(n.hap.perLocInd=n_distinct(haplo),
         read.ratio.new=depth/depth[1],
         rank=row_number(),
         coverage_perLocInd=sum(depth)) %>% 
  ungroup() %>% 
  select(id, locus, haplo, locus_haplo, depth, coverage_perLocInd, n.hap.perLocInd, read.ratio, rank)

length(unique(mock.filter$id))
length(unique(mock.filter$locus))

summary(mock.filter) 

mock.filter <- mock.filter %>% 
  filter(coverage_perLocInd>=c)

length(unique(mock.filter$id))
length(unique(mock.filter$locus))


write.table(mock.filter, paste0("./data/haplotypes/mockmix.haplot.pF.", d, ".", r, ".", c, ".txt"), quote=F, col.names = TRUE, row.names = FALSE)


#############################
###Step 3: MLE estimation ###
#############################
##you need two information:
##y. : maximum putative NOC. The function will calculate from 1 to y. and get the one with MLE;
##loci: observed alleles, and their allele frequency for each locus. 

chnk_haplo.af <- read.table("./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt",header=TRUE)
mock.filter.af <- left_join(mock.filter,chnk_haplo.af,by= c("locus_haplo"="haplo_info", "locus", "haplo")) %>% 
  select(id, locus, haplo_id, af)


###MLE magic starts from here###
samples <- mock.filter.af %>% 
  distinct(id) %>% 
  pull()

mle <- data.frame(id=NA,seith_ori=NA)
max.n <- 20
j <- 1

for (i in samples) {
  pv <- mock.filter.af %>% 
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

summary(mle.full)
str(mle.full)

write.table(mle.full, "./data/NOC/mock_NOCestimate_0_0.02_0.txt", quote = F, row.names = F)

#plot it
chnk_haplo.af <- read.table("./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt",header=TRUE)
mle.full <- read.table("./data/NOC/mock_NOCestimate_0_0.02_0.txt", header=TRUE)
               
p1 <- chnk_haplo.af %>% 
  ggplot(aes(x=af))+
  geom_histogram(fill="lightgray", color="black", alpha=0.8)+
  #geom_histogram(aes(y=..count..),fill="lightgray", color="black", alpha=0.8)+
  theme_bw(base_size = 25)+
  labs(x="Haplotype frequency", y = "Count")+
  scale_x_continuous(breaks=c(seq(0,1,0.2)))+ 
  theme(panel.grid.minor = element_blank())
               
p2 <- mle.full %>% 
  group_by(n_ind) %>% 
  ggplot(aes(x=n_ind, y=mle_bias))+
  geom_jitter(color="darkgray", width = 0.25, alpha=0.5, size=4)+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red", alpha=0.8, size=1)+
  scale_x_continuous(breaks=c(2,3,5,7,9,10,15,20),
                     limits=c(0,20))+
  scale_y_continuous(breaks=c(-10,-5,-2.5,0,2.5,5))+
  theme_bw(base_size = 25)+
  xlab("True N")+
  ylab("Bias (Estimate - True)")+
  theme(panel.grid.minor = element_blank())

p1
p2

pdf(paste0("./figures/fig2_chnk_haplo_af_mockmix_est_", d, "_", r, "_", c, ".pdf"), width = 12, height = 9)
p1/p2+
  plot_layout(heights = c(1,2))+
  plot_annotation(tag_levels = "a")
dev.off()

mle.full %>% 
  group_by(n_ind) %>% 
  summarise(sd=sd(mle_bias))

mle.full %>% 
  filter(n_ind<11) %>% 
  summarise_at(vars(mle_bias), list(sd))


mle.full %>% 
  mutate(mle_bias_ab=case_when(
  mle_bias<0 ~ abs(mle_bias),
  TRUE ~ mle_bias
  )) %>% 
  filter(n_ind <= 10) %>% 
  summary()

