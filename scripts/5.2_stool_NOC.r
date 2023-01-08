library(foreach)
library(doParallel)
library(doRNG) # needed for 'set seed' across parallelization capabilities
library(RcppAlgos) # faster combinations?
library(forensim)
library(gtools)
library(tidyverse)
library(patchwork)
library(wesanderson)
library(scales)
library(broom)
library(ggsci)

rm(list=ls())
source("./scripts/SeithSA_DNAMixture.r")

##################
###stool ###
##################

haplo_noNX_full <- read.table("./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", header=TRUE)
unique(haplo_noNX_full$group)

stool <- haplo_noNX_full %>% 
  filter(group=="stool")

stool$id <- str_replace_all(stool$id, "-","_")

stool.meta <- read.csv("./data/meta/stool_meta.csv", header=T)

stool <- left_join(stool, stool.meta, by=c("id"="Sample_ID"))

head(stool)

length(unique(stool$id))
length(unique(stool$locus))
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

stool.filter <- stool %>% 
  mutate(locus_haplo=paste0(locus, "_", haplo)) %>% 
  filter(locus %in% loci.whitelist) %>% 
  filter(locus_haplo %in% haplo_key$haplo_info)

table(stool.filter$locus_haplo %in% haplo_key$haplo_info)
table(stool.filter$locus %in% chnk.hap$locus)

length(unique(stool.filter$id))
length(unique(stool.filter$locus))


###############################
###Step 2: filter by quality###
###############################

stool.filter <- stool.filter %>% 
  group_by(locus, id) %>% 
  mutate(read.ratio=depth/depth[1]) %>% 
  filter(depth>=d & read.ratio>=r)

#update information
stool.filter <- stool.filter %>%
  group_by(locus,id) %>% 
  dplyr::mutate(n.hap.perLocInd=n_distinct(haplo),
                read.ratio.new=depth/depth[1],
                rank=row_number(),
                coverage_perLocInd=sum(depth)) %>% 
  ungroup() %>% 
  select(id, locus, haplo, locus_haplo, depth, coverage_perLocInd, n.hap.perLocInd, read.ratio, rank)

length(unique(stool.filter$id))
length(unique(stool.filter$locus))

#filter by c
stool.filter <- stool.filter %>% 
  filter(coverage_perLocInd>=c)

length(unique(stool.filter$id))
length(unique(stool.filter$locus))
length(unique(stool.filter$locus_haplo))

write.table(stool.filter, paste0("./data/haplotypes/stoolmix.haplot.pF.", d, ".", r,".", c, ".txt"), quote=F, col.names = TRUE, row.names = FALSE)


#############################
###Step 3: MLE estimation ###
#############################
##you need two information:
##y. : maximum putative NOC. The function will calculate from 1 to y. and get the one with MLE;
##loci: observed alleles, and their allele frequency for each locus. 

chnk_haplo.af <- read.table("./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt",header=TRUE)

stool.filter.af <- left_join(stool.filter,chnk_haplo.af,by= c("locus_haplo"="haplo_info", "locus", "haplo")) %>% 
  select(id, locus, haplo_id, af)


###MLE magic starts from here###
samples <- stool.filter.af %>% 
  distinct(id) %>% 
  pull()

mle <- data.frame(id=NA,seith_ori=NA)
max.n <- 20
j <- 1

for (i in samples) {
  pv <- stool.filter.af %>% 
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
mle.full <- left_join(mle,stool.meta,by=c("id"="Sample_ID"))

mle.full <- mle.full %>% 
  select(id, n_ind=smolt, mle=seith_ori, species, temp, hour) %>% 
  mutate(mle=as.numeric(mle)) %>% 
  mutate(mle_bias=mle-n_ind)

write.table(mle.full, file=paste0("./data/NOC/stool_NOCestimate_",d,"_",r,"_", c, ".txt"), quote=F, row.names=F)


#plot it
rm(list=ls())
align <- read.csv("./data/alignment/alignment_summary_filtered.csv", header = T)
stool_meta <- read.table("./data/meta/stool_meta_n3.txt", header=T)

align <- align %>% 
  mutate(on_target_rate=total_ontarget/total_passed) 

align_filter_stool <- align %>% 
  filter(group=="GI tract content") %>% 
  left_join(stool_meta, by=c("Sample"="Sample_ID"))

align_filter_stool$hour <- as.factor(align_filter_stool$hour)
align_filter_stool$hour <- factor(align_filter_stool$hour, levels = c("6","12","24","36","48", "60", "72","84","96","120"))


mle.full <- read.table("./data/NOC/stool_NOCestimate_0_0.02_0.txt", header=T) %>% 
  mutate(id=str_replace_all(id, "_", "-"))
stool_sp <- read.table("./data/meta/stool_n3_sample_nloc15.list", header = T)


mle.full_filtered <- mle.full %>% 
  filter(id %in% stool_sp$x)

#boxplot
mle.full_filtered$hour <- as.factor(mle.full_filtered$hour)
mle.full_filtered$hour <- factor(mle.full_filtered$hour, levels = c("6","12","24","36","48", "60", "72","84","96","120"))


p1 <- align_filter_stool %>% 
  ggplot(aes(x=hour, y=on_target_rate, fill=species))+
  geom_boxplot(alpha=0.8)+
  scale_fill_manual(values=c("#E64B35FF", "#00A087FF"))+
  labs(x="", y="On-target rate", fill="Species")+
  facet_wrap(~temp)+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = c(0.95,0.8))

p2 <- mle.full_filtered %>% 
  na.omit() %>% 
  ggplot(aes(x=hour,y=mle, fill=species))+
  geom_boxplot(alpha=0.8)+
  scale_fill_manual(values=c("#E64B35FF", "#00A087FF"))+
  labs(x="Hours post-ingestion", y="Estimated NOC", fill="Species")+
  facet_wrap(~temp)+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = "none")

pdf("./figures/fig5_stool_ontarget_noc.pdf", width = 12, height = 9)
p1 / p2 +
  plot_annotation(tag_levels = "a")

dev.off()

mle.full_filtered %>% 
  na.omit() %>% 
  group_by(species, temp, hour) %>% 
  summarise(mean_mle=mean(mle)) %>% 
  spread(hour, mean_mle)



###exponential decay

mycol <- c("#F39B7FFF", "#E64B35FF","#91D1C2FF", "#00A087FF")

pdf("./figures/figS5_stool_NOCestimate_0_0.02_0_exponential_decay.pdf", width = 12, height = 9)

mle.full_filtered %>% 
  na.omit() %>% 
  mutate(hour=as.numeric(hour)) %>% 
  group_by(species, temp, hour) %>% 
  summarize(mean_mle=mean(mle)) %>% 
  mutate(temp=as.character(temp)) %>% 
  mutate(group=paste0(species, " / ", temp)) %>% 
  ggplot(aes(x=hour, y=mean_mle, color=group))+
  geom_smooth(method="nls",se=F,formula=y~a*b**x)+
  scale_x_continuous(breaks=c(6,12,25,36,48,60,72,84,96,120))+
  scale_color_manual(values=mycol)+
  labs(x="Hours post-ingestion", y="Mean NOC estimate", col=expression('Species / Temperature ('*~degree*C*')'))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank(),
        legend.position = c(0.85,0.9))

dev.off()

mle.full_filtered %>% 
  na.omit() %>% 
  group_by(species, temp, hour) %>% 
  summarize(mean_mle=mean(mle)) %>% 
  mutate(temp=as.character(temp)) %>% 
  mutate(group=paste0(species, " / ", temp)) %>% 
  group_by(group) %>% 
  group_map(~broom::tidy(nls(mean_mle~a*b**hour, data=.x)))

