rm(list=ls())
library(microhaplot)
library(tidyverse)
##Establish shinny app

shiny.path <- "./shiny/"
microhaplot::mvShinyHaplot(shiny.path)


run.label="DNAmixture"
sam.path="./sam"

#prepare label file
#label <- read.table("./data/microhaplot_input/sam.list", header=F)
#label$V1 <- basename(as.character(label$V1))
#head(label)

#label <- label %>% 
#  mutate(sample=sapply(strsplit(V1,"_"), "[", 1)) %>% 
#  mutate(group=case_when(
#    str_detect(sample, "Chin-FinClp") ~ "chin", 
#    str_detect(sample, "LMB-FinClp") ~ "lmb",
#    str_detect(sample, "CCF-FinClp") ~ "ccf", 
#    str_detect(sample, "H12-NTC") ~ "ntc",
#    str_detect(sample, "MockMixture") ~ "mock",
#    str_detect(sample, "Stool-Chin") ~ "stool", 
#    TRUE ~ "ntc")) %>% 
#  select(sam=V1, sample, group)

#table(label$group)

#write.table(label, "./data/microhaplo_input/CA_label_group.txt", sep="\t", quote=F, col.names = F, row.names = F)


label <- read.table("./data/microhaplot_input/CA_label_group.txt",header=FALSE)
label.path <- file.path("./data/microhaplot_input/CA_label_group.txt")
vcf.path <- file.path("./data/microhaplot_input/chinook_ultima3_q30dp10mac3.recode.vcf")
app.path <- file.path(shiny.path, "microhaplot")


haplo.read.tbl <- prepHaplotFiles(run.label=run.label,
                                  sam.path = sam.path,
                                  label.path = label.path,
                                  vcf.path = vcf.path,
                                  app.path = app.path)

write.table(haplo.read.tbl, "./data/haplotypes/microhaplotype_readcount_CA.txt",quote=FALSE, row.names=FALSE) #this is a big file. 


###Filtering
haplo.read.tbl <- read.table("./data/haplotypes/microhaplotype_readcount_CA.txt", header=TRUE)
head(haplo.read.tbl)
length(unique(haplo.read.tbl$id))
length(unique(haplo.read.tbl$locus))

haplo.read.tbl %>% 
  filter(group=="chin") %>% 
  pull(id) %>% 
  unique() %>% 
  length()


##Step 1: Remove hapolotypes with N and X
haplo.read.tbl %>% 
  filter(grepl("X|N", haplo))

haplo_noNX <- haplo.read.tbl%>% 
  filter(!grepl("X|N",haplo)) 


num.id <- length(unique(haplo_noNX$id)) #2286
num.loci <- length(unique(haplo_noNX$locus)) #115

haplo_noNX$id %>% unique() %>% sort()


###Step 2: Recalculate read ratio and rank, since some haplotypes with N or X could have rank of 1.

haplo_noNX_recal <- haplo_noNX %>% 
  group_by(locus,id) %>% 
  mutate(n.haplot.perLocInd=n(),
         read.ratio=depth/depth[1],
         rank_new=row_number(),
         coverage_perLocInd=sum(depth)) %>% ###depth[1] is the depth of the allele with highest read count. depth[2] is the depth for the second, etc
  ungroup() %>% 
  select(group, id, locus, haplo, depth, coverage_perLocInd, n.haplot.perLocInd, read.ratio, rank_new)


###Step 3: Add some other stats

haplo_noNX_full <- haplo_noNX_recal %>% 
  group_by(locus) %>% 
  mutate(n.ind.perLoc=length(unique(id)), minn.hap.perLoc=min(n.haplot.perLocInd), maxn.hap.perLoc=max(n.haplot.perLocInd), 
         avgCoverage_perLoc=sum(depth)/num.id) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  mutate(n.loci.perInd=length(unique(locus)), avgCoverage_perInd=sum(depth)/num.loci) %>% 
  ungroup()

write.table(haplo_noNX_full, "./data/haplotypes/microhaplotype_readcount_CA_noNX.txt", quote=FALSE, row.names=FALSE)


  
  
