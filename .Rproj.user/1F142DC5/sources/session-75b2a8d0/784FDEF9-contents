library(tidyverse)
library(scales)
rm(list=ls())

chnk.hap <- read.table("./data/haplotypes/chinook.haplot.pF.20.0.2.txt", header=TRUE)
haplo_key <- read.table("./data/haplotypes/chinook_haplot_key.txt", header=TRUE)
chnk_haplo.af <- read.table("./data/haplotypes/chinook.haplot_addHOMO.allelefreq.txt",header=TRUE)

summary(chnk_haplo.af)

length(unique(chnk_haplo.af$locus))
length(unique(chnk_haplo.af$haplo_id))

chnk_haplo.af %>% 
  filter(af<0.1) %>% 
  #filter(af<0.05) %>% 
  pull(locus) %>% 
  unique() %>% 
  length()

pdf("./figures/figS2_distribution_num_haplo.pdf", width = 9, height = 6)
chnk_haplo.af %>% 
  group_by(locus) %>% 
  summarise(n=n_distinct(haplo_info)) %>% 
  ggplot(aes(x=n))+
  geom_histogram(binwidth = 0.5)+
  scale_x_continuous(breaks = pretty_breaks())+
  labs(x="Number of haplotypes per locus", 
       y="Frequency")+
  theme_bw(base_size=15)
dev.off()



