library(tidyverse)
library(patchwork)

rm(list=ls())

df <- read.csv("./data/alignment/alignment_summary_filtered.csv", header=T)
stool_pf <- read.table("./data/haplotypes/stoolmix.haplot.pF.0.0.02.0.filtered.txt", header=T)
mock_pf <- read.table("./data/haplotypes/mockmix.haplot.pF.0.0.02.0.txt", header=T) %>% 
  mutate(id=str_replace_all(id, "-", "_"))

mock_pf %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(n_loc=n_distinct(locus),
            n_hap=n_distinct(locus_haplo)) %>% 
  left_join(df, by=c("id"="Sample")) %>% 
  summary()


p1 <- mock_pf %>% 
  group_by(id) %>% 
  summarise(n_loc=n_distinct(locus),
            n_hap=n_distinct(locus_haplo)) %>% 
  left_join(df, by=c("id"="Sample")) %>% 
  #ggplot(aes(x=total_passed, y=total_ontarget))+
  ggplot(aes(x=total_ontarget, y=n_loc))+
  geom_point(col="darkred", alpha=0.5)+
  scale_y_continuous(breaks = c(73,74,75))+
  labs(x="Total on-target reads",
       y="Number of loci genotyped",
       title="Mock DNA mixture samples")+
  theme_bw(base_size = 15)

p1


stool_pf %>% 
  group_by(id) %>% 
  summarise(n_loc=n_distinct(locus),
            n_hap=n_distinct(locus_haplo)) %>% 
  left_join(df, by=c("id"="Sample")) %>% 
  #filter(n_loc>67) %>% 
  summary()




p2 <- stool_pf %>% 
  group_by(id) %>% 
  summarise(n_loc=n_distinct(locus),
            n_hap=n_distinct(locus_haplo)) %>% 
  left_join(df, by=c("id"="Sample")) %>% 
  ggplot(aes(x=total_ontarget, y=n_loc))+
  geom_point(alpha=0.5,col="darkred")+
  gghighlight::gghighlight(n_loc>=15, unhighlighted_params = list(col="lightgray"))+
  scale_y_continuous(limits=c(0, 74), 
                     breaks=c(0, 20, 40, 60, 74))+
  labs(x="Total on-target reads",
       y="Number of loci genotyped",
       title="GI tract samples from the feeding trial")+
  theme_bw(base_size = 15)

p2
pdf("./figures/fig3_ontarget_mockmixVSstool.pdf", width = 12, height = 6)
p1 + p2 +
  plot_annotation(tag_levels = "a")
dev.off()