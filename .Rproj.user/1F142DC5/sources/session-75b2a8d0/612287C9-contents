library(tidyverse)
library(formattable)
rm(list=ls())

align <- read.csv("./data/alignment/alignment_summary_filtered.csv", header = T)
mock_meta <- read.table("./data/meta/mock_meta.txt", header=T)
stool_meta <- read.table("./data/meta/stool_meta_n3.txt", header=T)

align %>% 
  count(group)

align <- align %>% 
  mutate(on_target_rate=total_ontarget/total_passed) 

align_filter <- align %>% 
  filter(group %in% c("Tissue", "Mock DNA mixture", "GI tract content")) %>% 
  droplevels()

table(align_filter$group)

align_filter$group <- as.factor(align_filter$group)
levels(align_filter$group)
align_filter$group <- factor(align_filter$group, levels=c("Tissue", "Mock DNA mixture", "GI tract content"))

median <- align_filter %>% 
  group_by(group) %>% 
  summarise(median=median(on_target_rate, na.rm = TRUE))

median

align_filter %>% 
  filter(group=="Tissue") %>% 
  filter(on_target_rate>0) %>% 
  summary()


pdf("./figures/fig1_on_target_rate.pdf", width = 12, height = 9)
align_filter %>% 
  ggplot(aes(x=on_target_rate, fill=group))+
  geom_histogram(alpha=0.8)+
  geom_vline(data=median, aes(xintercept=median), col="darkred", size=1)+
  geom_text(data=median, aes(x=0.02, y=150, label=paste0("Median = ", percent(round(median,4)))), size=5, col="darkred")+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~group, ncol=1)+
  labs(x="On-target rate",
       y="Frequency")+
  theme_bw(base_size = 20)+
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background  = element_blank(),
        strip.text = element_text(hjust=0))
dev.off()


#LMB and CCF
summary(align)
pdf("./figures/figS3_ontarget_LMB_CCF.pdf", width = 9, height = 6)
align %>% 
  filter(group %in% c("LMB", "CCF")) %>% 
  ggplot(aes(x=total_ontarget))+
  geom_histogram()+
  geom_vline(xintercept = 1100, col="darkred")+
  theme_bw(base_size=20)+
  xlab("Total number of on-target reads per sample")+
  ylab("Frequency")+
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))

dev.off()





