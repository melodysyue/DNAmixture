library(tidyverse)
library(gghighlight)
library(scales)


rm(list=ls())

map <- read.csv("./data/alignment/mapReadCount_CHIN_fin_initialTesting.txt", header=T)

panel <- read.table("./data/panel/panel_N125.list", header=T)

map.loci <- map %>% 
  gather(loci, mapped, starts_with("tag_id"),-Sample) %>% 
  spread(Sample, mapped) %>% 
  mutate(loci_sum = rowSums(select(., starts_with("Chin-FinClp")))) %>% 
  select(loci, loci_sum) %>% 
  filter(loci %in% panel$x)


#make loci rank plot based primer reads and primer.probe reads
map.loci <- map.loci[order(map.loci$loci_sum, decreasing = TRUE),] #sort
map.loci$loci <- factor(map.loci$loci, levels=map.loci$loci) #retain the order

max=59000
min=7000

loci.keep <- map.loci %>% 
  filter(loci_sum>=min & loci_sum<=max) %>% 
  pull(loci)

loci.remove <- panel %>% 
  filter(! x %in% loci.keep) %>% 
  pull(x)

write.table(loci.keep, "./data/panel/DNAMixture_loci2keep_N114.list", quote=F, row.names = F)

pdf("./figures/figS1_panel_optimization_N114.pdf", width = 9, height = 6)
map.loci %>% 
  ggplot(aes(x=loci, y=loci_sum))+
  geom_point(col="darkred", alpha=0.8)+
  gghighlight(loci_sum>= min & loci_sum <= max, unhighlighted_colour = "darkgray", use_direct_label = FALSE)+
  theme_bw(base_size = 15)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  xlab("Microhaplotype loci")+
  ylab("Sum of on-target reads across 377 tested individuals")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank()
  )

dev.off()