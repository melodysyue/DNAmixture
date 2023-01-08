library(tidyverse)
library(patchwork)

rm(list=ls())


###input###
###look through the parameter space of r (read depth ratio of a haplotype within an individual) and the proportion of loci out of 125 loci. 

files <- list.files(path="./data/mock_subsampling_output/", pattern="*_sim_50.txt", full.names = TRUE) %>% 
  str_sort(numeric = FALSE)


mock_df <- files %>% 
  map_df(read.table, header=T)


r_list <- list.files(path="./data/mock_subsampling_output/", pattern="*_sim_50.txt", full.names = TRUE) %>% 
  str_sort(numeric = FALSE) %>% 
  basename() %>% 
  str_remove_all(., "_sim_50.txt")

r_list <- sapply(strsplit(r_list, "_"), '[', 2)



prop_list <- list.files(path="./data/mock_subsampling_output/", pattern="*_sim_50.txt", full.names = TRUE) %>% 
  str_sort(numeric = FALSE) %>% 
  basename() %>% 
  str_remove_all(., "_sim_50.txt")

prop_list <- sapply(strsplit(prop_list, "_"), '[', 4)


head(mock_df) #282 mock samples * 50 simulations * 27 parameters

n <- length(files)
mock_df$r <- rep(r_list, each=nrow(mock_df)/n)
mock_df$prop <- rep(prop_list, each=nrow(mock_df)/n)

head(mock_df)

table(mock_df$r, mock_df$prop)
table(mock_df$r, mock_df$n_ind)


###loci list
n.loci <- 74

df_prop <- NULL
for(p in seq(0.1, 0.9, 0.1)){
  print(paste0("Proportion ", p, " has ", round(p*n.loci), " loci"))
  df <- data.frame(p=p, n_loci=round(p*n.loci))
  df_prop <- rbind(df_prop, df)
}

df_prop <- df_prop %>% 
  mutate(label=paste0(p*100, "% (", n_loci, " loci)"))


###plot it###
label <- df_prop$label
names(label) <- seq(0.1,0.9,0.1)





p <- mock_df %>% 
  ggplot(aes(x=n_ind, y=mle_bias))+
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), geom="pointrange", color="red", alpha=0.8, size=1)+
  scale_x_continuous(breaks=c(2,3,5,7,9,10,15,20),
                     limits=c(0,20))+
  scale_y_continuous(breaks=c(-10,-5,-2.5,0,2.5,5))+
  theme_bw(base_size = 20)+
  xlab("True N")+
  ylab("Bias (Estimate - True)")+
  theme(panel.grid.minor = element_blank())+
  facet_grid(r ~ prop, labeller = labeller(prop=label))


ggsave(filename = "./figures/fig4_chnk_haplo_af_mockmix_subsampling_byProp.png", p,  width=36, height=12, units="in")


#variance and mean

p_var <- mock_df %>% 
  group_by(r, prop, n_ind) %>% 
  summarize(var_mle=var(mle)) %>% 
  ggplot(aes(x=prop, y=as.factor(n_ind), fill=var_mle))+
  geom_tile()+
  scale_fill_gradient(low="aliceblue", high="navy", name="Variance of \nNOC estimate")+
  scale_x_discrete(labels=label)+
  labs(x="Simulated Genotyping Rate",
       y="True NOC")+
  facet_grid(.~r)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.background = element_rect(fill="white"))


p_var

p_bias <- mock_df %>% 
  group_by(r, prop, n_ind) %>% 
  summarize(mean_bias=mean(mle_bias))%>% 
  ggplot(aes(x=prop, y=as.factor(n_ind), fill=mean_bias))+
  geom_tile()+
  scale_fill_gradient(low="salmon", high="white",
                      name="Mean Bias\n(Estimate - True)")+
  scale_x_discrete(labels=label)+
  labs(x="Simulated Genotyping Rate",
       y="True NOC")+
  facet_grid(.~r)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.background = element_rect(fill="white"))


p_bias

ylab <- p_var$labels$y
p_var$labels$y <- p_bias$labels$y <- " "

pdf("./figures/figS4_subsampling_var_mean.pdf", width = 12, height = 9)
((p_var+theme(axis.text.x = element_blank(), 
             axis.ticks.x = element_blank(),
             axis.title.x = element_blank()))/p_bias)
  grid::grid.draw(grid::textGrob(ylab, x=0.02, y=0.55,rot = 90))
dev.off()

mock_df %>% 
  group_by(r, prop, n_ind) %>% 
  summarize(var_mle=var(mle)) %>% 
  filter(r==0.002) %>% 
  spread(n_ind, var_mle)


mock_df %>% 
  group_by(r, prop, n_ind) %>% 
  summarize(mean_bias=mean(mle_bias))%>% 
  filter(r==0.002) %>% 
  spread(n_ind, mean_bias)



