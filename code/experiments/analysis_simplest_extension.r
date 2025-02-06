library(tidyverse)

dat2 = read.table("../../data/results_simplest_extension_n_2.csv", sep = "")
dat3 = read.table("../../data/results_simplest_extension_n_3.csv", sep = "")
dat4 = read.table("../../data/results_simplest_extension_n_4.csv", sep = "")
dat5 = read.table("../../data/results_simplest_extension_n_5.csv", sep = "")


dat = rbind(dat2, dat3, dat4, dat5)
colnames(dat) = c("seed", "n", "alpha", "pert_dir", "delta")
datplot = dat %>% filter(delta < 2)

ggplot(datplot, aes(x = delta,
                fill = as.factor(alpha),
                after_stat(density))) +
  geom_histogram(position="identity", alpha = 0.5)+
  theme(aspect.ratio = 1)+
  facet_wrap(~n,
             nrow = 2,
             scales = "free")

