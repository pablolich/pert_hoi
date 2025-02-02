library(tidyverse)

dat2 = read.table("../../data/results_simplest_extension_n_2.csv", sep = "")
dat3 = read.table("../../data/results_simplest_extension_n_3.csv", sep = "")
dat4 = read.table("../../data/results_simplest_extension_n_4.csv", sep = "")

dat = rbind(dat2, dat3, dat4)
colnames(dat) = c("seed", "n", "alpha", "pert_dir", "delta")

ggplot(dat) +
  geom_density(aes(x = delta,
                   color = as.factor(alpha)))+
  facet_wrap(~n,
             nrow = 2,
             scales = "free")
