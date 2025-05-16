library(tidyverse)
library(gridExtra)

dat = read.table("../data/results_simplest_extension_n_2.csv", sep = "")
# dat3 = read.table("../data/results_simplest_extension_n_3.csv", sep = "")
# dat4 = read.table("../data/results_simplest_extension_n_4.csv", sep = "")
# dat5 = read.table("../data/results_simplest_extension_n_5.csv", sep = "")
# dat6 = read.table("../data/results_simplest_extension_n_6.csv", sep = "")
# dat7 = read.table("../data/results_simplest_extension_n_7.csv", sep = "")
# dat8 = read.table("../../data/results_simplest_extension_n_8.csv", sep = "")

seed_i, n, d, pert_size, pert_i, alpha_i, δₚ, δₚₗᵢₙ, δₓ, δₓₗᵢₙ
#dat = rbind(dat2, dat3, dat4, dat5, dat6, dat7, dat8)
colnames(dat) = c("seed", "n", "d", "pert_size", "pert", "alpha", "delta", "delta_lin", "delta_x", "delta_x_lin")


datplot = dat %>% filter(delta < 4) %>% 
  mutate(alpha = factor(alpha, levels = c(0.1, 0.9)))

#load parameters

datpar = read.table("../data/parameters_n_2", sep = "")

###########################################################################
