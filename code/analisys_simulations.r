library(tidyverse)

data_all = read.table("../data/results_critical_radius_n_eq_all.csv")
data_follow = read.table("../data/results_critical_radius_n_eq_follow.csv")
colnames(data_all)<-c("sim", "n", "alpha", "rho", "n_eq_av", "mode")
colnames(data_follow)<-c("sim", "n", "alpha", "rho", "n_eq_av", "mode")
data_merged = rbind(data_all, data_follow)


toplot = data_merged %>% group_by(n, alpha, mode) %>% 
  summarize(rhoav = mean(rho),
            neq_av=mean(n_eq_av))

ggplot(toplot, aes(x = alpha, y = rhoav,
                   color = as.factor(n),
                   shape = as.factor(mode)))+
         geom_point()
