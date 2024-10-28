library(tidyverse)

dat1 = read.table("../data/prop_feasible_states_1.csv")
dat2 = read.table("../data/prop_feasible_states_2.csv")
# dat3 = read.table("../data/feasibility_boundary_radius_3.csv")

dat = rbind(dat1, dat2)

colnames(dat) = c("prop_feas", "rho", "alpha", "n", "sim")

dataplot = dat %>% 
  group_by(n, alpha, rho) %>% 
  mutate(prop_feas_mean = mean(prop_feas))

ggplot(dataplot, aes(x = rho, y = alpha, fill = prop_feas_mean))+
  geom_tile()+
  facet_wrap(~n)+
  theme(aspect.ratio = 1)+
  scale_fill_distiller(palette = "RdBu")
