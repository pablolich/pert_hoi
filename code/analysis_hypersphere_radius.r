library(tidyverse)

dat1 = read.table("../data/feasibility_boundary_radius_1.csv")
dat2 = read.table("../data/feasibility_boundary_radius_2.csv")
dat3 = read.table("../data/feasibility_boundary_radius_3.csv")

dat = rbind(dat1, dat2, dat3)

colnames(dat) = c("sim", "n", "alpha", "rmax")

dataplot = dat %>% 
  group_by(n, alpha) %>% 
  summarise(rmean = median(rmax))

ggplot(dataplot)+
  geom_point(aes(x = alpha, y = rmean,
                 color = as.factor(n)))
