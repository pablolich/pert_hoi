library(tidyverse)
dat = read.table("../data/feasibility_radius_merged.csv")
colnames(dat) = c("sim", "n", "alpha", "rmax")

dataplot = dat %>% 
  group_by(n, alpha) %>% 
  summarise(rmean = mean(rmax))

ggplot(dataplot)+
  geom_point(aes(x = alpha, y = rmean,
                 color = as.factor(n)))
