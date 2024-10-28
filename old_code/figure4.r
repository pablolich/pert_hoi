library(tidyverse)

datamaxorbit = read.table("../data/feasibility_boundary_radius.csv")
colnames(datamaxorbit) = c("sim", "n", "alpha", "rmax")

toplot = datamaxorbit %>% group_by(alpha) %>% mutate(rmaxav = mean(rmax))
ggplot(toplot)+
  geom_point(aes(x = alpha, y = rmaxav))
