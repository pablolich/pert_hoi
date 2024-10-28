library(tidyverse)

dat = read.table("../data/portrait.csv")
colnames(dat) = c("sim", "n", "constr", "alpha", "spp", "eqID", "xRe", "lambdaRe")

datplot = dat %>% 
  filter(spp == 1) %>% 
  group_by(sim, n, constr, eqID, alpha) %>% 
  mutate(feasible = all(xRe > 0),
         stable = lambdaRe < 0) %>% 
  filter(feasible == T) %>% 
  group_by(sim, n, constr, alpha, stable) %>% 
  tally(name = "neq") %>% 
  group_by(n, constr, alpha, stable) %>% 
  summarise(neqav = mean(neq))

ggplot(datplot)+
  geom_point(aes(x = alpha, y = neqav,
                 color = as.factor(stable)))+
  scale_color_manual(values = c("gray", "orange"))+
  scale_shape_manual(values = c(20,3))+
  facet_grid(constr~n)

