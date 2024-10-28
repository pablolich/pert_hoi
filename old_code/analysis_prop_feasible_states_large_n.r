library(tidyverse)

dat = read.table("../data/prop_feasible_states_large_n.csv")
colnames(dat) = c("prop_feas", "rho", "alpha", "n", "sim")

datplot  = dat %>% 
  group_by(n, alpha) %>% 
  summarize(av_prop_feas = mean(prop_feas))

ggplot(datplot, aes(x = n, 
                y = av_prop_feas))+
  geom_line(aes(group = alpha,
                color = as.factor(alpha)))
