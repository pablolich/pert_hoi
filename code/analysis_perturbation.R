library(tidyverse)

dat = read.table("../data/perturbation.csv")
colnames(dat) = c("n", "alphavalue", "r0", "pertmag", "sim", "rpert", "ratioeqs")

#plot ratio as a function of alpha for different n

datplot = dat %>% 
  group_by(n, alphavalue, pertmag) %>% 
  filter(magpert == 0.01) %>% 
  summarise(avratio = mean(ratioeqs))

ggplot(datplot, aes(x = alphavalue, 
                    y = avratio,
                    color = as.factor(n))) +
  geom_point()+
  scale_y_continuous(trans = "log10")+
  theme(aspect.ratio = 1)+
  facet_wrap(~pertmag)

