library(tidyverse)

dat = read.table("stabresults.csv", sep = ",")
colnames(dat) = c("n", "sim", "d", "nfinal", "subcom", "spp", "xstar", "lambda")

dataplot = dat %>% 
  group_by(n, sim, d, subcom) %>% 
  filter(spp == 0) %>% 
  mutate(stable = all(lambda < 0)) %>% 
  group_by(n, d, nfinal) %>% 
  count(name = "nfreq")

ggplot(dataplot, aes(x = d, fill = as.factor(nfinal))) +
  geom_bar(stat = "identity", 
           position = "fill",
           aes(y = nfreq))+
  scale_fill_brewer(palette = "YlOrRd")+
  facet_wrap(~n)
