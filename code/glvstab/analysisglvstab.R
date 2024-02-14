library(tidyverse)

dat = read.table("stabresults.csv", sep = ",")
colnames(dat) = c("n", "sim", "d", "nfinal", "subcom", "spp", "xstar", "lambda",
                  "invasible")

dataplot = dat %>% 
  group_by(n, sim, d, subcom) %>% 
  mutate(stable = all(lambda < 0)) %>% 
  filter(spp == 0,
         invasible == 0) %>% 
  group_by(n, d, nfinal) %>% 
  count(name = "nfreq")

#plots of proportion of uninvadable eq of each size.
ggplot(dataplot, aes(x = d, fill = nfinal/n)) +
  geom_bar(stat = "identity", 
           position = "fill",
           aes(y = nfreq))+
  scale_fill_distiller(palette = "YlOrRd",
                       direction = 1)+
  facet_wrap(~n)

#histograms of number of equilibria for each size

dat = read.table("histdresults.csv", sep = ",")
colnames(dat) = c("n", "sim", "d", "nfinal", "subcom", "spp", "xstar", 
                  "lambda","invasible")

dataplothist = dat %>% 
  group_by(n, sim, d, subcom) %>% 
  mutate(stable = all(lambda < 0)) %>% 
  filter(spp == 0,
         invasible == 0) %>% 
  group_by(n, d, nfinal)

ggplot(dataplothist, aes(x = nfinal))+
  geom_bar(aes(fill = as.factor(n)))+
  scale_fill_brewer(palette = "YlOrRd")+
  facet_wrap(~interaction(n, d),
              nrow = length(unique(dataplothist$d)))
