library(tidyverse)

dat = read.table("../data/perturbationdroptractpres.csv")
colnames(dat) = c("n", "alphavalue", "pertmag", "sim", "rpert","xpert", "ratioeqs")

#plot ratio as a function of alpha for different n

datplot = dat %>% 
  group_by(n, alphavalue, pertmag) %>% 
  #filter(pertmag == 0.01) %>% 
  summarise(avratio = mean(ratioeqs),
            avratiotheo = mean(xpert))

#Numerical simulations of perturbation respnse as a function of alpha
ggplot(datplot, aes(x = alphavalue, 
                    y = avratio,
                    color = as.factor(n))) +
  geom_point()+
  scale_y_continuous(trans = "log10")+
  theme(aspect.ratio = 1)+
  facet_wrap(~pertmag)

#Agreement with theory as a function of the perturbation
ggplot(datplot, aes(x = avratiotheo/pertmag, 
                    y = avratio/pertmag,
                    color = as.factor(pertmag)))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme(aspect.ratio = 1)


  