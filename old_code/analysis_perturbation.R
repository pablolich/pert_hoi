library(tidyverse)

dat = read.table("../data/og_polar_pert.csv")
colnames(dat) = c("n", "alphavalue", "pert", "rpert","xpert", "ratioeqs")
plot(dat$alphavalue, dat$ratioeqs)
#plot ratio as a function of alpha for different n

datplot = dat %>% 
  group_by(alphavalue) %>% 
  summarise(avratio = mean(ratioeqs),
            avratiotheo = mean(xpert))

#Numerical simulations of perturbation respnse as a function of alpha
ggplot(datplot, aes(x = alphavalue, 
                    y = avratio)) +
  geom_point()+
  scale_y_continuous(trans = "log10")+
  theme(aspect.ratio = 1)

#Agreement with theory as a function of the perturbation
ggplot(datplot, aes(x = avratiotheo,
                    y = avratio))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  theme(aspect.ratio = 1)


  