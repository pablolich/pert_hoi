library("tidyverse")

dat = read.table("../data/solutionevoloution.csv")
colnames(dat) = c("sim", "n", "nfinal", "alpha", "spp", "xRe", "xIm", "lambdaRe")

datplot = dat %>% group_by(sim, n) %>% 
  filter(xRe > 0.01,
         nfinal == 2)%>%
  ggplot(aes(x=alpha, y=xRe, 
             color = as.factor(spp), 
             group = as.factor(spp)))+
  geom_point()+
  scale_y_continuous(trans='log10')+
  facet_wrap(~n)

datplot
