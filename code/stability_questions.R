library(tidyverse)

#load data
dat = read.table("../data/portraitsigned.csv")
colnames(dat) = c("sim", "n", "rsign", "Asign", "Bsign", "alpha", "spp", "eqID",
                  "xRe", "lambdaRe")

datplot = dat %>% group_by(sim, n, eqID, alpha, rsign, Asign, Bsign) %>% 
  mutate(stable = lambdaRe<0,
         feasible = all(xRe > 0))

nmax = max(datplot$n)
nsim = max(datplot$sim)

###############################################################################

#1. Number of stable equilibria?

dataq1 = datplot %>%
  ungroup() %>% 
  group_by(sim, n, alpha, rsign, Asign, Bsign, stable) %>% 
  filter(feasible == T,
         spp == 1, #to make sure I only count each community once
         n ==2) %>% 
  tally(name = "neq") %>%
  group_by(n, alpha, rsign, Asign, Bsign, stable) %>% 
  summarize(neqav = mean(neq))


ggplot(dataq1)+
  geom_point(aes(x = alpha, y = neqav,
                 color = stable,
                 shape = stable))+
  scale_color_manual(values = c("gray", "orange"))+
  scale_shape_manual(values = c(20,3))+
  facet_wrap(~interaction(rsign, Asign, Bsign))

###############################################################################
# Stability robustness
#to quantify the domain, do all possible domains and plot in gray
#then overlay the 27 points of interest

#allwindows = function(alphavec)
#loop through window size
alphavec = seq(0,1,0.05)
pointsnull = data.frame()
for (wsize in 1:length(alphavec)){
  #loop through window position
  #get how many possible positions there are.
  npos = length(alphavec)-wsize
  if (npos == 0){
    pointsnull = rbind(pointsnull, c(mean(alphavec), sd(alphavec)))
  }
  else{
    for (wpos in 0:npos){
      window = seq(1,wsize) + wpos
      alphawind = alphavec[window]
      pointostore = c(mean(alphawind), sd(alphawind))
      if (is.na(pointostore[2])){
        pointostore[2] = 0
      }
      pointsnull = rbind(pointsnull, pointostore)
    }
  }
}
colnames(pointsnull) = c("alphaav", "alphasd")

#get the empirical scatter plot for the alpha's distribution
#get each parametrization
#get each simulation at that parametrization
#calculate the average alpha, and the standard deviation

dataq2 = datplot %>% 
  filter(spp == 1,
         feasible == T,
         n==4) %>% 
  group_by(sim, n, rsign, Asign, Bsign, stable) %>%
  summarise(alphaav = mean(alpha),
            alphasd = sd(alpha))

empiricalplot = ggplot()+
  geom_point(data = pointsnull,
             aes(x = alphaav,
                 y = alphasd),
             alpha = 0.1)+
  geom_point(data = dataq2, 
             aes(x = alphaav,
                 y = alphasd,
                 color = stable,
                 shape = stable),
             alpha = 0.5,
             size = 1.5)+
  scale_color_manual(values = c("black", "red"))+
  scale_shape_manual(values = c(3, 20))+
  facet_wrap(~interaction(rsign, Asign, Bsign))


empiricalplot  

###############################################################################

#what is the magnitude of equilibria on average, and as a function of 
#parametrizations, alpha, species how does it increase/decrease as alpha moves?

dataq3 = datplot %>% 
  filter(n == 4, 
         feasible == T) %>% 
  group_by(n, alpha, rsign, Asign, Bsign, stable) %>% 
  summarise(xav = mean(xRe))

plotq3 = ggplot(dataq3) +
  geom_point(aes(x = alpha,
                 y = xav,
                 color = stable,
                 shape = stable)) +
  scale_color_manual(values = c("gray", "orange"))+
  scale_shape_manual(values = c(20,3))+
  scale_y_continuous(trans='log10')+
  facet_wrap(~interaction(rsign, Asign, Bsign))

plotq3
