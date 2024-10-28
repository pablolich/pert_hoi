library("tidyverse")

dat = read.table("../data/portrait.csv")
colnames(dat) = c("sim", "n", "alpha", "spp", "eqID", "xRe", "lambdaRe")

datplot = dat %>% group_by(sim, n, eqID, alpha) %>% 
  mutate(feasible = all(xRe > 0),
         stable = lambdaRe<0) %>% 
  filter(feasible == T,
         stable==T)

nsims = max(datplot$sim)
nmax = max(datplot$n)

plotlist = list()
ind = 1
for (ncurr in 1:nmax){
  for (i in 1:nsims){
    data_i = datplot %>% 
      ungroup(alpha) %>% 
      filter(sim == i, 
             n == ncurr)
    if (dim(data_i)[1]==0){
    }
    else{
      ploti= ggplot(data_i, aes(x=alpha, y=xRe, 
                                color = as.factor(spp), 
                                size = stable))+
        geom_point()+
        scale_y_continuous(trans='log10')+
        scale_size_manual(values = c(0.1, 1))+
        theme(aspect.ratio = 0.8)
      plotlist[[ind]] = ploti
      ind = ind + 1
    }
  }
}

# Another option: create pdf where each page is a separate plot.
pdf("../data/portraitsigned.pdf")
for (i in 1:length(plotlist)) {
  print(plotlist[[i]])
}
dev.off()

#questions about curves vs parameters (need to save them)
#do similar phase portraits correspond to similar parametrizations, or to the fact that
#certain sign combinations are dynamically equivalent?

#questions to ask about feasibitliy
#is there any correlation between robustness and number of equlibria?
#what is the alpha regime for which feasibility is most common in each combination? and for what diversity

#questions to ask about stability
#for what combinations is stability most robust to parameter variation?
#does stability correlate with number of equilibria?
#how does stability evolve as I increase n?
#are there situations in which equilibrium are the same, but stability are different?
#are unstable equilibria more robust than stable equilibria?
#what is the magnitude of stable equilibria vs non-stable equilibria?