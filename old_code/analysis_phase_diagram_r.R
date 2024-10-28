library(tidyverse)
library(gridExtra)
library(RColorBrewer)
#library(tune)

dat = read.table("../data/feasibility_allcons_niceviz.csv")
colnames(dat) = c("sim", "alpha", 
                  "dtheta", "theta", "rho", 
                  "r01", "r02",
                  "r1", "r2",
                  "x1", "x2",
                  "x01", "x02",
                  "x1l", "x2l",
                  "area")

dat = dat %>% 
  mutate(xmag = sqrt(((x1/x01-1))^2 + ((x2/x02-1))^2),
         xmagnorm = xmag/rho,
         xmaglin = sqrt((x1l/x01-1)^2 + (x2l/x02-1)^2),
         xmaglinnorm = xmaglin/rho,
         homo = dtheta<1e-6,
         orbit = rho > 0 & rho < 2e-1)

pdf("../data/phaseportraitallcons.pdf")
#less alphas for plotting purposes
datred = dat %>% filter(alpha == 0.4| alpha == 0.7 | alpha == 0.8)
for (i in 1:max(dat$sim)){
  datsim = datred %>% filter(sim == i, 
                             rho < 6e-1)
  #datsim = dat %>% filter(sim == i)
  
  p =
    ggplot(datsim, aes(x = x1/(x01), y = x2/(x02),
                     color = orbit                     ))+
    geom_point()+
    geom_point(data = datsim,
               aes(x = x1l/(x01), y = x2l/(x02)),
               color = "gray50",
               size = 0.3)+
    scale_color_manual(values = c("black", "blue"))+
    scale_size_manual(values = c(1, 2))+
    geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
    #coord_cartesian(xlim = c(0, NA), ylim = c(0,NA))+
    facet_wrap(~alpha)+
    #scale_color_viridis_c(option = "plasma")+
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.position = "none")
  ggsave("../data/example_manuscript.pdf", width = 6, height = 2.5)
  q =
    ggplot(datsim,
             aes(r1, r2,
                 color = xmag))+
    geom_point()+
    #labs(title = paste0(i))+
    geom_point(data = dataround,
               aes(x = r1, y = r2),
               color = "green")+
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
    scale_color_viridis_c(option = "plasma")+
    facet_wrap(~alpha)+
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill = NA),
          legend.position = "none")
  print(p)
  print(q)
  #grid.arrange(q, p, ncol= 2)
}
dev.off()

 ###############################################################################
# dataxmin = dat %>% 
#    group_by(alpha, sim, rho) %>% 
#    slice_min(xmag) %>% 
# #plot of deltatheta versus alpha for the minimum distances
#    ungroup() %>% 
#    group_by(alpha) %>% 
#    summarise(dtheta_av = median(dtheta))
# 
#  ggplot(dataxmin, aes(x = alpha, y=dtheta_av)+
#    geom_point()
 
 #plot the distances in x for the cone
 
 datcone = dat %>% filter(dtheta > 0, 
                          rho == 0.1) %>% 
   group_by(alpha, rho) %>% 
   summarize(avmag = median(xmag),
             avmaglin = median(xmaglin)) %>% 
  filter(avmaglin<10)
 
 ggplot(datcone, aes(x = alpha, y = avmag,
                     #color = dtheta
                     ))+
   geom_point()+
   facet_wrap(~rho,
              scales = "free") %>% 
   theme(aspect.ratio = 1,
         panel.grid = element_blank(),
         panel.background = element_blank(),
         panel.border = element_rect(fill = NA))

 ###############################################################################
 
 
 #when there are no constraints, do the change of variables to get the
 #equilibrium to 1 and look at the phase space on the new variables.
 
 
