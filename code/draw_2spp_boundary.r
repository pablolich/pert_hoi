library(tidyverse)
library(ggforce)


alldata = read.table("../data/results_2spp_boundaries.csv")
data_rhocrit = read.table("../data/results_2spp_rho_crits.csv")
colnames(alldata) = c("sim", "alpha", "r1", "r2")
colnames(data_rhocrit) = c("sim", "alpha", "rhodense", "rhosparse")

n_sims = max(alldata$sim)

pdf("../data/boundary_plots.pdf", width = 8, height = 6.5)
for (sim_i in 1:n_sims){
  data_sim = alldata %>% filter(sim == sim_i)
  data_rho_sim = data_rhocrit %>% filter(sim == sim_i) %>% 
    pivot_longer(c(rhodense, rhosparse),
                 names_to = "type",
                 values_to = "rho")
  data_rho_sim$x0 = rep(0, nrow(data_rho_sim))
  data_rho_sim$y0 = rep(0, nrow(data_rho_sim))
  p <- 
    ggplot(data = data_sim)+
    geom_point(aes(x = r1, y = r2,
                   color = as.factor(alpha)),
               size = 0.7)+
    geom_circle(data = data_rho_sim,
                aes(x0 = 0, y0 = 0, r = rho,
                    linetype = as.factor(type),
                    #color = as.factor(alpha)
                    )
                )+
    scale_color_brewer(palette = "YlOrRd")+
    xlim(-2, 2)+
    ylim(-2, 2)+
    facet_wrap(~alpha)+
    theme(aspect.ratio = 1)
  print(p)
}
dev.off()