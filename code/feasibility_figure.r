library(tidyverse)
library(gridExtra)


datmax = read.table("../data/maximum_fully_feasible43.csv",
                    nrow = 1e6)
colnames(datmax) = c("sim", "n", "alpha", "idrow", "rho", "eqid", "rvec", "sppid", "xeq", "npert")

 datmax43 = datmax %>% filter(alpha == 0.8)

dat_fully_feas = datmax43%>% 
  group_by(sim, n, alpha, rho) %>% 
  mutate(full = any(npert==200)) %>% 
  filter(full == T) %>% 
  ungroup() %>% 
  group_by(sim, alpha, idrow, eqid) %>% #get only feasible equilibrium
  filter(all(xeq>0),
         all(xeq<4)
  ) %>%
  ungroup() %>%
  group_by(sim, alpha, idrow, sppid) %>%
  mutate(neq = n()) %>%
  ungroup() %>%
  group_by(sim, alpha, idrow) #%>%
  #filter(neq == 1)

dat_box_red_wide = dat_fully_feas %>% 
  mutate(rowsid = idrow+eqid) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec))  %>% 
  ungroup() %>% 
  group_by(sim, n, alpha, eqid) %>% 
  mutate(area = volume_hypersphere(n, rho))



ggplot(dat_box_red_wide)+
  geom_point(aes(x = xeq_1, y = xeq_2),
             alpha = 0.5, size = 0.5)+
  facet_wrap(~alpha,
             scales = "free")+
  theme(aspect.ratio = 1)

#above I got the orbits I wanted for the maximum feasible orbit

#now I need the exploration of simulation 43

dat_exploration = read.table("../data/feasibility_allcons_niceviz43.csv")
colnames(dat_exploration) = c("sim", "alpha", 
                                              "dtheta", "theta", "rho", 
                                              "r01", "r02",
                                              "r1", "r2",
                                              "x1", "x2",
                                              "x01", "x02",
                                              "x1l", "x2l",
                                              "area")
dat_exploration43 = dat_exploration %>% 
  mutate(rhoap = round(rho, digits = 2),
         orbit = rho == 0.2) %>% 
  filter(
         alpha == 0.8,
       rhoap == 0 |
         rhoap == 0.2|
         rhoap == 0.4 |
         rhoap == 0.6 |
         rhoap == 0.8 |
         rhoap == 1.2       )

#now get the orbits
rsplot = ggplot(dat_exploration43, aes(x = r1, y = r2,
                              color = orbit))+
  geom_point(size = 1.2)+
  geom_point(data = dat_box_red_wide,
             aes(x = rvec_1, y = rvec_2),
             color = "red",
             size = 1.2)+
  scale_color_manual(values = c("black", "blue"))+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

xplot = ggplot(dat_exploration43, aes(x = x1, y = x2,
                              color = orbit))+
  geom_point(size = 1.2)+
  scale_color_manual(values = c("black", "blue"))+
  geom_point(data = dat_box_red_wide,
             aes(x = xeq_1, y = xeq_2),
             color = "red",
            size = 1.2)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

volume_hypersphere <- function(d, r) {
  pi_value <- pi
  gamma_value <- gamma(d/2 + 1)
  volume <- pi_value^(d/2) / gamma_value * r^d
  return(volume)
}


datmax = read.table("../data/maximum_fully_feasible.csv" ,
                    nrows = 20e6)
colnames(datmax) = c("sim", "n", "alpha", "idrow", "rho", "eqid", "rvec", "sppid", "xeq", "npert")

dat_fully_feas = datmax %>% 
  group_by(sim, n, alpha, rho) %>% 
  mutate(full = any(npert==100)) %>% 
  filter(full == T) %>% 
  ungroup() %>% 
  group_by(sim, alpha, idrow, eqid) %>% #get only feasible equilibrium
  filter(all(xeq>0),
         all(xeq<10)
  ) %>%
  ungroup() %>%
  group_by(sim, alpha, idrow, sppid) %>%
  mutate(neq = n()) %>%
  ungroup() %>%
  group_by(sim, alpha, idrow) %>%
  filter(neq == 1)

dat_box_red_wide = dat_fully_feas %>% 
  mutate(rowsid = idrow+eqid) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec))  %>% 
  ungroup() %>% 
  group_by(sim, n, alpha, eqid) %>% 
  mutate(area = volume_hypersphere(n, rho))

#average area
dat_area = dat_box_red_wide %>% 
  ungroup() %>% 
  group_by(alpha) %>% 
  summarise(av_area = mean(area))

feasareaplot = ggplot(dat_area) +
  geom_point(aes(x=alpha, y = av_area),
             size = 1.2)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

pdf("../data/feasibilityareaplot.pdf", width = 6.2, height = 2
)
p = grid.arrange(rsplot, xplot, feasareaplot, nrow =1)
print(p)
dev.off()
ggsave("../data/feasibilityareaplot.pdf", width = 6.2, height = 2.5)

#example sims figure

dat_exploration43_example = dat_exploration %>% 
  mutate(rhoap = round(rho, digits = 2),
         orbit = rhoap == 0.2) %>% 
  filter(alpha == 0.6 |
           alpha == 0.8,
         rhoap == 0 |
           rhoap == 0.2|
           rhoap == 0.4 |
           rhoap == 0.6 |
           rhoap == 0.8 |
           rhoap == 1 |
           rhoap == 1.2
  )

alpha06 = ggplot(dat_exploration43_example %>% filter(alpha == 0.6), aes(x = x1/(x01), y = x2/(x02),
                   color = orbit))+
  geom_point(size = 1.2)+
  geom_point(data = dat_exploration43_example %>% filter(alpha == 0.6),
             aes(x = x1l/(x01), y = x2l/(x02)),
             color = "gray50",
             size = 0.2)+
  scale_color_manual(values = c("black", "blue"))+
  scale_size_manual(values = c(1, 2))+
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
  coord_cartesian(xlim = c(0, NA), ylim = c(0,3.25))+
  #facet_wrap(~alpha)+
  #scale_color_viridis_c(option = "plasma")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

alpha08 = ggplot(dat_exploration43_example %>% filter(alpha == 0.8), aes(x = x1/(x01), y = x2/(x02),
                                                                         color = orbit))+
  geom_point(size = 1.2)+
  geom_point(data = dat_exploration43_example %>% filter(alpha == 0.8),
             aes(x = x1l/(x01), y = x2l/(x02)),
             color = "gray50",
             size = 0.2)+
  scale_color_manual(values = c("black", "blue"))+
  scale_size_manual(values = c(1, 2))+
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
  coord_cartesian(xlim = c(0, NA), ylim = c(0,3.25))+
  #facet_wrap(~alpha)+
  #scale_color_viridis_c(option = "plasma")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")


#decreasing xmag for theta = 0
dat_exploration_proc = dat_exploration %>% 
  mutate(xmag = sqrt(((x1/x01-1))^2 + ((x2/x02-1))^2),
         xmagnorm = xmag/rho,
         xmaglin = sqrt((x1l/x01-1)^2 + (x2l/x02-1)^2),
         xmaglinnorm = xmaglin/rho,
         homo = dtheta<1e-6,
         orbit = rho > 0 & rho < 2e-1)

datcone = dat_exploration_proc %>% filter(dtheta > 0, 
                                          rho == 0.1) %>% 
  group_by(alpha, rho) %>% 
  summarize(avmag = median(xmag),
            avmaglin = median(xmaglin)) %>% 
  filter(avmaglin<10)

increasingxmagplot = ggplot(datcone, aes(x = alpha, y = avmag,
                                         #color = dtheta
))+
  geom_point()+
  facet_wrap(~rho,
             scales = "free") %>% 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))


pdf("../data/largenonhomogeneous_example.pdf", width = 6.2, height = 2.5)
p = grid.arrange(alpha06, alpha08, increasingxmagplot, nrow = 1)
print(p)
dev.off()

#####################################################################

#initial example figure


dat_exploration43_initial_example = dat_exploration %>% 
  mutate(rhoap = round(rho, digits = 2),
         orbit = rhoap == 0.2) %>% 
  filter(alpha == 0.6 |
           alpha == 0.8,
         rhoap == 0 |
           rhoap == 0.2
  )

rinitialex = ggplot(dat_exploration43_initial_example, aes(x = r1, y = r2))+
  geom_point(size = 1.2)+
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

xinitialex = ggplot(dat_exploration43_initial_example, aes(x = x1/(x01), y = x2/(x02),
                                      color = orbit))+
  geom_point(size = 1.2)+
  geom_point(data = dat_exploration43_initial_example,
             aes(x = x1l/(x01), y = x2l/(x02)),
             color = "gray50",
             size = 1.2, 
             alpha = 0.5)+
  scale_color_manual(values = c("black", "blue"))+
  scale_size_manual(values = c(1, 2))+
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
  #coord_cartesian(xlim = c(0, NA), ylim = c(0,NA))+
  #scale_color_viridis_c(option = "plasma")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")


#decreasing xmag for theta = 0
dat_exploration_proc = dat_exploration %>% 
  mutate(xmag = sqrt(((x1/x01-1))^2 + ((x2/x02-1))^2),
         xmagnorm = xmag/rho,
         xmaglin = sqrt((x1l/x01-1)^2 + (x2l/x02-1)^2),
         xmaglinnorm = xmaglin/rho,
         homo = dtheta<1e-6,
         orbit = rho > 0 & rho < 2e-1)

datcone = dat_exploration_proc %>% filter(dtheta == 0, 
                                          rho == 0.1) %>% 
  group_by(alpha, rho) %>% 
  summarize(avmag = median(xmag),
            avmaglin = median(xmaglin)) %>% 
  filter(avmaglin<10)

decreasingxmagplot = ggplot(datcone, aes(x = alpha, y = avmag,
                                         #color = dtheta
))+
  geom_point(size = 1.2)+
  facet_wrap(~rho,
             scales = "free") %>% 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))



pdf("../data/initial_example.pdf", width = 6.2, height = 2.5)
p = grid.arrange(rinitialex, xinitialex, decreasingxmagplot, nrow = 1)
print(p)
dev.off()


 #####################################################################

#decreasing xmag for theta = 0
dat_exploration_proc = dat_exploration %>% 
  mutate(xmag = sqrt(((x1/x01-1))^2 + ((x2/x02-1))^2),
         xmagnorm = xmag/rho,
         xmaglin = sqrt((x1l/x01-1)^2 + (x2l/x02-1)^2),
         xmaglinnorm = xmaglin/rho,
         homo = dtheta<1e-6,
         orbit = rho > 0 & rho < 2e-1)

datcone = dat_exploration_proc %>% filter(dtheta > 0, 
                         rho == 0.1) %>% 
  group_by(alpha, rho) %>% 
  summarize(avmag = median(xmag),
            avmaglin = median(xmaglin)) %>% 
  filter(avmaglin<10)

decreasingxmagplot = ggplot(datcone, aes(x = alpha, y = avmag,
                    #color = dtheta
))+
  geom_point()+
  facet_wrap(~rho,
             scales = "free") %>% 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

#####################################################################

#plot of feasibility versus reactivity

feasvsreac = tibble(reac = datcone$avmag, 
                    feas = dat_area$av_area)

ggplot(feasvsreac, 
       aes(x = reac, 
           y = feas))+
  geom_point()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA))

ggsave("../data/feasvsreac.pdf", width = 2, height = 2.5
)  
