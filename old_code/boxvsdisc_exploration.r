library(tidyverse)
library(geometry)

#load data and name columns
dat_disc = read.table("../data/disc_exploration.csv")
colnames(dat_disc) = c("sim", "alpha", 
                       "dtheta", "theta", "rho", 
                       "r01", "r02",
                       "r1", "r2", 
                       "x1", "x2", 
                       "x01", "x02", 
                       "xl1", "xl2", 
                       "area")

dat_box = read.table("../data/box_exploration.csv")
colnames(dat_box) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                      "cubeind", "dtheta", "drho", "sppid", "xeq")

#select appropriate columns

dat_disc_red = dat_disc %>% 
  select(c(sim, alpha, dtheta, rho, r1, r2, x1, x2))

dat_box_red_feasible = dat_box %>% 
  select(sim, alpha, idrowvec, eqid, rvec, cubeind, dtheta, drho, sppid, xeq) %>% 
  group_by(sim, alpha, idrowvec, eqid) %>% #get only feasible equilibria
  filter(all(xeq>0),
         sqrt(sum(rvec^2)) < 1)

dat_box_red_full_feasible = dat_box %>% 
  group_by(sim, alpha, cubeind) %>% 
  filter(all(xeq>0))

#transform box data into wider
dat_box_red_wide = dat_box_red_full_feasible %>% 
  mutate(rowsid = idrowvec + eqid) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec))
###########################################
#plot overlaying to see if they are the same

#plot the rs

ggplot(dat_disc_red)+
  geom_point(aes(x = r1, y = r2),
             alpha = 0.5, size = 0.5)+
  geom_point(data = dat_box_red_wide, 
             aes(x = rvec_1, 
                 y = rvec_2,),
             color = "red")+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1)


ggplot(dat_disc_red)+
  geom_point(aes(x = x1, y = x2),
             alpha = 0.5, size = 0.5)+
  geom_point(data = dat_box_red_wide, 
             aes(x = xeq_1, 
                 y = xeq_2,),
             color = "red")+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1)


#calculate area inside convex hull
pointsx = dat_box_red_wide %>%ungroup() %>%  filter(alpha == 1) %>% select(xeq_1, xeq_2)

area = convhulln(pointsx, output.options = T)$area