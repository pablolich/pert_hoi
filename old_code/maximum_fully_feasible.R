library(tidyverse)

polygon_area <- function(vertices) {
  n <- nrow(vertices)
  if (n < 3) {
    return(0.0)  # Not a polygon
  }
  
  area <- 0.0
  for (i in 1:n) {
    x1 <- vertices[i, 1]
    y1 <- vertices[i, 2]
    x2 <- vertices[(i %% n) + 1, 1]  # Next vertex, with wrap-around
    y2 <- vertices[(i %% n) + 1, 2]
    area <- area + x1 * y2 - y1 * x2
  }
  
  area <- abs(area) / 2.0
  return(area)
}

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

pdf("../data/maximum_full_feasible_letgo.pdf")
nsimvec = unique(dat_box_red_wide$sim)
for (i in nsimvec){
  datsim= dat_box_red_wide %>%
    filter(sim == i)

  p = ggplot(datsim)+
    geom_point(aes(x = xeq_1, y = xeq_2),
               alpha = 0.5, size = 0.5)+
    facet_wrap(~alpha)+
    theme(aspect.ratio = 1)
  print(p)
}
dev.off()

#average area
dat_area = dat_box_red_wide %>% 
  ungroup() %>% 
  group_by(alpha) %>% 
  summarise(av_area = mean(area))

ggplot(dat_area) +
  geom_point(aes(x=alpha, y = av_area))

