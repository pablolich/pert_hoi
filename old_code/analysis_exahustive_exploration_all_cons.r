library(tidyverse)

dat_exhaustive = read.table("../data/exhaustive_exploration_all_cons.csv")
colnames(dat_exhaustive) = c("sim", "n", "alpha", "idrow", "eqid", "rvec", "dtheta", "rho", "sppid", "xeq")

dat_exhaustive_feasible = dat_exhaustive %>% 
  mutate(xeqround = round(xeq, digits = 4)) %>% 
  group_by(sim, alpha, idrow, eqid) %>% #get only feasible equilibrium
  filter(all(xeq>0)) %>% 
  mutate(deltaxmag = sqrt(sum((1-xeq)^2)),
         relresp = deltaxmag/rho) %>% 
  ungroup() %>% 
  distinct(sim, n, alpha, idrow, sppid, xeqround, .keep_all = T) %>% 
  select(!xeqround) %>% 
  group_by(sim, alpha, idrow, sppid) %>%
  mutate(neq = n()) %>%
  ungroup() %>%
  filter(neq == 1) #get only situations leading to one equilibrium

dat_exhaustive_feasible_wide = dat_exhaustive_feasible %>% 
  mutate(rowsid = idrow+eqid) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec))  #make wider

dat_magresponse = dat_exhaustive_feasible_wide %>% 
  filter(dtheta < 1e-1,
         rho < 0.11) %>% 
  group_by(n, alpha) %>% 
  summarize(relrespav = median(relresp))

ggplot(dat_magresponse, 
       aes(x = as.factor(alpha), y = relrespav))+
  geom_point()
  # facet_wrap(~rho,
  #            scales = "free")

ggplot(dat_exhaustive_feasible_wide,
       aes(x = rvec_1,
           y = rvec_2)) +
  geom_point()+
  facet_wrap(~alpha)


ggplot(dat_exhaustive_feasible_wide,
       aes(x = xeq_1,
           y = xeq_2)) +
  geom_point()+
  facet_wrap(~alpha)

