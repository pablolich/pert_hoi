library(tidyverse)

###############################################################################
#FIGURE 3
###############################################################################

#load data
datarbitrarypertssim43 = read.table("../data/simulations_arbitrary_perturbationssim1.csv")
colnames(datarbitrarypertssim43) = c("n", "sim", "rho", "ismax",
                                     "alpha", "pertid", "sppid",
                                     "r0", "rpert", "dtheta", "xstar", "xstarlin")

datarbitrarypertssim43 = datarbitrarypertssim43 %>% mutate(deltaxstar = xstar - 1,
                                                           deltaxstarlin = xstarlin - 1)
######## ######## ########  PANELS 1, 2, 3 ######## ######## ########
####################################################################

datar = datarbitrarypertssim43 %>% 
  mutate(deltar = r0 - rpert) %>% 
  select(n, sim, rho, alpha, pertid, sppid, deltar, xstar, xstarlin) %>% 
  group_by(n, sim, rho, alpha, sppid, pertid) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(deltar, xstar, xstarlin)) %>% 
  select(n, sim, rho, alpha, pertid, deltar_1, deltar_2)

dataralpha1 = datar %>% 
  filter(alpha == 0.01) %>% 
  mutate(set = if_else(rho == 0.11, 1, 2)) %>% 
  ungroup() %>% 
  select(deltar_1, deltar_2, set) %>% 
  rename(x = deltar_1,
         y = deltar_2)

ggplot(dataralpha1,
       aes(x = x, 
           y = y))+
  geom_point(size = 0.5)+
  theme(aspect.ratio = 1)

dataralpha2 = datar %>% 
  filter(alpha == 0.61) %>% 
  mutate(set = if_else(rho == 0.11, 1, 2)) %>%
  ungroup() %>% 
  select(deltar_1, deltar_2, set) %>% 
  rename(x = deltar_1,
         y = deltar_2)

dataralpha3 = datar %>% 
  filter(alpha == 0.81) %>% 
  mutate(set = if_else(rho == 0.11, 1, 2)) %>% 
  ungroup() %>% 
  select(deltar_1, deltar_2, set) %>% 
  rename(x = deltar_1,
         y = deltar_2)

write.table(dataralpha1, 
            file = "../data/dataralpha1.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

write.table(dataralpha2, 
            file = "../data/dataralpha2.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

write.table(dataralpha3, 
            file = "../data/dataralpha3.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

######## ######## ########  PANELS 4, 5, 6 ######## ######## ########
#####################################################################


datax = datarbitrarypertssim43 %>% 
  mutate(deltar = r0 - rpert) %>% 
  select(n, sim, rho, alpha, pertid, sppid, deltar, xstar, xstarlin) %>% 
  group_by(n, sim
           , rho, alpha, sppid, pertid) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(deltar, xstar, xstarlin)) %>% 
  select(n, sim, rho, alpha, pertid, xstar_1, xstar_2) %>% 
  rename(xstar1 = xstar_1,
         xstar2 = xstar_2) %>% 
  mutate(type = "num")

dataxlin = datarbitrarypertssim43 %>% 
  mutate(deltar = r0 - rpert) %>% 
  select(n, sim, rho, alpha, pertid, sppid, deltar, xstar, xstarlin) %>% 
  group_by(n, sim, rho, alpha, sppid, pertid) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(deltar, xstar, xstarlin)) %>% 
  select(n, sim, rho, alpha, pertid, xstarlin_1, xstarlin_2) %>% 
  rename(xstar1 = xstarlin_1,
         xstar2 = xstarlin_2) %>% 
  mutate(type = "lin")

datamerged = rbind(datax, dataxlin)

dataxalpha1 = datamerged %>% filter(alpha == 0.01) %>% 
  group_by(type) %>% 
  mutate(set = case_when((type == "num" & rho == 0.11) ~ 1,
                         (type == "num" & rho != 0.011) ~ 2,
                         type == "lin" ~ 3)) %>% 
  ungroup() %>% 
  select(xstar1, xstar2, set) %>% 
  rename(x = xstar1,
         y = xstar2)

ggplot(dataxalpha1, aes(x, y))+
  geom_point(aes(color = set))+
  theme(aspect.ratio = 1)

dataxalpha2 = datamerged %>% filter(alpha == 0.61) %>% 
  group_by(type) %>% 
  mutate(set = case_when((type == "num" & rho == 0.11) ~ 1,
                         (type == "num" & rho != 0.011) ~ 2,
                         type == "lin" ~ 3)) %>% 
  ungroup() %>% 
  select(xstar1, xstar2, set) %>% 
  rename(x = xstar1,
         y = xstar2)

dataxalpha3 = datamerged %>% filter(alpha == 0.81) %>% 
  group_by(type) %>% 
  mutate(set = cur_group_id()) %>% 
  mutate(set = case_when((type == "num" & rho == 0.11) ~ 1,
                         (type == "num" & rho != 0.011) ~ 2,
                         type == "lin" ~ 3)) %>% 
  ungroup() %>% 
  select(xstar1, xstar2, set) %>% 
  rename(x = xstar1,
         y = xstar2)

write.table(dataxalpha1, 
            file = "../data/dataxalpha1.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

write.table(dataxalpha2, 
            file = "../data/dataxalpha2.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

write.table(dataxalpha3, 
            file = "../data/dataxalpha3.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)


