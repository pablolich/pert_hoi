library(tidyverse)

###############################################################################
#FIGURE 2
###############################################################################

#load data
datarbitraryperts = read.table("../data/simulations_arbitrary_perturbationssim43.csv")
colnames(datarbitraryperts) = c("n", "sim", "rho",
                                "alpha", "pertid", "sppid",
                                "r0", "rpert", "dtheta",
                                "xstar", "xstarlin")

######## ######## ######## ######## PANEL 1 ######## ######## ######## ########
###############################################################################

datasmallorbs = datarbitraryperts %>% 
  mutate(deltar = r0 - rpert) %>% 
  select(n, sim, rho, alpha, pertid, sppid, deltar, xstar) %>% 
  filter(n==2,
         sim == 43, 
         rho == 0.11,
         alpha == 0.01 | alpha == 0.61 | alpha == 0.81
  ) %>% 
  group_by(n, sim, rho, alpha, sppid, pertid) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(deltar, xstar))

datasmallorbsrvectex =
  datasmallorbs %>% 
  ungroup() %>% 
  group_by(alpha) %>% 
  mutate(set = cur_group_id()) %>% 
  ungroup() %>% 
  select(deltar_1, deltar_2, set) %>% 
  rename(x = deltar_1,
         y = deltar_2)

write.table(datasmallorbsrvectex, 
            file = "../data/smallorbsr.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

ggplot(datasmallorbsrvectex, aes(x, y, color = set))+
  geom_point()

######## ######## ######## ######## PANEL 2 ######## ######## ######## ########
###############################################################################

datasmallorbs = datarbitraryperts %>% 
  mutate(deltar = r0 - rpert) %>% 
  select(n, sim, rho, alpha, pertid, sppid, deltar, xstar) %>% 
  filter(n==2,
         sim == 43, 
         rho == 0.11,
         alpha == 0.01 | alpha == 0.61 | alpha == 0.81
  ) %>% 
  group_by(n, sim, rho, alpha, sppid, pertid) %>% 
  pivot_wider(names_from = sppid,
            values_from = c(deltar, xstar)) 

datasmallorbstex =
  datasmallorbs %>% 
  ungroup() %>% 
  group_by(alpha) %>% 
  mutate(set = cur_group_id()) %>% 
  ungroup() %>% 
  mutate(deltax1 =1-xstar_1,
         deltax2 =1-xstar_2) %>% 
  select(deltax1, deltax2, set) %>% 
  rename(x = deltax1,
         y = deltax2)

write.table(datasmallorbstex, 
            file = "../data/smallorbsx.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)
#check
ggplot(datasmallorbstex, aes(x, y, color = set))+
         geom_point()


######## ######## ######## ######## PANEL 3 ######## ######## ######## ########
###############################################################################



dataavradius = datarbitraryperts %>%
  filter(n == 2,
         rho == 0.11) %>% 
  mutate(deltar = rpert-r0) %>% 
  group_by(n, sim, rho, alpha, pertid) %>% 
  mutate(magx = sum((xstar-1)^2),
         magr = sum(deltar^2)) %>% 
  ungroup() %>% #calculate average distance for each n
  group_by(n, alpha) %>% 
  filter(dtheta < 1e-2) %>% 
  ungroup() %>% 
  group_by(n, alpha) %>% 
  summarize(avradius = mean(magx)) %>% 
  ungroup() %>% 
  rename(x = alpha, 
         y = avradius,
         set = n) %>% 
  relocate(set, .after = y) %>% 
  mutate(set = case_when(x == 0.01 ~ 1, 
                         x == 0.61 ~ 2, 
                         x == 0.81 ~ 3)) %>% 
  replace_na(list(set = 4))  

write.table(dataavradius, 
            file = "../data/decreasingradius.dat", sep = "\t", 
            row.names = F,
            quote = FALSE)

ggplot(dataavradius,
       aes(x , 
           y ))+
  geom_point(aes(color= set))


