library(tidyverse)

dat = read.table("../data/simulations2spp.csv")
colnames(dat) = c("sim", "alpha", 
                  "dtheta", "theta", "rho", 
                  "r01", "r02",
                  "r1", "r2",
                  "x1", "x2",
                  "x01", "x02",
                  "x1l", "x2l",
                  "area")

datorbav = dat %>%
  group_by(alpha, theta, rho) %>% 
  summarize(x1av = median(x1),
            x2av = median(x2)) %>% 
  filter(rho == 0.5)

datxmagav = dat %>% 
  mutate(xmag = sqrt(x1^2 + x2^2)) %>% 
  group_by(alpha, rho) %>% 
  summarise(xmagav= median(xmag))

ggplot(datorbav, aes(x = x1av, y = x2av))+
  geom_point()+
  facet_wrap(~alpha,
             scales = "free")

ggplot(datxmagav, aes(x = alpha, y = xmagav))+
  geom_point()+
  facet_wrap(~rho,
             scales = "free")

#FOR THIS FIGURE TO BE BETTER i SHOULD RUN MORE SIMULATIONS


#check also the better exploration case;

datgen = read.table("../data/simulation2sppgen.csv")
colnames(datgen) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                         "dtheta", "drho", "sppid", "xeq")

datagenwide = datgen %>%
  mutate(rowsid = idrowvec + eqid) %>% 
  mutate(binrho = cut(drho, seq(min(drho), max(drho), (max(drho)-min(drho))/20), 
                      right = FALSE),
         bintheta = cut(dtheta, seq(0, pi, 0.1), right = FALSE)) %>% 
  group_by(sim, alpha, binrho, bintheta) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec)) %>% 
  filter(xeq_1 > 0 & xeq_2 > 0) %>% 
  ungroup() %>% 
  group_by(alpha, binrho, bintheta) %>% 
  summarise(x1av = mean(xeq_1),
            x2av = mean(xeq_2)) %>% 
  filter(x1av < 10 & x2av < 10)

ggplot(datagenwide, aes(x = x1av, y = x2av))+
  geom_point()+
  facet_wrap(~alpha,
             scales = "free")
