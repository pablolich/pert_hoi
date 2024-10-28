library("tidyverse")

dat2sppgen = read.table("../data/simulation2sppgen.csv")
colnames(dat2sppgen) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                         "cubeind", "dtheta", "drho", "sppid", "xeq")

datagenwide = dat2sppgen %>%
  mutate(rowsid = idrowvec + eqid) %>% 
  mutate(bin = cut(drho, seq(min(drho), max(drho), (max(drho)-min(drho))/11), right = FALSE)) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec)) %>% 
  filter(xeq_1 > 0 & xeq_2 > 0) %>% 
  mutate(xmag = sqrt(xeq_1^2 + xeq_2^2)) %>% 
  filter(sim == 16)

#plot orbit
ggplot(datagenwide, aes(x = xeq_1, y = xeq_2))+
  geom_point()+
    theme(aspect.ratio = 1)+
  facet_wrap(~alpha)

datamags = datagenwide %>% 
  group_by(alpha, bin) %>%
  filter(dtheta < 1e-1) %>% 
  summarize(avmag = median(xmag))

ggplot(datamags, aes(x = alpha, y = avmag))+
  geom_point()+
  facet_wrap(~bin)

###############################################################################

dat2spp = read.table("../data/simulations2spp.csv")

colnames(dat2spp) = c("sim", "alpha", 
                      "dtheta", "theta", "rho", 
                      "r01", "r02",
                      "r1", "r2",
                      "x1", "x2",
                      "x01", "x02",
                      "x1l", "x2l",
                      "area")

dat2spp = dat2spp %>%
  filter(x1 > 0 & x2 > 0) %>% 
  mutate(xmag = sqrt(x1^2 + x2^2))

datamags2spp = dat2spp %>% 
  filter(dtheta < 1e-1) %>% 
  group_by(alpha, rho) %>% 
  summarize(avmag = median(xmag))


ggplot(datamags2spp, aes(x = alpha, y = avmag))+
  geom_point()+
  facet_wrap(~rho)
