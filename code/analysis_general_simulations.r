library(tidyverse)


dat2sppgen = read.table("../data/simulation2sppgen.csv")
colnames(dat2sppgen) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                         "dtheta", "drho", "sppid", "xeq")

datagen = dat2sppgen %>%
  mutate(rowsid = idrowvec + eqid) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec)) %>% 
  filter(sim == 2,
         xeq_1 > 0 & xeq_2 > 0)
 

############################################################################

dat2spp = read.table("../data/simulations2spp.csv")
colnames(dat2spp) = c("sim", "alpha", 
                  "dtheta", "theta", "rho", 
                  "r01", "r02",
                  "r1", "r2",
                  "x1", "x2",
                  "x01", "x02",
                  "x1l", "x2l",
                  "area")

dat = dat2spp%>% 
  mutate(xmag = sqrt(((x1/x01-1))^2 + ((x2/x02-1))^2),
         xmagnorm = xmag/rho,
         xmaglin = sqrt((x1l/x01-1)^2 + (x2l/x02-1)^2),
         xmaglinnorm = xmaglin/rho,
         homo = dtheta<2e-1)


data = dat %>% filter(sim == 2)

#rs 2spp
ggplot(data, aes(x = r1, y = r2))+
  geom_point()+
  geom_point(data = datagen, 
             aes(x = rvec_1, 
                 y = rvec_2),
             color = "gray",
             alpha = 0.5)+
  facet_wrap(~alpha)

#eqs 2spp
ggplot(data, aes(x = x1/(x01), y = x2/(x02),
                       color = xmag))+
  geom_point(size = 1)+
  geom_point(data = datagen, 
             aes(x = xeq_1,
                 y = xeq_2),
             color = "gray", 
             alpha = 0.5)+
  # geom_point(data = data,
  #            aes(x = x1l/(x01), y = x2l/(x02)),
  #            color = "gray",
  #            size = 0.5,
  #            alpha = 0.7)+
  coord_cartesian(xlim = c(0, NA), ylim = c(0,NA))+
  facet_wrap(~alpha,
             #scales = "free"
  )+
  scale_color_viridis_c(option = "plasma")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

##############################################################

#check the decreasing magnitude

datmin = dat2sppgen %>%
  filter(dtheta < 1e-6) %>% 
  group_by(sim, alpha, idrowvec, eqid) %>% 
  mutate(xmag = sqrt(sum(xeq^2))) %>% 
  ungroup()

datminbin = datmin %>% 
  mutate(bin = cut(drho, seq(min(drho), max(drho), (max(drho)-min(drho))/5), right = FALSE)) %>% 
  group_by(alpha, bin) %>% 
  summarize(avmag = median(xmag))

ggplot(datminbin, aes(x = alpha, y = avmag,
                      color = bin))+
  geom_point()

#now with the other dataset

datmin = dat2spp %>% 
  filter(dtheta < 1e-6) %>% 
  group_by(sim , alpha) %>% 
  mutate(xmag = sqrt(x1^2 + x2^2)) %>% 
  ungroup() %>% 
  group_by(alpha, rho) %>% 
  summarise(avmag = median(xmag))

ggplot(datmin, aes(x = alpha, y = avmag)) +
  geom_point()+
  facet_wrap(~rho)

datmags = dat2sppgen %>%
  filter(dtheta < 1e-6) %>% 
  group_by(sim, alpha, idrowvec, eqid) %>% 
  summarise(mag = sqrt(sum(xeq^2))-sqrt(2),
            rho = sqrt(sum(rvec))) %>% 
  ungroup() %>% 
  group_by(alpha, idrowvec, eqid) %>% 
  summarise(avmag = median(mag))

ggplot(data = datmags, 
       aes(x = alpha, y = avmag))+
  geom_point()
  
