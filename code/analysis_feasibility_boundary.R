library(tidyverse)

circleFun <- function(center = c(0,0), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

datcirc <- circleFun(c(0,0),10,npoints = 100)
#geom_path will do open circles, geom_polygon will do filled circles


data2spp = read.table("../data/feasibility_boundary2sppconst.csv")

colnames(data2spp) = c("sim", "alpha",  "r1", "r2", "limit", "eqid", "area")

dataplot = data2spp %>% 
  mutate(rho = sqrt(r1^2 + r2^2),
         theta = atan2(r1,r2)*180/pi)


ggplot(dataplot,
       aes(theta, rho,
           color = as.factor(eqid)))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray")+
  geom_point(size = 0.7)+
  scale_color_brewer(palette = "Set3")+
  coord_polar()+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")
  
#############################################################

#get areas of each polygon

areasdata = data2spp %>% 
  group_by(alpha, eqid, sim) %>% 
  slice_max(area) %>% 
  group_by(alpha, eqid) %>% 
  filter(!(r1 == 0 & r2 == 0)) %>% 
  summarize(avarea = median(area))

area_small = ggplot(areasdata,
       aes(x = alpha,
           y = avarea,
           color = as.factor(eqid)))+
  geom_point()


###############################################################################


data2sppbig = read.table("../data/feasibility_boundary2sppconstbigdom.csv")

colnames(data2sppbig) = c("sim", "alpha",  "r1", "r2", "limit", "eqid", "area")

dataplotbig = data2sppbig %>% 
  mutate(rho = sqrt(r1^2 + r2^2),
         theta = atan2(r1,r2)*180/pi) %>% 
  filter(sim == 3)


ggplot(dataplotbig,
       aes(theta, rho,
           color = as.factor(eqid)))+
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray")+
  geom_point(size = 0.7)+
  scale_color_brewer(palette = "Set3")+
  coord_polar()+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

#############################################################

#get areas of each polygon

areasdatabig = data2sppbig %>% 
  group_by(alpha, eqid, sim) %>% 
  slice_max(area) %>% 
  group_by(alpha, eqid) %>% 
  filter(!(r1 == 0 & r2 == 0)) %>% 
  summarize(avarea = median(area))

area_big = ggplot(areasdatabig,
       aes(x = alpha,
           y = avarea,
           color = as.factor(eqid)))+
  geom_point()


###############################################################################


data2sppeqplanted = read.table("../data/feasibility_boundary2sppeqplanted.csv")

colnames(data2sppeqplanted) = c("sim", "alpha",  "r1", "r2", "limit", "eqid", "area")

dataploteqplanted = data2sppeqplanted %>% 
  mutate(rho = sqrt(r1^2 + r2^2),
         theta = atan2(r1,r2)*180/pi) %>% 
  filter(sim == 3)

ggplot(dataploteqplanted,
       aes(theta, rho,
           color = as.factor(eqid)))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray")+
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray")+
  geom_point(size = 0.7)+
  scale_color_brewer(palette = "Set3")+
  coord_polar()+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")

areasdataeqplanted = data2sppeqplanted %>% 
  group_by(alpha, eqid, sim) %>% 
  slice_max(area) %>% 
  group_by(alpha, eqid) %>% 
  filter(!(r1 == 0 & r2 == 0)) %>% 
  summarize(avarea = median(area))

area_eqplanted = ggplot(areasdataeqplanted,
                  aes(x = alpha,
                      y = avarea,
                      color = as.factor(eqid)))+
  geom_point()



area_big
area_small
area_eqplanted
