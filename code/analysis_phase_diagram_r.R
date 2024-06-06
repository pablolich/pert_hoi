library(tidyverse)

dat = read.table("../data/boundaryportrait.csv")
colnames(dat) = c("sim", "alpha", 
                  "theta", "rho", 
                  "r1", "r2", 
                  "limit", 
                  "x1", "x2", 
                  "eqid",
                  "area")

dat = dat %>% 
  filter(sim == 2)

ggplot(dat, aes(x = x1, y = x2,
                color = as.factor(eqid)))+
  geom_point(size = 0.6)+
  facet_wrap(~alpha,
             scales = "free"
             )+
  scale_color_brewer(palette = "Set2")+
  theme(legend.position = "none",
        aspect.ratio = 1)


ggplot(dat,
       aes(theta, rho,
           color = as.factor(eqid)))+
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray")+
  geom_point(size = 0.6)+
  scale_color_brewer(palette = "Set2")+
  coord_polar()+
  facet_wrap(~alpha)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none")


dat %>% 
  group_by(alpha, eqid, sim) %>% 
  slice_max(area) %>% 
  group_by(alpha, eqid) %>% 
  summarize(avarea = median(area))
