library(tidyverse)
library(gridExtra)
library(RColorBrewer)

library(extrafont) 
# link www.fontsquirrel.com/fonts/latin-modern-roman

# execute once to add fonts:
#font_import(pattern = "lmroman*") 


#library(tune)

dat = read.table("../data/feasibility_allcons_niceviz.csv")
colnames(dat) = c("sim", "alpha", 
                  "dtheta", "theta", "rho", 
                  "r01", "r02",
                  "r1", "r2",
                  "x1", "x2",
                  "x01", "x02",
                  "x1l", "x2l",
                  "area")

dat_small = dat %>% filter(sim == 43,
                           rho > 0 & rho < 2e-1,
                             alpha == .4| alpha == .7 | alpha == .8)

ggplot(dat_small)+
  geom_point(aes(x = r1, y = r2))+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")
ggsave("../data/linear_map_example_r.pdf", width = 4, height = 4)

ggplot(dat_small)+
  geom_point(aes(x = x1, y = x2),
             color = "red")+
  geom_point(data = dat_small,
             aes(x = x1l, y = x2l),
             color = "black",
             size = 0.3)+
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed")+
  facet_wrap(~alpha)+
  scale_color_viridis_c(option = "plasma")+
  theme(aspect.ratio = 1,
        #text = element_text(size=10, family="LM Roman 10"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        legend.position = "none")
  ggsave("../data/linear_map_example.pdf", width = 6, height = 2.5)
  
