library("tidyverse")

datcube = read.table("../data/cubenumbertest.csv")
colnames(datcube) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                         "cubeind", "dtheta", "drho", "sppid", "xeq")

datcubewide = datcube %>%
  mutate(rowsid = idrowvec + eqid) %>% 
  mutate(bin = cut(drho, seq(min(drho), max(drho), (max(drho)-min(drho))/11), right = FALSE)) %>% 
  group_by(sim, alpha) %>% 
  pivot_wider(names_from = sppid,
              values_from = c(xeq, rvec)) %>% 
  filter(xeq_1 > 0 & xeq_2 > 0) %>% 
  mutate(xmag = sqrt(xeq_1^2 + xeq_2^2))


ggplot(datcubewide, aes(x = rvec_1, y = rvec_2,
                        color = as.factor(cubeind)))+
  geom_point()+
  facet_wrap(~alpha)
