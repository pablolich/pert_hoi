library(tidyverse)

dat2 = read.table("../../data/results_simplest_extension_n_2.csv", sep = "")
dat3 = read.table("../../data/results_simplest_extension_n_3.csv", sep = "")
dat4 = read.table("../../data/results_simplest_extension_n_4.csv", sep = "")
dat5 = read.table("../../data/results_simplest_extension_n_5.csv", sep = "")


dat = rbind(dat2, dat3, dat4, dat5)
colnames(dat) = c("seed", "n", "alpha", "pert_dir", "delta")


datplot = dat %>% filter(delta <= 2) %>% 
  mutate(alpha = factor(alpha, levels = c(0.1, 0.9)))

dathead = head(dat, n = 1e4)
datwide = 
  dat %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = alpha,
              values_from = delta,
              values_fill = 0,) %>% 
  rename(low = `0.1`,
         high = `0.9`) %>% 
  group_by(seed, n, pert_dir) %>% 
    summarise(low = max(low),
              high = max(high),
              diff = high - low)

ggplot(datwide,
       aes(x = high - low))+
  geom_histogram()+
  facet_wrap(~n)

h = ggplot(datplot, aes(x = delta,
                fill = alpha,
                after_stat(density))) +
  geom_histogram(position="fill", 
                 alpha = 0.5,
                 bins = 20)+
  theme(aspect.ratio = 1)+
  facet_wrap(~n,
             nrow = 2,
             scales = "free")

plot_data = ggplot_build(h)$data[[1]] %>% 
  select(c(x, ymin, PANEL)) %>% 
  filter(ymin>0) %>% 
  rename(pan = PANEL)

ggplot(plot_data,
       aes(x = x,
           y = ymin)) +
  geom_point(aes(color = pan))+
  geom_line(aes(group = pan, col = pan))+
  theme(aspect.ratio = 0.5)+
  labs(x = "Distance",
       y = "")
