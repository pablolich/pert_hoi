library("tidyverse")

dat = read.table("../data/feasibility_allcons_niceviz.csv")
colnames(dat) = c("sim", "n", "alpha", "idrowvec", "eqid", "rvec",
                  "cubeind", "dtheta", "drho", "sppid", "xeq")

###########################################################################

#get largest feasible cube for each alpha

dat_feasible_cubes = dat %>% 
  group_by(alpha, sim, n, cubeind, eqid) %>% 
  filter(all(xeq > 0))

#get maximum value of drho for each feasible group

dat_largest_cube = dat_feasible_cubes %>% 
  group_by(sim, alpha, n, cubeind) %>% 
  mutate(maxcube = max(cubeind)) %>% 
  ungroup() %>% 
  group_by(alpha, n) %>% 
  summarise(cubeindav = mean(maxcube))

#plot 

ggplot(dat_largest_cube, 
       aes(x = alpha, y = cubeindav))+
  geom_point()

############################################################################

#get a small sample of the data

#get statistics on the radii

#get all rows containing at least one feasible equilibrium in the small radius
dat_feasible = dat %>%
  #filter(dtheta < 1e-6) %>% 
  group_by(idrowvec, eqid) %>% 
  filter(all(xeq > 0))

#get average length of equilibria and pertrubations
dat_radii = dat_feasible %>% 
  group_by(idrowvec, eqid) %>% 
  mutate(xmag = sqrt(sum(xeq^2)),
         rmag = sqrt(sum(rvec^2)),
         deltaxmag = sqrt(sum((xeq-1)^2))) %>% 
  ungroup() %>% 
  group_by(idrowvec) %>% 
  slice_min(deltaxmag) %>% 
  mutate(neq = max(eqid)) %>% 
  group_by(alpha, n, cubeind) %>% 
  summarise(xmagav = median(xmag),
            rmagav = median(rmag),
            neqav = median(neq))

ggplot(dat_radii,
       aes(x = alpha, 
           y = xmagav)) +
  geom_point() +
  facet_wrap(~cubeind)

ggplot(dat_radii,
       aes(x = alpha, 
           y = rmagav)) +
  geom_point() 

ggplot(dat_radii, 
       aes(x = alpha, 
           y = neqav))+
  geom_point()