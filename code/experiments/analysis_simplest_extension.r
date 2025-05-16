library(tidyverse)
library(gridExtra)

dat2 = read.table("../../data/results_simplest_extension_n_2.csv", sep = "")
dat3 = read.table("../../data/results_simplest_extension_n_3.csv", sep = "")
dat4 = read.table("../../data/results_simplest_extension_n_4.csv", sep = "")
dat5 = read.table("../../data/results_simplest_extension_n_5.csv", sep = "")
dat6 = read.table("../../data/results_simplest_extension_n_6.csv", sep = "")
dat7 = read.table("../../data/results_simplest_extension_n_7.csv", sep = "")

dat = rbind(dat2, dat3, dat4, dat5, dat6, dat7, dat8)
colnames(dat) = c("seed", "n", "alpha", "pert_dir", "delta")


datplot = dat %>% filter(delta < 4) %>% 
  mutate(alpha = factor(alpha, levels = c(0.1, 0.9)))

###########################################################################

#Analysis histogram by histogram, n by n

datplot %>%
  group_by(n) %>%
  do({
    # For each value of n, create a separate PDF file
    n_value <- unique(.$n)  # Get the current 'n' value
    
    # Open a new PDF file for the current 'n'
    pdf_file <- paste0("histograms_n_", n_value, ".pdf")
    pdf(pdf_file, height = 3, width = 8)  # Open the PDF device with width to fit both plots
    
    # Iterate over each unique 'seed' for the current 'n'
    lapply(unique(.$seed), function(seed_value) {
      # Filter data for the current seed
      seed_data <- filter(., seed == seed_value)
      
      # Create the histogram plot
      histogram_plot <- ggplot(seed_data, aes(x = delta, fill = factor(alpha))) +
        geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
        labs(title = paste("Histogram of delta by alpha for seed =", seed_value), 
             x = "Delta", 
             y = "Frequency") +
        theme_minimal()
      
      if (seed == 838){
        browser()
      }
      seed_data_wide <- seed_data %>% 
        mutate(row = row_number()) %>% 
        pivot_wider(names_from = alpha, values_from = delta, names_prefix = "alpha_") %>%
        group_by(seed, n, pert_dir) %>%
        mutate(alpha_0.9 = c(alpha_0.9[-1], alpha_0.9[1])) %>%
        na.omit() %>% 
        mutate(diff = alpha_0.9 - alpha_0.1,
               sign = ifelse(diff >= 0, "Positive", "Negative"))
      
      histogram_dist_var = ggplot(seed_data_wide, aes(x = alpha_0.9 - alpha_0.1,
                                                      fill = sign))+
        geom_histogram(bins = 30,
                       position = "identity",
                       alpha = 0.5)
      
      # Arrange both plots side by side
      grid.arrange(histogram_plot, histogram_dist_var, ncol = 2)
    })
    
    # Close the PDF file after adding all the plots
    dev.off()
    
    # Return the grouped data (optional)
    .
  })

############################################################################
#ANALYSIS OF DISTANCES

#analysis of compounded plots (all simulations by n)

h =
  ggplot(datplot, aes(x = delta,
                fill = alpha)) +
  geom_histogram(position="fill", 
                 alpha = 0.5,
                 bins = 50)+
  theme(aspect.ratio = 1)+
  facet_wrap(~n,
             nrow = 2,
             scales = "free")

print(h)

plot_data = ggplot_build(h)$data[[1]] %>% 
  select(c(x, ymin, PANEL)) %>% 
  filter(ymin>0) %>% 
  rename(n = PANEL)

ggplot(plot_data %>% filter(x < 2),
       aes(x = x,
           y = ymin)) +
  geom_point(aes(color = n))+
  geom_line(aes(group = n, col = n))+
  geom_hline(yintercept = 0.5, size = 1.2, linetype = "dashed")+
  theme(aspect.ratio = 0.5)+
  labs(x = "Distance",
       y = "")

#histograms of distance change of pairs of systems with two values of alpha

datplot_wide <- datplot %>% 
  mutate(row = row_number()) %>% 
  pivot_wider(names_from = alpha, values_from = delta, names_prefix = "alpha_") %>%
  group_by(seed, n, pert_dir) %>%
  mutate(alpha_0.9 = c(alpha_0.9[-1], alpha_0.9[1])) %>%
  na.omit() %>% 
  mutate(diff = alpha_0.9 - alpha_0.1,
         sign = ifelse(diff >= 0, "Positive", "Negative"))

ggplot(datplot_wide,
       aes(x = alpha_0.9 - alpha_0.1,
           fill = as.factor(sign)))+
  geom_histogram(position="identity",
                 alpha = 0.5)+
  facet_wrap(~n,
             nrow = 2)

############################################################################
#ANALYSIS OF MOMENTS OF DISTRIBUTION OF DISTANCES
#scatterplots of the moments: mean and variance of all distributions compounded

moments_data = datplot %>% 
  group_by(seed, n, alpha) %>% 
  summarise(mean_dist = median(delta),
            var_dist = sd(delta))

ggplot(moments_data)+
  geom_bin_2d(aes(x = mean_dist,
                  y = var_dist,
                  fill = as.factor(alpha)),
              alpha = 0.5)+
  theme(aspect.ratio = 1)+
  facet_wrap(~n)

#scatter plot of the change in moments between pairs of models with weak and
#strong HOIs

moments_diff <- moments_data %>%
  # Spread data so that alpha 0.1 and 0.9 become separate columns
  pivot_wider(names_from = alpha, values_from = c(mean_dist, var_dist)) %>%
  # Calculate the differences for each group (seed, n)
  mutate(
    mean_dist_diff = `mean_dist_0.9` - `mean_dist_0.1`,
    var_dist_diff = `var_dist_0.9` - `var_dist_0.1`
  ) %>%
  # Select relevant columns for the output
  select(seed, n, mean_dist_diff, var_dist_diff)

  ggplot(moments_diff)+
  geom_bin_2d(aes(x = mean_dist_diff,
              y = var_dist_diff,
              bins = 50))+
  geom_vline(xintercept = 0, color = "grey")+
  geom_hline(yintercept = 0, color = "grey")+
  theme(aspect.ratio = 1)+
  facet_wrap(~n,
             nrow = 3)
  
###############################################################################
#ANALYSIS OF ABUNDANCES AT FEASIBILITY BOUNDARY
  
#distributions of abundances for all simulations when HOIs are weak and strong.

#Change in abundance distribution as alpha increases.

###############################################################################
#ANALYSIS OF ROBUSTNESS CHANGE 
  
#Here I look at the robustness change in two ways: 
  #1. Increase in mean of distribution of distances (susceptible to outliers)
  #2. Increase in median of distribution of distances 
  #(not susceptible to outliers)
  
##############################################################################
#ANALYSIS OF THE PARAMETERS
  #Analyze parameters of model conditioned on the robustness increase
  #as per either of the cases.