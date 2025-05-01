library(tidyverse)
library(gridExtra)

# dat2 = read.table("../data/results_simplest_extension_n_2.csv", sep = "")
dat = read.table("../data/results_simplest_extension_n_3.csv", sep = "")
# dat4 = read.table("../data/results_simplest_extension_n_4.csv", sep = "")
# dat5 = read.table("../data/results_simplest_extension_n_5.csv", sep = "")
# dat6 = read.table("../data/results_simplest_extension_n_6.csv", sep = "")
# dat7 = read.table("../data/results_simplest_extension_n_7.csv", sep = "")
# dat8 = read.table("../../data/results_simplest_extension_n_8.csv", sep = "")

#dat = rbind(dat2, dat3, dat4, dat5, dat6, dat7, dat8)
colnames(dat) = c("seed", "n", "alpha", "pert_dir", "delta")

#colnames(dat) = c("seed", "n", "d", "pert_size", "pert_dir", "alpha", "delta", "delta_lin", "delta_x", "delta_x_lin")


datplot = dat %>% filter(#delta < 4,
                         seed >= 1) %>% 
  mutate(alpha = factor(alpha, levels = c(0.1, 0.9)))

#load parameters

#datpar = read.table("../data/parameters_n_2", sep = "")
###########################################################################

#Analysis histogram by histogram, n by n
pdf("histograms_n_3.pdf", height = 3, width = 8)

# Loop over each seed for n = 3
unique(datplot$seed) %>%
  walk(function(seed_value) {
    seed_data <- filter(datplot, seed == seed_value)
    
    # Compute mean and median per alpha
    summary_lines <- seed_data %>%
      group_by(alpha) %>%
      summarise(
        mean = mean(delta, na.rm = TRUE),
        median = median(delta, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Plot histogram with mean/median lines
    print(seed_value)
    histogram_plot <- ggplot(seed_data, aes(x = delta, fill = factor(alpha))) +
      geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
      geom_vline(data = summary_lines, aes(xintercept = median, color = factor(alpha)),
                 linetype = "solid", size = 1, show.legend = FALSE) +
      geom_vline(data = summary_lines, aes(xintercept = mean, color = factor(alpha)),
                 linetype = "dashed", size = 1, show.legend = FALSE) +
      scale_color_manual(values = c("0.1" = "#F8766D", "0.9" = "#00BFC4")) +
      scale_fill_manual(values = c("0.1" = "#F8766D", "0.9" = "#00BFC4")) +
      labs(
        title = paste("Histogram of delta by alpha for seed =", seed_value),
        x = "Delta",
        y = "Frequency"
      ) +
      theme_minimal()
    
    print(histogram_plot)
  })

dev.off()

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