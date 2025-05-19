library(tidyverse)

# Read and name columns
data <- read.csv("../../data/boundary_distances/n_3.csv", header = FALSE, sep = " ")
dat_param = read.csv("../../data/all_parameter_sets_n3.csv")

data <- read.csv("../../data/boundary_distances/n_4.csv", header = FALSE, sep = " ")
dat_param = read.csv("../../data/all_parameter_sets_n4.csv")

# data <- read.csv("../../data/boundary_distances/n_5.csv", header = FALSE, sep = " ")
# dat_param = read.csv("../../data/all_parameter_sets_n5.csv")

plot_all_hists = function(data){
  # Group and split
  data_plot <- data %>%
    group_by(seed_i, n, pert_size, d) %>%
    group_split()
  
  # Output PDF
  pdf("all_histograms_by_seed.pdf", width = 8, height = 6)
  
  # Loop through groups and plot
  for (group in data_plot) {
    seed <- unique(group$seed_i)
    n <- unique(group$n)
    p <- unique(group$pert_size)
    
    p_hist <- ggplot(group) +
      geom_histogram(aes(x = deltap, fill = as.factor(alpha_i)), 
                     bins = 20, 
                     color = "black", 
                     alpha = 0.5,
                     position = "identity") +
      labs(
        title = paste("Seed =", seed, "| n =", n, "| Perturbation =", p),
        x = "deltap", y = "Count", fill = "alpha"
      ) +
      theme_minimal()
    
    print(p_hist)  # adds to current PDF page
  }
  
  dev.off()
}
library(tidyverse)

u_thresh = 0.5
l_thresh = 0.5

process_n_data <- function(n_val) {
  # Load data
  data <- read.csv(paste0("../../data/boundary_distances/n_", n_val, ".csv"), header = FALSE, sep = " ")
  dat_param <- read.csv(paste0("../../data/all_parameter_sets_n", n_val, ".csv"))
  
  # Name columns
  colnames(data) <- c("seed_i", "n", "d", "pert_size", "pert_i", "alpha_i", "deltap", "deltalin", "deltax", "deltaxlin", "flag")
  
  # Pivot wider and classify
  data_wide <- data %>%
    mutate(alpha_label = paste0("alpha_", alpha_i)) %>%
    select(seed_i, n, d, pert_size, pert_i, alpha_label, deltap) %>%
    pivot_wider(names_from = alpha_label, values_from = deltap) %>%
    mutate(diff_delta = alpha_0.5 - alpha_0.1) %>%
    filter(diff_delta != 0) %>%
    mutate(robust = diff_delta > 0) %>%
    ungroup() %>%
    group_by(seed_i) %>%
    summarise(av_robust = mean(robust)) %>%
    mutate(threshold_robust = ifelse(av_robust >= u_thresh, 1, ifelse(av_robust < l_thresh, -1, 0))) %>%
    filter(threshold_robust != 0)
  
  # Merge and annotate n
  dat_param %>%
    left_join(data_wide, by = "seed_i") %>%
    drop_na(threshold_robust) %>%
    mutate(n = n_val)
}

# Process for all desired n
all_data <- bind_rows(
  process_n_data(3),
  process_n_data(4),
  process_n_data(5)
)

# Count -1 and 1 per n
robust_counts <- all_data %>%
  group_by(n, threshold_robust, seed_i) %>%
  slice_head() %>% 
  group_by(n, threshold_robust) %>% 
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("class", threshold_robust, "count=", count),
         y = ifelse(threshold_robust == 1, 2, 1.2))

# Base plot
p <- ggplot(all_data) +
  geom_histogram(aes(x = value, fill = as.factor(threshold_robust), after_stat(density)),
                 position = "identity", alpha = 0.5) +
  geom_density(aes(x = value, color = as.factor(threshold_robust))) +
  facet_grid(n ~ order) +
  xlim(-1, 1) +
  labs(x = "Coefficient Value", fill = "Threshold Robust", color = "Threshold Robust") +
  theme_minimal()

# Annotate counts on plot — one label per n + class, placed at a fixed x, y position
p + geom_text(data = robust_counts, 
              aes(x = -0.8, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 3.5)

table(count_each$threshold_robust)
