library(tidyverse)
setwd("/home/srilenak/Desktop/pert_hoi/code/experiments")  # Set working directory to the parent directory

u_thresh = 0.5
l_thresh = 0.5

get_robustness_parameters <- function(n_val) {
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
  get_robustness_parameters(3),
  get_robustness_parameters(4),
  get_robustness_parameters(5), 
  get_robustness_parameters(6)
)

# Count -1 and 1 per n
robust_counts <- all_data %>%
  group_by(n, threshold_robust, seed_i) %>%
  slice_head() %>% 
  group_by(n, threshold_robust) %>% 
  summarise(count = n(), .groups = "drop") %>%
  mutate(label = paste0("class ", threshold_robust, " count=", count),
         y = ifelse(threshold_robust == 1, 2.5, 2.1))

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
              aes(x = -1, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 3.5)

###############################################################################

n_val <- 6
ecosystem_class_df <- read_csv("ecosystem_class_summary.csv")

dat_param <- read.csv(paste0("../../data/all_parameter_sets_n", n_val, ".csv"))
# Rename and merge
dat_param <- dat_param %>%
  rename(seed = seed_i)

merged_df <- inner_join(dat_param, ecosystem_class_df, by = "seed")

# Plot density plots for each 'order'
ggplot(merged_df, aes(x = value, color = ecosystem_class, fill = ecosystem_class)) +
  geom_density(alpha = 0.3) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~ order, scales = "free", nrow = 3) +
  labs(title = paste("Distribution of Parameter Values by Ecosystem Class (n =", n_val, ")"),
       x = "Parameter Value", y = "Density", color = "Class", fill = "Class") +
  theme_minimal(base_size = 12)
  #xlim(-0.5, 0.5)

