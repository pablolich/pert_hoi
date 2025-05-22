library(tidyverse)

# Read and name columns
# data <- read.csv("../../data/boundary_distances/n_3.csv", header = FALSE, sep = " ")
# dat_param = read.csv("../../data/all_parameter_sets_n3.csv")
# 
# data <- read.csv("../../data/boundary_distances/n_4.csv", header = FALSE, sep = " ")
# dat_param = read.csv("../../data/all_parameter_sets_n4.csv")
# 
data <- read.csv("../../data/boundary_distances/n_5.csv", header = FALSE, sep = " ")
dat_param = read.csv("../../data/all_parameter_sets_n5.csv")

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
    
    group$flag <- factor(group$flag, levels = c("boundary", "complex", "negative", "nonconverged"))
    
    # Define color mapping
    flag_colors <- c(
      "negative"     = "#377EB8",  # blue
      "complex"      = "#984EA3",  # purple
      "boundary"     = "#4DAF4A",  # green
      "nonconverged" = "#E41A1C"   # red
    )
    
    # Plot
    p_hist <- ggplot(group) +
      geom_histogram(aes(x = deltap, fill = flag),
                     bins = 20,
                     color = "black",
                     alpha = 0.7) +
      facet_wrap(~ alpha_i, ncol = 1, labeller = label_both) +
      scale_fill_manual(values = flag_colors, drop = FALSE) +
      labs(
        title = paste("Seed =", seed, "| n =", n, "| Perturbation =", p),
        x = "Δp", y = "Count", fill = "Flag"
      ) +
      theme_minimal()
    
    print(p_hist)
    
  }
  
  dev.off()
}
library(tidyverse)

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

get_boundary_composition_change <- function(n_val, l_thresh = 0.5, u_thresh = 0.5) {

  # Load data
  data <- read.csv(paste0("../../data/boundary_distances/n_", n_val, ".csv"), header = F, sep=" ")
  # Name columns
  colnames(data) <- c("seed_i", "n", "d", "pert_size", "pert_i", "alpha_i", 
                      "deltap", "deltalin", "deltax", "deltaxlin", "flag")
  
  # Prepare wide format to compute difference in delta_p
  data_classified <- data %>%
    mutate(alpha_label = paste0("alpha_", alpha_i)) %>%
    select(seed_i, n, d, pert_size, pert_i, alpha_label, deltap) %>%
    pivot_wider(names_from = alpha_label, values_from = deltap) %>%
    mutate(diff_delta = alpha_0.5 - alpha_0.1) %>%
    filter(!is.na(diff_delta) & diff_delta != 0) %>%
    mutate(robust = diff_delta > 0)
  
  # Compute average robustness and threshold decision per system
  robust_summary <- data_classified %>%
    group_by(seed_i) %>%
    summarise(av_robust = mean(robust), .groups = "drop") %>%
    mutate(threshold_robust = ifelse(av_robust >= u_thresh, 1,
                                     ifelse(av_robust < l_thresh, -1, 0))) %>%
    filter(threshold_robust != 0)
  
  # Compute flag proportions per seed and alpha
  flag_props <- data %>%
    group_by(seed_i, alpha_i, flag) %>%
    summarise(count = n(), .groups = "drop_last") %>%
    mutate(total = sum(count)) %>%
    ungroup() %>%
    mutate(prop = count / total) %>%
    select(seed_i, alpha_i, flag, prop) %>%
    pivot_wider(names_from = flag, values_from = prop, values_fill = 0)
  
  # Join everything
  final_output <- robust_summary %>%
    left_join(flag_props, by = c("seed_i")) %>% 
    mutate(n = n_val)
  
  return(final_output)
}


# Process for all desired n
all_data <- bind_rows(
  get_robustness_parameters(3),
  get_robustness_parameters(4),
  get_robustness_parameters(5)
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
  xlim(-5, 5) +
  labs(x = "Coefficient Value", fill = "Threshold Robust", color = "Threshold Robust") +
  theme_minimal()

# Annotate counts on plot — one label per n + class, placed at a fixed x, y position
p + geom_text(data = robust_counts, 
              aes(x = -1, y = y, label = label),
              inherit.aes = FALSE,
              hjust = 0, size = 3.5)

# Process for all desired n
all_data_comp_change <- bind_rows(
  get_boundary_composition_change(3),
  get_boundary_composition_change(4),
  get_boundary_composition_change(5)
)

library(tidyverse)

# Your fixed color palette
flag_colors <- c(
  "negative"     = "#377EB8",  # blue
  "complex"      = "#984EA3",  # purple
  "boundary"     = "#4DAF4A",  # green
  "nonconverged" = "#E41A1C"   # red
)

# Prepare long-format data
df_long <- df %>%
  pivot_longer(
    cols = c(boundary, complex, negative, nonconverged),
    names_to = "edge_type",
    values_to = "proportion"
  ) %>%
  mutate(
    robust_label = ifelse(threshold_robust == 1, "More Robust", "Less Robust"),
    alpha_label = paste0("alpha = ", alpha_i)
  )

# Plot stacked histograms (percentage of systems)
ggplot(df_long, aes(x = proportion, fill = edge_type)) +
  geom_histogram(position = "fill", bins = 10, color = "black", alpha = 0.85) +
  facet_grid(rows = vars(robust_label), cols = vars(alpha_label)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = flag_colors, drop = FALSE) +
  labs(
    x = "Prevalence",
    y = "Boundary composition",
    fill = "Boundary Type",
    title = "Boundary Composition Distributions by Robustness and Alpha"
  ) +
  theme_minimal()
