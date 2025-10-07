library(tidyverse)

u_thresh = 0.5
l_thresh = 0.5


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
all_data_comp_change <- bind_rows(
  get_boundary_composition_change(3),
  get_boundary_composition_change(4),
  get_boundary_composition_change(5),
  get_boundary_composition_change(6)
)


# Your fixed color palette
flag_colors <- c(
  "negative"     = "#377EB8",  # blue
  "complex"      = "#984EA3",  # purple
  "boundary"     = "#4DAF4A",  # green
  "nonconverged" = "#E41A1C"   # red
)

# Prepare long-format data
df_long <- all_data_comp_change %>%
  filter(n==6) %>% 
  pivot_longer(
    cols = c(boundary, complex, negative, nonconverged),
    names_to = "edge_type",
    values_to = "proportion"
  ) %>%
  filter(edge_type != "nonconverged", 
         edge_type != "boundary") %>% 
  mutate(
    robust_label = ifelse(threshold_robust == 1, "More Robust", "Less Robust"),
    alpha_label = paste0("alpha = ", alpha_i)
  )

# Plot stacked histograms (percentage of systems)
ggplot(df_long, aes(x = proportion, fill = edge_type)) +
  geom_histogram(#position = "fill", 
                 bins = 20, color = "black", alpha = 0.85) +
  facet_grid(rows = vars(robust_label), cols = vars(alpha_label)) +
  scale_fill_manual(values = flag_colors, drop = FALSE) +
  labs(
    x = "Prevalence",
    y = "Boundary composition",
    fill = "Boundary Type",
    title = "Boundary Composition Distributions by Robustness and Alpha"
  ) +
  theme_minimal()
