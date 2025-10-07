library(dplyr)
library(tidyr)

# Extract A and B for each seed
dat_param <- read.csv("../../data/all_parameter_sets_n6.csv") # adapt to your path

# Assume order == 1 is A, order == 2 is B
A_stats <- dat_param %>%
  filter(order == 1) %>%
  group_by(seed_i) %>%
  summarise(
    mean_A = mean(value),
    sd_A   = sd(value),
    min_A  = min(value),
    max_A  = max(value),
    frac_pos_A = mean(value > 0),
    frac_neg_A = mean(value < 0)
  )

B_stats <- dat_param %>%
  filter(order == 2) %>%
  group_by(seed_i) %>%
  summarise(
    mean_B = mean(value),
    sd_B   = sd(value),
    min_B  = min(value),
    max_B  = max(value),
    frac_pos_B = mean(value > 0),
    frac_neg_B = mean(value < 0)
  )

# Read transitions/robustness summary (already has p_robust, p_reactive, ecosystem_class)
transitions_summary <- read.csv("transitions_summary.csv")

# Join all together by seed
main_df <- transitions_summary %>%
  left_join(A_stats, by = c("seed" = "seed_i")) %>%
  left_join(B_stats, by = c("seed" = "seed_i")) 

# If you want to include the entire transition matrix per seed, you can just keep the transition columns
# (e.g., transition_boundary_to_negative, etc.), or add any transition-specific summary stat you like

# View head/main columns
head(main_df)

# List of variables to color by
color_vars <- c(
  "mean_A", "sd_A", "min_A", "max_A", "frac_pos_A",
  "frac_neg_A", "mean_B", "sd_B", "min_B", "max_B",
  "frac_pos_B", "frac_neg_B"
)

# Generate plots
plots <- lapply(color_vars, function(var) {
  ggplot(main_df, aes(x = p_robust, y = p_reactive, color = !!sym(var))) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_viridis_c(option = "plasma") +
    labs(title = paste("Colored by", var), color = var) +
    theme_bw(base_size = 12)
})

# Arrange plots in a grid and save
pdf("p_robust_vs_p_reactive_colored.pdf", width = 12, height = 14)
wrap_plots(plots, ncol = 3)
dev.off()


library(ggplot2)
library(patchwork)

# List of summary statistics for x-axis
x_vars <- c(
  "mean_A", "sd_A", "min_A", "max_A", "frac_pos_A",
  "frac_neg_A", "mean_B", "sd_B", "min_B", "max_B",
  "frac_pos_B", "frac_neg_B"
)

# Generate plots: x-axis = summary variable, y-axis = p_robust, color = p_reactive
plots <- lapply(x_vars, function(var) {
  ggplot(main_df, aes_string(x = var, y = "p_robust", color = "p_reactive")) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = paste("p_robust vs", var, "| Color: p_reactive"),
      x = var, y = "p_robust", color = "p_reactive"
    ) +
    theme_bw(base_size = 12)
})

# Save all plots to a PDF in a grid layout
pdf("robustness_vs_summary_stats.pdf", width = 14, height = 16)
wrap_plots(plots, ncol = 3)
dev.off()


library(ggplot2)
library(patchwork)

# Vector of statistic names
stats <- c("mean", "sd", "min", "max", "frac_pos", "frac_neg")

# Generate plots colored by p_robust
plots_robust <- lapply(stats, function(stat) {
  x_var <- paste0(stat, "_A")
  y_var <- paste0(stat, "_B")
  ggplot(main_df, aes_string(x = x_var, y = y_var, color = "p_robust")) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = paste(x_var, "vs", y_var, "| Color: Robustness"),
      x = x_var, y = y_var, color = "Robustness"
    ) +
    theme_bw(base_size = 12)
})

# Save robustness-colored plots
pdf("A_vs_B_stats_colored_by_robustness.pdf", width = 14, height = 10)
wrap_plots(plots_robust, ncol = 3)
dev.off()

# Generate plots colored by p_reactive
plots_reactive <- lapply(stats, function(stat) {
  x_var <- paste0(stat, "_A")
  y_var <- paste0(stat, "_B")
  ggplot(main_df, aes_string(x = x_var, y = y_var, color = "p_reactive")) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = paste(x_var, "vs", y_var, "| Color: Reactivity"),
      x = x_var, y = y_var, color = "Reactivity"
    ) +
    theme_bw(base_size = 12)
})

# Save reactivity-colored plots
pdf("A_vs_B_stats_colored_by_reactivity.pdf", width = 14, height = 10)
wrap_plots(plots_reactive, ncol = 3)
dev.off()


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tidyverse)

# List transition pairs (in your data, 'From' and 'To')
transition_pairs <- main_df %>%
  select(From, To) %>%
  distinct() %>%
  arrange(From, To) %>%
  unite("transition", From, To, sep = "_") %>%
  pull(transition)

# List of B and A stats to use
B_stats <- c("frac_pos_B", "frac_neg_B")
A_stats <- c("frac_pos_A", "frac_neg_A")

# Add a transition column for easier plotting
main_df <- main_df %>%
  mutate(transition = paste(From, To, sep = "_"))

plot_transition_vs_stat <- function(stat, stat_label = NULL) {
  plots <- lapply(transition_pairs, function(trans) {
    ggplot(main_df %>% filter(transition == trans),
           aes_string(x = paste0("factor(", stat, ")"), y = "Proportion")) +
      geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
      geom_boxplot(fill = "lightblue", alpha = 0.5, outlier.shape = NA) +
      labs(
        title = paste("Transition:", trans),
        x = ifelse(is.null(stat_label), stat, stat_label),
        y = "Proportion"
      ) +
      theme_bw(base_size = 12)
  })
  return(plots)
}

for (stat in B_stats) {
  plots <- plot_transition_vs_stat(stat)
  g <- wrap_plots(plots, ncol = 3)
  ggsave(
    filename = paste0("transition_vs_", stat, "_B.png"),
    plot = g,
    width = 14,
    height = 10,
    dpi = 300
  )
}


for (stat in A_stats) {
  plots <- plot_transition_vs_stat(stat)
  g <- wrap_plots(plots, ncol = 3)
  ggsave(
    filename = paste0("transition_vs_", stat, "_A.png"),
    plot = g,
    width = 14,
    height = 10,
    dpi = 300
  )
}
