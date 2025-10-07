library(tidyverse)
library(patchwork)

setwd("~/Desktop/pert_hoi/code/experiments")  # Set your working directory
# Load
transitions_summary_df <- read_csv("transitions_summary.csv") %>%
  select(-n, -pert_size, -d)

avg_transition_df <- transitions_summary_df %>%
  filter(From != "nonconverged", To != "nonconverged") %>%
  group_by(ecosystem_class, From, To) %>%
  summarise(Proportion = mean(Proportion), .groups = "drop")
plot_avg_transition_matrix <- function(class_df, class_name) {
  ggplot(class_df, aes(x = To, y = From, fill = Proportion)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = class_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      aspect.ratio = 1
    )
}
library(patchwork)

plot_list <- avg_transition_df %>%
  split(.$ecosystem_class) %>%
  imap(~ plot_avg_transition_matrix(.x, .y))

g <- wrap_plots(plot_list, ncol = 2)

ggsave("avg_transition_matrix_2x2.png", g, width = 10, height = 10, dpi = 300)


library(tidyverse)

# Define color palette (same as before)
ecosystem_colors <- c(
  "non_robust_non_reactive" = "#E41A1C",
  "non_robust_reactive"     = "#377EB8",
  "robust_non_reactive"     = "#4DAF4A",
  "robust_reactive"         = "#984EA3"
)

# Add new columns for grouping
class_counts_grouped <- transitions_summary_df %>%
  distinct(seed, ecosystem_class) %>%  # one row per ecosystem
  mutate(
    robustness = if_else(str_starts(ecosystem_class, "robust"), "robust", "non_robust"),
    reactivity = if_else(str_ends(ecosystem_class, "reactive"), "reactive", "non_reactive")
  ) %>%
  count(robustness, reactivity)

# Plot stacked bar chart
ggplot(class_counts_grouped, aes(x = robustness, y = n, fill = reactivity)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "non_reactive" = ecosystem_colors["robust_non_reactive"],
    "reactive"     = ecosystem_colors["robust_reactive"]
  )) +
  labs(
    x = "Robustness Class",
    y = "Number of Ecosystems",
    fill = "Reactivity"
  ) +
  theme_minimal(base_size = 12)

# Save the plot
ggsave("ecosystem_class_counts.png", width = 8, height = 6, dpi = 300)