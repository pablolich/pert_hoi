setwd("/home/srilenak/Desktop/pert_hoi/code/experiments")

library(tidyverse)
library(patchwork)

# ---- Robustness/reactivity classifier ----
evaluate_response_type <- function(group, type = c("robustness", "reactivity")) {
  type <- match.arg(type)
  varname <- if (type == "robustness") "deltap" else "deltax"

  df_wide <- group %>%
    select(pert_i, alpha_i, !!sym(varname)) %>%
    pivot_wider(names_from = alpha_i, values_from = !!sym(varname), names_prefix = "alpha_")

  df_wide <- df_wide %>%
    mutate(diff = .data[["alpha_0.5"]] - .data[["alpha_0.1"]],
           increased = diff > 0)

  prop_increased <- mean(df_wide$increased, na.rm = TRUE)

  classification <- if (prop_increased > 0.5) {
    if (type == "robustness") "Robust" else "Reactive"
  } else {
    if (type == "robustness") "Non-Robust" else "Non-Reactive"
  }

  return(list(proportion = prop_increased, classification = classification))
}

make_hist_grid <- function(df_wide, var_name, title_prefix) {
  plots <- list()
  for (f1 in levels(df_wide$flag_0.1)) {
    for (f2 in levels(df_wide$flag_0.5)) {
      df_sub <- df_wide %>%
        filter(flag_0.1 == f1, flag_0.5 == f2) %>%
        mutate(sign = ifelse(.data[[var_name]] < 0, "negative", "positive"))

      p <- ggplot(df_sub, aes_string(x = var_name, fill = "sign")) +
        geom_histogram(bins = 30, color = "black") +
        scale_fill_manual(values = c("negative" = "#E41A1C", "positive" = "#4DAF4A")) +
        #theme_minimal(base_size = 7) +
        labs(
          x = NULL, y = NULL,
          title = paste(f1, "\u2192", f2)  # e.g., "complex → boundary"
        ) +
        theme(
          plot.margin = margin(1, 1, 1, 1, "pt"),
          axis.text = element_text(size = 6),
          axis.title = element_blank(),
          axis.ticks = element_line(size = 0.2),
          # panel.grid.minor = element_blank(),
          # panel.grid.major = element_blank(),
          aspect.ratio = 1,
          legend.position = "none",
          plot.title = element_text(size = 8)
        )

      plots[[paste(f1, f2, sep = "_")]] <- p
    }
  }
  
  
  (
    plots[["negative_boundary"]]   | plots[["negative_complex"]]   | plots[["negative_negative"]]  ) /
    (plots[["complex_boundary"]] | plots[["complex_complex"]] | plots[["complex_negative"]]) /
    (plots[["boundary_boundary"]] | plots[["boundary_complex"]] | plots[["boundary_negative"]])
}

# ---- Load and classify all ecosystems ----
n_val <- 6
data <- read.csv(paste0("../data/boundary_distances/n_", n_val, ".csv"),
                 header = FALSE, sep = " ")
colnames(data) <- c("seed_i", "n", "d", "pert_size", "pert_i", "alpha_i", 
                    "deltap", "deltalin", "deltax", "deltaxlin", "flag")

# Classify each ecosystem (grouped by seed, n, pert_size, d)
classified <- data %>%
  group_by(seed_i, n, pert_size, d) %>%
  group_split() %>%
  map_dfr(function(group) {
    seed <- unique(group$seed_i)
    result_robust <- evaluate_response_type(group, "robustness")
    result_react  <- evaluate_response_type(group, "reactivity")

    tibble(
      seed = seed,
      p_robust = result_robust$proportion,
      p_reactive = result_react$proportion,
      ecosystem_class = case_when(
        p_robust > 0.5 & p_reactive > 0.5  ~ "robust_reactive",
        p_robust > 0.5 & p_reactive <= 0.5 ~ "robust_non_reactive",
        p_robust <= 0.5 & p_reactive <= 0.5 ~ "non_robust_non_reactive",
        TRUE ~ "non_robust_reactive"
      )
    )
  })

# ---- Prepare main data: normalize deltax ----
normalized_data <- data %>%
  select(seed_i, pert_i, alpha_i, deltap, deltax, flag) %>%
  group_by(seed_i, alpha_i) %>%
  mutate(deltax = deltax / max(deltax, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(seed = seed_i)

# ---- Merge classifications with normalized data ----
merged_data <- inner_join(normalized_data, classified, by = "seed")
merged_data <- merged_data %>% filter(flag != "nonconverged")


transitions_summary_df <- read_csv("transitions_summary.csv")
average_matrices <- transitions_summary_df %>%
  group_by(ecosystem_class, From, To) %>%
  summarise(Avg_Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

plot_avg_matrix <- function(df, class_name) {
  ggplot(df, aes(x = To, y = From, fill = Avg_Proportion)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    coord_fixed() +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = margin(4, 4, 4, 4, "pt"),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# ---- Plot by ecosystem_class (4 pages) ----
pdf("combined_histograms_by_ecosystem_class.pdf", width = 15, height = 7)

for (class_name in unique(merged_data$ecosystem_class)) {
  df_class <- merged_data %>% filter(ecosystem_class == class_name)

  flag_colors <- c(
    "negative"     = "#377EB8",
    "complex"      = "#984EA3",
    "boundary"     = "#4DAF4A",
    "nonconverged" = "#E41A1C"
  )

p_deltap = ggplot(df_class, aes(x = deltap)) +
  geom_histogram(data = subset(df_class, alpha_i == 0.1),
                 aes(fill = "α = 0.1"), bins = 30,
                 color = "black", fill = NA,
                 position = "identity") +
  geom_histogram(data = subset(df_class, alpha_i == 0.5),
                 aes(fill = "α = 0.5"), bins = 30,
                 color = "#E41A1C", fill = NA,
                 position = "identity") +
  facet_wrap(~ flag, labeller = label_both) +
  scale_fill_manual(values = c("α = 0.1" = "black", "α = 0.5" = "#E41A1C")) +
  labs(x = "Δx (normalized)", y = "Count", fill = "α") +
  theme_minimal()

p_deltax = ggplot(df_class, aes(x = deltax)) +
  geom_histogram(data = subset(df_class, alpha_i == 0.1),
                 aes(fill = "α = 0.1"), bins = 30,
                 color = "black", fill = NA,
                 position = "identity") +
  geom_histogram(data = subset(df_class, alpha_i == 0.5),
                 aes(fill = "α = 0.5"), bins = 30,
                 color = "#E41A1C", fill = NA,
                 position = "identity") +
  facet_wrap(~ flag, labeller = label_both) +
  scale_fill_manual(values = c("α = 0.1" = "black", "α = 0.5" = "#E41A1C")) +
  labs(x = "Δx (normalized)", y = "Count", fill = "α") +
  theme_minimal()
    

    # Compute average matrix for this class
avg_mat_class <- average_matrices %>% filter(ecosystem_class == class_name)
p_transition <- plot_avg_matrix(avg_mat_class, class_name)

# Combine and add inset
main_plot <- p_deltap / p_deltax +
  plot_annotation(title = paste("Ecosystem Class:", class_name))

full_plot <- main_plot +
  inset_element(p_transition, left = -.5, right = 0.7, top = 2.3, bottom = 1.6)

print(full_plot)
}

dev.off()



library(tidyverse)
library(patchwork)

# Read data as before...

# Wide format for flag and distances
wide <- data %>%
  select(seed_i, pert_i, alpha_i, deltap, deltax, flag) %>%
  pivot_wider(
    names_from = alpha_i,
    values_from = c(flag, deltap, deltax),
    names_sep = "_"
  )

# List all transitions found in the data
transitions <- wide %>%
  filter(flag_0.1 != "nonconverged", flag_0.5 != "nonconverged") %>%
  count(flag_0.1, flag_0.5) %>%
  filter(n > 0) %>%
  select(flag_0.1, flag_0.5)

plot_transition_histograms <- function(dist_name, title_prefix = "") {
    plots <- list()
    for (i in seq_len(nrow(transitions))) {
    flag_from <- transitions$flag_0.1[i]
    flag_to <- transitions$flag_0.5[i]
    
    subset <- wide %>% filter(flag_0.1 == flag_from, flag_0.5 == flag_to)
    if (nrow(subset) == 0) {
        p <- ggplot() + ggtitle(paste(flag_from, "→", flag_to, "\n(no data)"))
    } else {
        plot_data <- bind_rows(
        tibble(
            dist = subset[[paste0(dist_name, "_0.1")]],
            alpha = "α = 0.1"
        ),
        tibble(
            dist = subset[[paste0(dist_name, "_0.5")]],
            alpha = "α = 0.5"
        )
        )
        p <- ggplot(plot_data, aes(x = dist, fill = alpha, color = alpha)) +
        geom_histogram(position = "identity", bins = 30, alpha = 0.5) +
        scale_fill_manual(values = c("α = 0.1" = "black", "α = 0.5" = "#E41A1C")) +
        scale_color_manual(values = c("α = 0.1" = "black", "α = 0.5" = "#E41A1C")) +
        labs(
            title = paste(flag_from, "→", flag_to),
            x = title_prefix,
            y = "Count",
            fill = "α", color = "α"
        ) +
        theme_minimal(base_size = 11)
    }
    plots[[i]] <- p
    }

    # To display or save as a PDF:
    filename = paste0("transition_histograms_", dist_name, ".pdf")
    pdf(filename, width = 15, height = 10)
    print(wrap_plots(plots, ncol = 3) + plot_annotation(title = paste(title_prefix, "by transition (flag α=0.1 → α=0.5)")))
    dev.off()
}

# For parameter distances
dist_name <- "deltap"
title_prefix <- "Parameter distance (Δp)"
plot_transition_histograms(dist_name, title_prefix)

# For solution distances
dist_name <- "deltax"
title_prefix <- "Solution distance (Δx)"
plot_transition_histograms(dist_name, title_prefix)



library(ggplot2)
library(patchwork)
library(dplyr)

# Set up color scheme
flag_colors <- c(
  "negative" = "#377EB8",
  "complex"  = "#984EA3",
  "boundary" = "#4DAF4A"
)

plot_transition_histories <- function(wide, dist_name = "deltap", title_prefix = "Parameter distance (Δp)") {
  flag_colors <- c(
    "negative" = "#377EB8",
    "complex"  = "#984EA3",
    "boundary" = "#4DAF4A"
  )
  all_flags <- names(flag_colors)

  # Build each row: one boundary type, two plots (alpha = 0.1 and alpha = 0.5)
  row_plots <- lapply(all_flags, function(flag) {
    # α = 0.1 plot (left): All with flag_0.1 == flag, color by To (flag_0.5)
    plot_0.1 <- wide %>%
      filter(flag_0.1 == flag, flag_0.5 %in% all_flags) %>%
      select(flag_0.5, !!sym(paste0(dist_name, "_0.1"))) %>%
      rename(To = flag_0.5, dist = paste0(dist_name, "_0.1")) %>%
      mutate(To = factor(To, levels = all_flags)) %>%
      ggplot(aes(x = dist, fill = To, color = To)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.5) +
      scale_fill_manual(values = flag_colors) +
      scale_color_manual(values = flag_colors) +
      labs(
        title = paste("α = 0.1,", flag),
        x = title_prefix, y = "Count", fill = "To", color = "To"
      ) +
      theme_minimal(base_size = 12)

    # α = 0.5 plot (right): All with flag_0.5 == flag, color by From (flag_0.1)
    plot_0.5 <- wide %>%
      filter(flag_0.5 == flag, flag_0.1 %in% all_flags) %>%
      select(flag_0.1, !!sym(paste0(dist_name, "_0.5"))) %>%
      rename(From = flag_0.1, dist = paste0(dist_name, "_0.5")) %>%
      mutate(From = factor(From, levels = all_flags)) %>%
      ggplot(aes(x = dist, fill = From, color = From)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.5) +
      scale_fill_manual(values = flag_colors) +
      scale_color_manual(values = flag_colors) +
      labs(
        title = paste("α = 0.5,", flag),
        x = title_prefix, y = "Count", fill = "From", color = "From"
      ) +
      theme_minimal(base_size = 12)

    plot_0.1 | plot_0.5
  })

  final_plot <- row_plots[[1]] / row_plots[[2]] / row_plots[[3]] +
    plot_annotation(title = paste(title_prefix, ": 3 rows × 2 columns, α = 0.1 left, α = 0.5 right"))
  return(final_plot)
}

# Save for both metrics
pdf("transition_histories_deltap.pdf", width = 14, height = 15)
print(plot_transition_histories(wide, "deltap", "Parameter distance (Δp)"))
dev.off()

pdf("transition_histories_deltax.pdf", width = 14, height = 15)
print(plot_transition_histories(wide, "deltax", "Solution distance (Δx)"))
dev.off()



plot_transition_histories_v2 <- function(wide, dist_name = "deltap", title_prefix = "Parameter distance (Δp)") {
  flag_colors <- c(
    "negative" = "#377EB8",
    "complex"  = "#984EA3",
    "boundary" = "#4DAF4A"
  )
  all_flags <- names(flag_colors)
  
  # Top row: by From (flag_0.1), distances at alpha=0.1, color by To (flag_0.5)
  plots_top <- lapply(all_flags, function(from_flag) {
    plot_data <- wide %>%
      filter(flag_0.1 == from_flag, flag_0.5 %in% all_flags) %>%
      select(flag_0.5, !!sym(paste0(dist_name, "_0.1"))) %>%
      rename(To = flag_0.5, dist = paste0(dist_name, "_0.1")) %>%
      mutate(To = factor(To, levels = all_flags))
    ggplot(plot_data, aes(x = dist, fill = To, color = To)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.5) +
      scale_fill_manual(values = flag_colors) +
      scale_color_manual(values = flag_colors) +
      labs(
        title = paste("α = 0.1, From", from_flag),
        x = title_prefix, y = "Count", fill = "To", color = "To"
      ) +
      theme_minimal(base_size = 12)
  })
  
  # Bottom row: by To (flag_0.5), distances at alpha=0.5, color by From (flag_0.1)
  plots_bottom <- lapply(all_flags, function(to_flag) {
    plot_data <- wide %>%
      filter(flag_0.5 == to_flag, flag_0.1 %in% all_flags) %>%
      select(flag_0.1, !!sym(paste0(dist_name, "_0.5"))) %>%
      rename(From = flag_0.1, dist = paste0(dist_name, "_0.5")) %>%
      mutate(From = factor(From, levels = all_flags))
    ggplot(plot_data, aes(x = dist, fill = From, color = From)) +
      geom_histogram(position = "identity", bins = 30, alpha = 0.5) +
      scale_fill_manual(values = flag_colors) +
      scale_color_manual(values = flag_colors) +
      labs(
        title = paste("α = 0.5, To", to_flag),
        x = title_prefix, y = "Count", fill = "From", color = "From"
      ) +
      theme_minimal(base_size = 12)
  })
  
  final_plot <- (plots_top[[1]] | plots_top[[2]] | plots_top[[3]]) /
                (plots_bottom[[1]] | plots_bottom[[2]] | plots_bottom[[3]]) +
                plot_annotation(title = paste(title_prefix, ": Transitions by Boundary History (v2)"))
  return(final_plot)
}

# Save for both deltap and deltax
pdf("transition_histories_deltap_v2.pdf", width = 15, height = 10)
print(plot_transition_histories_v2(wide, "deltap", "Parameter distance (Δp)"))
dev.off()

pdf("transition_histories_deltax_v2.pdf", width = 15, height = 10)
print(plot_transition_histories_v2(wide, "deltax", "Solution distance (Δx)"))
dev.off()


# Filter and pivot all systems together
df_wide_all <- data %>%
  filter(flag != "nonconverged") %>%
  select(seed_i, pert_i, alpha_i, deltap, deltax, flag) %>%
  pivot_wider(
    names_from = alpha_i,
    values_from = c(deltap, deltax, flag),
    names_sep = "_"
  ) %>%
  mutate(
    delta_p_diff = deltap_0.5 - deltap_0.1,
    delta_x_diff = deltax_0.5 - deltax_0.1,
    flag_0.1 = factor(flag_0.1, levels = c("complex", "boundary", "negative")),
    flag_0.5 = factor(flag_0.5, levels = c("complex", "boundary", "negative"))
  ) %>%
  drop_na(delta_p_diff, delta_x_diff, flag_0.1, flag_0.5)


# Create histogram grids
p_grid_deltap_all <- make_hist_grid(df_wide_all, "delta_p_diff", "Δp differences")
p_grid_deltax_all <- make_hist_grid(df_wide_all, "delta_x_diff", "Δx differences")

grid_deltap_wrapped <- wrap_elements(full = p_grid_deltap_all)
grid_deltax_wrapped <- wrap_elements(full = p_grid_deltax_all)

# Stack side by side
right_stack_wrapped <- grid_deltap_wrapped | grid_deltax_wrapped + 
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "All systems — Δp and Δx differences by flag transition",
    theme = theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  )

ggsave("delta_diff_grid_all_systems.pdf", right_stack_wrapped, width = 14, height = 7)
