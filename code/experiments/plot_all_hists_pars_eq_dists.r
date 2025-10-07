library(tidyverse)
library(patchwork)

setwd("/home/srilenak/Desktop/pert_hoi/code")  # Set working directory to the parent directory

# ---- Function to classify robustness or reactivity ----
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

compute_flag_transition_matrix <- function(group) {
  # Ensure expected alpha values
  alpha_vals <- sort(unique(group$alpha_i))
  if (!all(c(0.1, 0.5) %in% alpha_vals)) return(NULL)

  # Create wide format data with flags
  df_wide <- group %>%
    select(pert_i, alpha_i, flag) %>%
    pivot_wider(names_from = alpha_i, values_from = flag, names_prefix = "alpha_")

  # Define all flag levels
  all_flags <- c("boundary", "complex", "negative", "nonconverged")

  # Cross-tabulate transitions
  transitions <- table(
    factor(df_wide$alpha_0.1, levels = all_flags),
    factor(df_wide$alpha_0.5, levels = all_flags)
  )

  # Convert to proportions
  transitions_prop <- transitions / sum(transitions)

  # Convert to data frame for plotting
  transition_df <- as.data.frame(as.table(transitions_prop))
  names(transition_df) <- c("From", "To", "Proportion")

  return(transition_df)
}

plot_transition_tile <- function(transition_df, seed, n, p) {
  # Filter out transitions involving nonconverged states
  filtered_df <- transition_df %>%
    filter(From != "nonconverged", To != "nonconverged")

  ggplot(filtered_df, aes(x = To, y = From, fill = Proportion)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "steelblue", guide = "none") +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_minimal(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(1, 1, 1, 1, "pt"),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    )
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


library(tidyverse)
library(patchwork)

plot_all_hists = function(data, data_param) {
  data_plot <- data %>%
    group_by(seed_i, n, pert_size, d) %>%
    group_split()

  all_transitions_df <- list()  # NEW

  pdf("all_histograms_by_seed.pdf", width = 15, height = 10)

  for (i in 1:50) {

        # --- Main plot for each group ---
    group = data_plot[[i]]
    seed <- unique(group$seed_i)
    n <- unique(group$n)
    p <- unique(group$pert_size)

    # --- Histograms of coefficients by order ---
    group_coeffs <- data_param %>%
      filter(seed_i == seed) %>%
      mutate(
        sign = ifelse(value < 0, "negative", "positive"),
        order_label = factor(order, levels = c(0, 1, 2), labels = c("r", "A", "B"))
      )

    coeffs_panel <- ggplot(group_coeffs) +
      geom_histogram(aes(x = value, fill = sign), bins = 30, color = "black") +
      scale_fill_manual(values = c("negative" = "#E41A1C", "positive" = "#4DAF4A")) +
      facet_wrap(~ order_label, scales = "free", nrow = 3) +
      theme_minimal(base_size = 9) +
      theme(legend.position = "none")

    # Classify robustness and reactivity
    result_robust <- evaluate_response_type(group, "robustness")
    result_react  <- evaluate_response_type(group, "reactivity")

    group$flag <- factor(group$flag, levels = c("boundary", "complex", "negative", "nonconverged"))
    flag_colors <- c(
      "negative"     = "#377EB8",
      "complex"      = "#984EA3",
      "boundary"     = "#4DAF4A",
      "nonconverged" = "#E41A1C"
    )

    # Define shapes based on flag types

    flag_shapes <- c(
      "complex" = 16,  # open circle
      "boundary" = 124,  # vertical bar |
      "negative" = 95    # horizontal bar —
    )


    p_deltap <- ggplot(group) +
      geom_histogram(aes(x = deltap, fill = flag),
                     bins = 20, color = "black", alpha = 0.7) +
      facet_grid(alpha_i ~ flag, labeller = label_both) +
      scale_fill_manual(values = flag_colors, drop = FALSE) +
      labs(
        title = paste(result_robust$classification,
                      sprintf("(%.2f)", result_robust$proportion)),
        x = "Δp", y = "Count", fill = "Flag"
      ) +
      theme_minimal()

    p_deltax <- ggplot(group) +
      geom_histogram(aes(x = deltax, fill = flag),
                     bins = 20, color = "black", alpha = 0.7) +
      facet_grid(alpha_i ~ flag, labeller = label_both) +
      scale_fill_manual(values = flag_colors, drop = FALSE) +
      labs(
        title = paste(result_react$classification,
                      sprintf("(%.2f)", result_react$proportion)),
        x = "Δx", y = "Count", fill = "Flag"
      ) +
      theme_minimal()

      # Prepare wide-format data for scatter plot
    df_wide <- group %>%
      filter(flag != "nonconverged") %>%
      select(pert_i, alpha_i, deltap, deltax, flag) %>%
      pivot_wider(names_from = alpha_i, values_from = c(deltap, deltax, flag),
                  names_sep = "_") %>%
      mutate(
        delta_p_diff = deltap_0.5 - deltap_0.1,
        delta_x_diff = deltax_0.5 - deltax_0.1,
        flag_0.1 = factor(flag_0.1, levels = names(flag_shapes)),
        flag_0.5 = factor(flag_0.5, levels = names(flag_colors))
      ) %>%
      drop_na(delta_p_diff, delta_x_diff, flag_0.1, flag_0.5)

    # Scatter plot
    p_scatter <- ggplot(df_wide, aes(x = delta_p_diff, y = delta_x_diff)) +
      geom_hline(yintercept = 0, color = "grey80", linetype = "dashed") +
      geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
      geom_point(aes(color = flag_0.5, shape = flag_0.1), size = 2, stroke = 0.6) +
      scale_color_manual(values = flag_colors, name = "Flag at α=0.5") +
      scale_shape_manual(values = flag_shapes, name = "Flag at α=0.1") +
      labs(
        x = expression(Delta * p[0.5] - Delta * p[0.1]),
        y = expression(Delta * x[0.5] - Delta * x[0.1]),
      ) +
      theme_minimal(base_size = 8) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank()
      )


    # Compute transition matrix and tile plot
    transition_df <- compute_flag_transition_matrix(group)
    # --- Annotate and collect transition_df ---
    transition_df <- transition_df %>%
      filter(From != "nonconverged", To != "nonconverged") %>%
      mutate(
        seed = seed,
        n = n,
        pert_size = p,
        d = unique(group$d),
        p_robust = result_robust$proportion,
        p_reactive = result_react$proportion
      )
    
    all_transitions_df[[length(all_transitions_df) + 1]] <- transition_df

    p_transition = plot_transition_tile(transition_df, seed, n, p)

    # Create the two 3x3 panels
    p_grid_deltap <- make_hist_grid(df_wide, "delta_p_diff", "Δp diff")
    p_grid_deltax <- make_hist_grid(df_wide, "delta_x_diff", "Δx diff")

    # Wrap the two 3x3 grids individually
    grid_deltap_wrapped <- wrap_elements(full = p_grid_deltap)
    grid_deltax_wrapped <- wrap_elements(full = p_grid_deltax)

    # Stack those (now atomic) panels into a vertical panel
    right_stack_wrapped <- grid_deltap_wrapped / grid_deltax_wrapped + 
    plot_layout(heights = c(1, 1))

    # Left panel with title
    left_panel <- (p_transition / p_scatter / coeffs_panel + 
                  plot_layout(heights = c(1, 1, 5)))

    # Right panel
    right_panel <- p_deltap / p_deltax + plot_layout(heights = c(1, 1))

    full_plot <- wrap_plots(left_panel, right_panel, right_stack_wrapped, widths = c(1, 2, 2))


    # Final layout with correct widths
    #full_plot <- wrap_plots(left_panel, right_panel, widths = c(1, 3.5))

    print(full_plot)
}

  dev.off()
  transitions_summary_df <- bind_rows(all_transitions_df)
  return(transitions_summary_df)
}

# Load data and run function
n_val = 6
data <- read.csv(paste0("../../data/boundary_distances/n_", n_val, ".csv"),
                 header = FALSE, sep = " ")
colnames(data) <- c("seed_i", "n", "d", "pert_size", "pert_i", "alpha_i", 
                    "deltap", "deltalin", "deltax", "deltaxlin", "flag")

#load parameters data

dat_param <- read.csv(paste0("../../data/all_parameter_sets_n", n_val, ".csv"))

data_small <-data %>% filter(seed_i < 100)
data_param_small <- data_param %>% filter(seed_i < 100) 
write_csv(data_small, "../../data/data_meeting/data_small.csv")
write_csv(data_param_small, "../../data/data_meeting/data_param_small.csv")

summary_df = plot_all_hists(data_small, data_param_small)


transitions_summary_df <- bind_rows(summary_df) %>%
  mutate(
    ecosystem_class = case_when(
      p_robust > 0.5 & p_reactive > 0.5  ~ "robust_reactive",
      p_robust > 0.5 & p_reactive <= 0.5 ~ "robust_non_reactive",
      p_robust <= 0.5 & p_reactive <= 0.5 ~ "non_robust_non_reactive",
      TRUE ~ "non_robust_reactive"
    )
  )

write_csv(transitions_summary_df, "transitions_summary.csv")

ecosystem_class_df <- transitions_summary_df %>%
  distinct(seed, n, pert_size, d, p_robust, p_reactive, ecosystem_class)

write_csv(ecosystem_class_df, "ecosystem_class_summary.csv")

library(tidyverse)
library(patchwork)

plot_scatter_grid_by_flags <- function(data_all) {
  # Define shapes and colors
  flag_colors <- c(
    "negative"     = "#377EB8",
    "complex"      = "#984EA3",
    "boundary"     = "#4DAF4A"
  )

  flag_shapes <- c(
    "complex"  = 16,  # closed circle
    "boundary" = 124, # vertical bar |
    "negative" = 95   # horizontal bar —
  )

  # Prepare wide-format data
  df_wide <- data_all %>%
    filter(flag != "nonconverged") %>%
    select(seed_i, pert_i, alpha_i, deltap, deltax, flag) %>%
    pivot_wider(names_from = alpha_i, values_from = c(deltap, deltax, flag),
                names_sep = "_") %>%
    mutate(
      delta_p_diff = deltap_0.5 - deltap_0.1,
      delta_x_diff = deltax_0.5 - deltax_0.1,
      flag_0.1 = factor(flag_0.1, levels = c("complex", "boundary", "negative")),
      flag_0.5 = factor(flag_0.5, levels = c("complex", "boundary", "negative"))
    ) %>%
    drop_na(delta_p_diff, delta_x_diff, flag_0.1, flag_0.5)

  # Global y limits
  y_min <- min(df_wide$delta_x_diff, na.rm = TRUE)
  y_max <- max(df_wide$delta_x_diff, na.rm = TRUE)

  # Generate plots for each combination
  plots <- list()
  for (f1 in levels(df_wide$flag_0.1)) {
    for (f2 in levels(df_wide$flag_0.5)) {
      df_sub <- df_wide %>%
        filter(flag_0.1 == f1, flag_0.5 == f2)

      p <- ggplot(df_sub, aes(x = delta_p_diff, y = delta_x_diff)) +
        geom_hline(yintercept = 0, color = "grey80", linetype = "dashed", linewidth = 0.3) +
        geom_vline(xintercept = 0, color = "grey80", linetype = "dashed", linewidth = 0.3) +
        geom_point(aes(color = flag_0.5, shape = flag_0.1), size = 1.5, stroke = 0.5) +
        scale_color_manual(values = flag_colors, drop = FALSE) +
        scale_shape_manual(values = flag_shapes, drop = FALSE) +
        coord_cartesian(xlim = c(-10, 10), ylim = c(y_min, y_max)) +
        theme_minimal(base_size = 8) +
        labs(
          title = paste("flag_0.1 =", f1, "| flag_0.5 =", f2),
          x = expression(Delta * p[0.5] - Delta * p[0.1]),
          y = expression(Delta * x[0.5] - Delta * x[0.1])
        ) +
        theme(
          legend.position = "none",
          plot.title = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text = element_text(size = 6),
          aspect.ratio = 1
        )

      plots[[paste(f1, f2, sep = "_")]] <- p
    }
  }

  # Assemble with patchwork into 3×3 grid
  p_grid <- (
    plots[["complex_complex"]]   | plots[["complex_boundary"]]   | plots[["complex_negative"]]  ) /
    (plots[["boundary_complex"]] | plots[["boundary_boundary"]] | plots[["boundary_negative"]]) /
    (plots[["negative_complex"]] | plots[["negative_boundary"]] | plots[["negative_negative"]])
  
  return(p_grid)
}


plot_out = plot_scatter_grid_by_flags(data)

# Save to PDF
ggsave("delta_diff_bin2d_grid.png", plot_out, width = 10, height = 10, units = "in")
