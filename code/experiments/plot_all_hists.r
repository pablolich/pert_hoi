library(tidyverse)

# Read and name columns
data <- read.csv("../../data/boundary_distances/n_3.csv", header = FALSE, sep = " ")
colnames(data) <- c("seed_i", "n", "d", "pert_size", "pert_i", "alpha_i", "deltap", "deltalin","deltax","deltaxlin","flag")

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

#now plot how many points in each simulation increase versus how many decrease.
data_wide <- data %>%
  mutate(alpha_label = paste0("alpha_", alpha_i)) %>%
  select(seed_i, n, d, pert_size, pert_i, alpha_label, deltap) %>%
  pivot_wider(
    names_from = alpha_label,
    values_from = deltap
  ) %>% 
  mutate(diff_delta = alpha_0.5 - alpha_0.1) %>% 
  filter(diff_delta != 0) %>% 
  mutate(robust = diff_delta > 0) %>% 
  ungroup() %>% 
  group_by(seed_i) %>% 
  summarise(av_robust = mean(robust))

#match with parameters and compare

dat_param = read.csv("../../data/all_parameter_sets_n3.csv")

u_thresh = 0.75
l_thresh = 0.25
# Merge data_wide with dat_param on seed_i
merged_data <- dat_param %>%
  left_join(data_wide, by = "seed_i") %>%
  drop_na() %>% 
  mutate(threshold_robust = ifelse(av_robust > u_thresh, 1, ifelse(av_robust < l_thresh, -1, 0))) %>% 
  filter(threshold_robust != 0)

ggplot(merged_data)+
  geom_histogram(aes(x = value, fill = as.factor(threshold_robust)),
                 position = "identity", alpha = 0.5)+
  facet_wrap(~order, ncol = 1)+
  xlim(-2,1.5)+
  geom_histogram(data = merged_data, aes(value),
                 alpha = 0.1) 
