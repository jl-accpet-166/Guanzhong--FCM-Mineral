# Libraries
library(MASS)
library(ggplot2)

# Data generation function
generate_traceability_data <- function(n_reference_per_source = 1500, n_test_points = 200, seed = 68) {
  set.seed(seed)
  
  source_parameters <- list(
    A = list(mu = c(1.0, 1.0), sigma = matrix(c(0.1, 0.05, 0.05, 0.1), nrow = 2)),
    B = list(mu = c(1.2, 1.1), sigma = matrix(c(0.1, 0.03, 0.03, 0.15), nrow = 2)),
    C = list(mu = c(2.5, 2.0), sigma = matrix(c(0.2, 0.1, 0.1, 0.2), nrow = 2))
  )
  
  reference_data <- data.frame()
  for (source_name in names(source_parameters)) {
    params <- source_parameters[[source_name]]
    source_data <- MASS::mvrnorm(n = n_reference_per_source, mu = params$mu, Sigma = params$sigma)
    ref_df <- data.frame(Sample_ID = paste0("Ref_", source_name, "_", 1:n_reference_per_source),
                         X = source_data[, 1], Y = source_data[, 2], Known_Source = source_name, Data_Type = "Reference", stringsAsFactors = FALSE)
    reference_data <- rbind(reference_data, ref_df)
  }
  
  test_points <- matrix(0, nrow = n_test_points, ncol = 2)
  true_sources <- character(n_test_points)
  
  n_A_B_overlap <- round(n_test_points * 0.3)
  n_A_only <- round(n_test_points * 0.2)
  n_B_only <- round(n_test_points * 0.2)
  n_C_only <- round(n_test_points * 0.3)
  
  overlap_center <- (source_parameters$A$mu + source_parameters$B$mu) / 2
  overlap_sigma <- (source_parameters$A$sigma + source_parameters$B$sigma) / 4
  test_points[1:n_A_B_overlap, ] <- MASS::mvrnorm(n = n_A_B_overlap, mu = overlap_center, Sigma = overlap_sigma)
  true_sources[1:n_A_B_overlap] <- "A_B_Overlap"
  
  A_center <- source_parameters$A$mu + (source_parameters$A$mu - source_parameters$B$mu) * 0.3
  start_idx <- n_A_B_overlap + 1
  end_idx <- n_A_B_overlap + n_A_only
  test_points[start_idx:end_idx, ] <- MASS::mvrnorm(n = n_A_only, mu = A_center, Sigma = source_parameters$A$sigma * 0.8)
  true_sources[start_idx:end_idx] <- "A"
  
  B_center <- source_parameters$B$mu + (source_parameters$B$mu - source_parameters$A$mu) * 0.3
  start_idx <- end_idx + 1
  end_idx <- end_idx + n_B_only
  test_points[start_idx:end_idx, ] <- MASS::mvrnorm(n = n_B_only, mu = B_center, Sigma = source_parameters$B$sigma * 0.8)
  true_sources[start_idx:end_idx] <- "B"
  
  start_idx <- end_idx + 1
  end_idx <- n_test_points
  test_points[start_idx:end_idx, ] <- MASS::mvrnorm(n = n_C_only, mu = source_parameters$C$mu, Sigma = source_parameters$C$sigma)
  true_sources[start_idx:end_idx] <- "C"
  
  test_data <- data.frame(Sample_ID = paste0("Test_", 1:n_test_points), X = test_points[, 1], Y = test_points[, 2], Data_Type = "Test", stringsAsFactors = FALSE)
  validation_info <- data.frame(Sample_ID = test_data$Sample_ID, True_Source = true_sources, stringsAsFactors = FALSE)
  
  list(reference_data = reference_data, test_data = test_data, validation_info = validation_info,
       parameters = list(reference_sources = c("A", "B", "C"), reference_size_per_source = n_reference_per_source,
                         test_distribution = c("A_B_Overlap" = n_A_B_overlap/n_test_points, "A" = n_A_only/n_test_points, "B" = n_B_only/n_test_points, "C" = n_C_only/n_test_points), random_seed = seed))
}

# Setup paths
save_path <- "./矿料"
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
scatter_plots_path <- file.path(save_path, "验证散点生成")
if (!dir.exists(scatter_plots_path)) dir.create(scatter_plots_path, recursive = TRUE)

# Generate data
trace_data <- generate_traceability_data(n_reference_per_source = 1500, n_test_points = 200, seed = 68)

# Split test data into groups
set.seed(68)
test_data <- trace_data$test_data
n_test <- nrow(test_data)
group_assignments <- sample(rep(1:4, length.out = n_test))
test_data$Group <- paste0("Group", group_assignments)

# Save data files
for (i in 1:4) {
  group_data <- test_data[test_data$Group == paste0("Group", i), ]
  write.csv(group_data, file.path(save_path, paste0("test_data_group", i, ".csv")), row.names = FALSE)
}
write.csv(trace_data$reference_data, file.path(save_path, "reference_data_known_sources.csv"), row.names = FALSE)
write.csv(test_data, file.path(save_path, "test_data_all_groups.csv"), row.names = FALSE)
write.csv(trace_data$validation_info, file.path(save_path, "validation_true_sources.csv"), row.names = FALSE)

# Theme setup
theme_custom <- theme_minimal() +
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        plot.margin = margin(20, 20, 20, 20))

# Reference data plot
p_reference <- ggplot(trace_data$reference_data, aes(x = X, y = Y, color = Known_Source)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c")) +
  labs(x = "Feature X", y = "Feature Y", color = "Source") +
  theme_custom

# Test data groups plot
p_test_groups <- ggplot(test_data, aes(x = X, y = Y, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Group1" = "#E41A1C", "Group2" = "#377EB8", "Group3" = "#4DAF4A", "Group4" = "#984EA3")) +
  labs(x = "Feature X", y = "Feature Y", color = "Group") +
  theme_custom

# Test data truth plot
test_with_truth <- merge(test_data, trace_data$validation_info, by = "Sample_ID")
p_test_truth <- ggplot(test_with_truth, aes(x = X, y = Y, color = True_Source)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c", "A_B_Overlap" = "#9467bd")) +
  labs(x = "Feature X", y = "Feature Y", color = "True Source") +
  theme_custom

# Combined plot
p_combined <- ggplot() +
  geom_point(data = trace_data$reference_data, aes(x = X, y = Y, color = Known_Source), alpha = 0.2, size = 1) +
  geom_point(data = test_data, aes(x = X, y = Y, color = Group), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c",
                                "Group1" = "#E41A1C", "Group2" = "#377EB8", "Group3" = "#4DAF4A", "Group4" = "#984EA3")) +
  labs(x = "Feature X", y = "Feature Y") +
  theme_custom

# Save plots
ggsave(file.path(scatter_plots_path, "reference_data_plot.png"), p_reference, width = 10, height = 8, dpi = 300)
ggsave(file.path(scatter_plots_path, "test_data_groups_plot.png"), p_test_groups, width = 10, height = 8, dpi = 300)
ggsave(file.path(scatter_plots_path, "test_data_truth_plot.png"), p_test_truth, width = 10, height = 8, dpi = 300)
ggsave(file.path(scatter_plots_path, "combined_data_plot.png"), p_combined, width = 10, height = 8, dpi = 300)