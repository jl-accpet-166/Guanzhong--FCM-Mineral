# Libraries
library(dplyr)
library(ggplot2)
library(MASS)
library(mvtnorm)
library(FNN)
library(tidyr)
library(gridExtra)
library(cluster)

# Setup paths
data_dir <- "./data"
validation_dir <- "./fcm_validation_results"
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)

# Read data
reference_data <- read.csv(file.path(data_dir, "reference_data_known_sources.csv"))
test_data_all <- read.csv(file.path(data_dir, "test_data_all_groups.csv"))
validation_info <- read.csv(file.path(data_dir, "validation_true_sources.csv"))

# Feature engineering function
enhanced_ABC_feature_engineering <- function(reference_data, test_data_all) {
  all_data <- rbind(
    reference_data[, c("X", "Y")] %>% mutate(Type = "Reference"),
    test_data_all[, c("X", "Y")] %>% mutate(Type = "Test")
  )
  
  median_vals <- apply(all_data[, c("X", "Y")], 2, median)
  iqr_vals <- apply(all_data[, c("X", "Y")], 2, IQR)
  iqr_vals[iqr_vals == 0] <- 1
  
  reference_data$X_std <- (reference_data$X - median_vals[1]) / iqr_vals[1]
  reference_data$Y_std <- (reference_data$Y - median_vals[2]) / iqr_vals[2]
  test_data_all$X_std <- (test_data_all$X - median_vals[1]) / iqr_vals[1]
  test_data_all$Y_std <- (test_data_all$Y - median_vals[2]) / iqr_vals[2]
  
  groups <- unique(test_data_all$Group)
  group_centers <- data.frame(Group = groups, Center_X = numeric(length(groups)), Center_Y = numeric(length(groups)))
  for (i in 1:length(groups)) {
    group_data <- test_data_all[test_data_all$Group == groups[i], ]
    group_centers$Center_X[i] <- median(group_data$X_std)
    group_centers$Center_Y[i] <- median(group_data$Y_std)
  }
  
  sources <- c("A", "B", "C")
  source_features <- list()
  
  for (source in sources) {
    source_data <- reference_data[reference_data$Known_Source == source, c("X_std", "Y_std")]
    n_samples <- nrow(source_data)
    
    if (n_samples > 5) {
      source_mean <- colMeans(source_data)
      source_cov <- cov(source_data)
      tryCatch({
        mcd_result <- cov.mcd(source_data, quantile.used = floor(0.75 * n_samples))
        robust_center <- mcd_result$center
        robust_cov <- mcd_result$cov
      }, error = function(e) {
        robust_center <- colMeans(source_data)
        robust_cov <- cov(source_data)
      })
      boundary_extension <- 0.22
      source_features[[source]] <- list(
        center = source_mean, robust_center = robust_center, cov = source_cov, robust_cov = robust_cov, n = n_samples,
        data_range = list(x = c(min(source_data$X_std) - boundary_extension, max(source_data$X_std) + boundary_extension),
                          y = c(min(source_data$Y_std) - boundary_extension, max(source_data$Y_std) + boundary_extension)), mcd_used = TRUE
      )
    } else if (n_samples > 1) {
      source_mean <- colMeans(source_data)
      source_cov <- cov(source_data)
      tryCatch({
        robust_cov <- cov.rob(source_data)$cov
        robust_center <- colMeans(source_data)
      }, error = function(e) {
        robust_cov <- source_cov
        robust_center <- source_mean
      })
      boundary_extension <- 0.25
      source_features[[source]] <- list(
        center = source_mean, robust_center = robust_center, cov = source_cov, robust_cov = robust_cov, n = n_samples,
        data_range = list(x = c(min(source_data$X_std) - boundary_extension, max(source_data$X_std) + boundary_extension),
                          y = c(min(source_data$Y_std) - boundary_extension, max(source_data$Y_std) + boundary_extension)), mcd_used = FALSE
      )
    } else {
      cov_factor <- 0.025
      boundary_extension <- 0.3
      source_features[[source]] <- list(
        center = as.numeric(source_data), robust_center = as.numeric(source_data),
        cov = diag(2) * cov_factor, robust_cov = diag(2) * cov_factor, n = 1,
        data_range = list(x = c(as.numeric(source_data[1]) - boundary_extension, as.numeric(source_data[1]) + boundary_extension),
                          y = c(as.numeric(source_data[2]) - boundary_extension, as.numeric(source_data[2]) + boundary_extension)), mcd_used = FALSE
      )
    }
  }
  
  test_data_all <- test_data_all %>%
    mutate(
      Dist_to_A = sqrt((X_std - source_features[["A"]]$robust_center[1])^2 + (Y_std - source_features[["A"]]$robust_center[2])^2),
      Dist_to_B = sqrt((X_std - source_features[["B"]]$robust_center[1])^2 + (Y_std - source_features[["B"]]$robust_center[2])^2),
      Dist_to_C = sqrt((X_std - source_features[["C"]]$robust_center[1])^2 + (Y_std - source_features[["C"]]$robust_center[2])^2)
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      Dist_ratio_AB = Dist_to_A / (Dist_to_B + 1e-6),
      Dist_ratio_AC = Dist_to_A / (Dist_to_C + 1e-6),
      Dist_ratio_BC = Dist_to_B / (Dist_to_C + 1e-6),
      Vector_A_X = X_std - source_features[["A"]]$robust_center[1],
      Vector_A_Y = Y_std - source_features[["A"]]$robust_center[2],
      Vector_B_X = X_std - source_features[["B"]]$robust_center[1],
      Vector_B_Y = Y_std - source_features[["B"]]$robust_center[2],
      Angle_AB = atan2(Vector_A_Y, Vector_A_X) - atan2(Vector_B_Y, Vector_B_X),
      Angle_AB_norm = (Angle_AB + pi) / (2 * pi)
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      Enhanced_Zone = case_when(
        X_std < -0.12 & Y_std > 0.25 ~ "A_Strong_Core",
        X_std < -0.05 & Y_std > 0.15 ~ "A_Moderate_Core",
        X_std < 0.02 & Y_std > 0.08 ~ "A_Extended_Zone",
        X_std > 0.25 & Y_std < -0.15 ~ "B_Strong_Core",
        X_std > 0.18 & Y_std < -0.08 ~ "B_Moderate_Core",
        X_std > 0.12 & Y_std < 0.02 ~ "B_Extended_Zone",
        (X_std > 0.35 & Y_std > 0.35) | (X_std < -0.25 & Y_std < -0.25) ~ "C_Strong_Core",
        (X_std > 0.28 & Y_std > 0.28) | (X_std < -0.18 & Y_std < -0.18) ~ "C_Moderate_Core",
        Dist_ratio_AB >= 0.85 & Dist_ratio_AB <= 1.18 & X_std >= -0.08 & X_std <= 0.15 & Y_std >= 0.02 & Y_std <= 0.18 ~ "AB_Core_Overlap_HighConfidence",
        Dist_ratio_AB >= 0.75 & Dist_ratio_AB <= 1.33 & X_std >= -0.15 & X_std <= 0.22 & Y_std >= -0.05 & Y_std <= 0.25 ~ "AB_Extended_Overlap_MediumConfidence",
        Dist_ratio_AB >= 0.65 & Dist_ratio_AB <= 1.54 & X_std >= -0.20 & X_std <= 0.28 & Y_std >= -0.10 & Y_std <= 0.30 ~ "AB_Potential_Overlap",
        Dist_ratio_AC >= 0.75 & Dist_ratio_AC <= 1.33 ~ "AC_Overlap_Zone",
        Dist_ratio_BC >= 0.75 & Dist_ratio_BC <= 1.33 ~ "BC_Overlap_Zone",
        TRUE ~ "Neutral_Zone"
      )
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      AB_Overlap_Score = case_when(
        Enhanced_Zone == "AB_Core_Overlap_HighConfidence" ~ 0.95,
        Enhanced_Zone == "AB_Extended_Overlap_MediumConfidence" ~ 0.80,
        Enhanced_Zone == "AB_Potential_Overlap" ~ 0.65,
        TRUE ~ 0
      ),
      Distance_Balance_AB = 1 - (abs(Dist_to_A - Dist_to_B) / (Dist_to_A + Dist_to_B + 1e-6)),
      Angle_Balance_AB = 1 - (abs(Angle_AB_norm - 0.5) * 2),
      Position_Centrality_AB = 1 - (sqrt((X_std - (source_features[["A"]]$robust_center[1] + source_features[["B"]]$robust_center[1])/2)^2 + 
                                         (Y_std - (source_features[["A"]]$robust_center[2] + source_features[["B"]]$robust_center[2])/2)^2) / 0.5),
      AB_Overlap_Confidence = (AB_Overlap_Score * 0.4 + Distance_Balance_AB * 0.25 + Angle_Balance_AB * 0.20 + Position_Centrality_AB * 0.15)
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      Is_AB_Overlap_HighConfidence = AB_Overlap_Confidence >= 0.85,
      Is_AB_Overlap_MediumConfidence = AB_Overlap_Confidence >= 0.70 & AB_Overlap_Confidence < 0.85,
      Is_AB_Overlap_Any = AB_Overlap_Confidence >= 0.60,
      AB_Overlap_Type = case_when(
        Is_AB_Overlap_HighConfidence ~ "High_Confidence_Overlap",
        Is_AB_Overlap_MediumConfidence ~ "Medium_Confidence_Overlap",
        Is_AB_Overlap_Any ~ "Low_Confidence_Overlap",
        TRUE ~ "Not_Overlap"
      )
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      A_Enhanced_Score = (1/(Dist_to_A + 0.05)) * 0.25 + (1 - pmin(Dist_to_A / 1.4, 1)) * 0.20 + exp(-Dist_ratio_AB) * 0.15 +
        (1 - abs(Angle_AB_norm - 0.5)) * 0.10 +
        case_when(Enhanced_Zone == "A_Strong_Core" ~ 0.20, Enhanced_Zone == "A_Moderate_Core" ~ 0.15, Enhanced_Zone == "A_Extended_Zone" ~ 0.08,
                  Enhanced_Zone %in% c("AB_Core_Overlap_HighConfidence", "AB_Extended_Overlap_MediumConfidence") ~ 0.12, TRUE ~ 0),
      B_Enhanced_Score = (1/(Dist_to_B + 0.05)) * 0.25 + (1 - pmin(Dist_to_B / 1.4, 1)) * 0.20 + exp(Dist_ratio_AB - 1) * 0.15 +
        abs(Angle_AB_norm - 0.5) * 0.10 +
        case_when(Enhanced_Zone == "B_Strong_Core" ~ 0.20, Enhanced_Zone == "B_Moderate_Core" ~ 0.15, Enhanced_Zone == "B_Extended_Zone" ~ 0.08,
                  Enhanced_Zone %in% c("AB_Core_Overlap_HighConfidence", "AB_Extended_Overlap_MediumConfidence") ~ 0.12, TRUE ~ 0),
      C_Enhanced_Score = (1/(Dist_to_C + 0.06)) * 0.30 + (1 - pmin(Dist_to_C / 1.8, 1)) * 0.25 + exp(-pmin(Dist_ratio_AC, Dist_ratio_BC)) * 0.15 +
        case_when(Enhanced_Zone == "C_Strong_Core" ~ 0.20, Enhanced_Zone == "C_Moderate_Core" ~ 0.15,
                  Enhanced_Zone %in% c("AC_Overlap_Zone", "BC_Overlap_Zone") ~ 0.10, TRUE ~ 0)
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      A_Protection_Enhanced = case_when(
        Enhanced_Zone == "A_Strong_Core" ~ 0.92, Enhanced_Zone == "A_Moderate_Core" ~ 0.78, Enhanced_Zone == "A_Extended_Zone" ~ 0.60,
        Enhanced_Zone %in% c("AB_Core_Overlap_HighConfidence", "AB_Extended_Overlap_MediumConfidence") ~ 0.45,
        Dist_ratio_AB < 0.7 ~ 0.35, TRUE ~ 0.15
      ),
      B_Protection_Enhanced = case_when(
        Enhanced_Zone == "B_Strong_Core" ~ 0.92, Enhanced_Zone == "B_Moderate_Core" ~ 0.78, Enhanced_Zone == "B_Extended_Zone" ~ 0.60,
        Enhanced_Zone %in% c("AB_Core_Overlap_HighConfidence", "AB_Extended_Overlap_MediumConfidence") ~ 0.45,
        Dist_ratio_AB > 1.4 ~ 0.35, TRUE ~ 0.15
      ),
      C_Protection_Enhanced = case_when(
        Enhanced_Zone == "C_Strong_Core" ~ 0.90, Enhanced_Zone == "C_Moderate_Core" ~ 0.75,
        Enhanced_Zone %in% c("AC_Overlap_Zone", "BC_Overlap_Zone") ~ 0.50,
        Dist_to_C < 0.3 ~ 0.40, TRUE ~ 0.15
      )
    )
  
  test_data_all <- test_data_all %>%
    mutate(
      AB_Score_Ratio_Enhanced = A_Enhanced_Score / (B_Enhanced_Score + 1e-6),
      AC_Score_Ratio_Enhanced = A_Enhanced_Score / (C_Enhanced_Score + 1e-6),
      BC_Score_Ratio_Enhanced = B_Enhanced_Score / (C_Enhanced_Score + 1e-6),
      Dominant_Source_Enhanced = case_when(
        A_Enhanced_Score > B_Enhanced_Score & A_Enhanced_Score > C_Enhanced_Score ~ "A",
        B_Enhanced_Score > A_Enhanced_Score & B_Enhanced_Score > C_Enhanced_Score ~ "B",
        C_Enhanced_Score > A_Enhanced_Score & C_Enhanced_Score > B_Enhanced_Score ~ "C",
        TRUE ~ "Tie"
      ),
      Max_Enhanced_Score = pmax(A_Enhanced_Score, B_Enhanced_Score, C_Enhanced_Score),
      Score_Uncertainty = (A_Enhanced_Score + B_Enhanced_Score + C_Enhanced_Score) / (3 * Max_Enhanced_Score + 1e-6)
    )
  
  return(list(reference_data = reference_data, test_data = test_data_all, group_centers = group_centers,
              source_features = source_features, norm_params = list(median = median_vals, iqr = iqr_vals)))
}

# FCM algorithm
enhanced_ABC_fcm_algorithm <- function(reference_data, test_data, sources, group_centers, source_features) {
  results <- data.frame(Sample_ID = test_data$Sample_ID, Group = test_data$Group, X = test_data$X, Y = test_data$Y,
                        X_std = test_data$X_std, Y_std = test_data$Y_std)
  
  max_iterations <- 120; tolerance <- 1e-7; fuzziness <- 1.4; movement_limit <- 0.0012; base_regularization <- 0.016
  
  source_specific_params <- list(
    A = list(mahalanobis_weight = 0.58, euclidean_weight = 0.42, confidence_boost = 1.16, regularization_factor = 0.35),
    B = list(mahalanobis_weight = 0.52, euclidean_weight = 0.48, confidence_boost = 1.16, regularization_factor = 0.35),
    C = list(mahalanobis_weight = 0.60, euclidean_weight = 0.40, confidence_boost = 1.18, regularization_factor = 0.30)
  )
  
  for (source in sources) results[[source]] <- 0
  groups <- unique(test_data$Group)
  
  for (group in groups) {
    group_test_data <- test_data[test_data$Group == group, ]
    if (nrow(group_test_data) == 0) next
    group_indices <- which(test_data$Group == group)
    test_points <- as.matrix(group_test_data[, c("X_std", "Y_std")])
    
    for (source in sources) {
      params <- source_specific_params[[source]]
      source_center <- source_features[[source]]$robust_center
      group_center_row <- group_centers[group_centers$Group == group, ]
      group_center <- c(group_center_row$Center_X, group_center_row$Center_Y)
      group_source_dist <- sqrt(sum((group_center - source_center)^2))
      similarity_weight <- exp(-group_source_dist / 0.30)
      adjusted_center <- 0.80 * source_center + (1 - 0.80) * group_center * similarity_weight
      current_center <- adjusted_center
      
      source_cov_reg <- tryCatch({
        eigen_vals <- eigen(source_features[[source]]$robust_cov)$values
        min_eigen <- min(eigen_vals)
        reg_factor <- base_regularization * params$regularization_factor
        if (min_eigen < 1e-6) source_features[[source]]$robust_cov + diag(2) * (reg_factor + 0.010)
        else source_features[[source]]$robust_cov + diag(2) * reg_factor
      }, error = function(e) diag(2) * 0.032)
      
      current_probs <- rep(45, nrow(test_points))
      
      for (iter in 1:max_iterations) {
        mahalanobis_dist <- tryCatch({
          md <- mahalanobis(test_points, current_center, source_cov_reg)
          md[is.infinite(md) | md < 0] <- 1e-6
          md
        }, error = function(e) sqrt(rowSums((test_points - matrix(current_center, nrow = nrow(test_points), ncol = 2, byrow = TRUE))^2)))
        
        euclidean_dist <- sqrt(rowSums((test_points - matrix(current_center, nrow = nrow(test_points), ncol = 2, byrow = TRUE))^2))
        euclidean_dist[is.infinite(euclidean_dist) | euclidean_dist < 0] <- 1e-6
        
        combined_dist <- params$mahalanobis_weight * mahalanobis_dist + params$euclidean_weight * euclidean_dist
        if (all(is.na(combined_dist))) combined_dist <- rep(1, length(combined_dist))
        
        scale_factor <- max(combined_dist, na.rm = TRUE) / 1.65
        if (scale_factor <= 0) scale_factor <- 1
        new_probs <- 100 * exp(-combined_dist / scale_factor)
        
        boost_indices <- if (source == "A") {
          which(test_data$A_Protection_Enhanced[group_indices] > 0.75 & 
                test_data$A_Enhanced_Score[group_indices] > test_data$B_Enhanced_Score[group_indices] * 1.25 &
                test_data$A_Enhanced_Score[group_indices] > test_data$C_Enhanced_Score[group_indices] * 1.25)
        } else if (source == "B") {
          which(test_data$B_Protection_Enhanced[group_indices] > 0.75 & 
                test_data$B_Enhanced_Score[group_indices] > test_data$A_Enhanced_Score[group_indices] * 1.25 &
                test_data$B_Enhanced_Score[group_indices] > test_data$C_Enhanced_Score[group_indices] * 1.25)
        } else {
          which(test_data$C_Protection_Enhanced[group_indices] > 0.75 & 
                test_data$C_Enhanced_Score[group_indices] > test_data$A_Enhanced_Score[group_indices] * 1.3 &
                test_data$C_Enhanced_Score[group_indices] > test_data$B_Enhanced_Score[group_indices] * 1.3)
        }
        
        if (length(boost_indices) > 0) new_probs[boost_indices] <- pmin(new_probs[boost_indices] * params$confidence_boost, 96)
        
        prob_change <- mean(abs(new_probs - current_probs))
        if (prob_change < tolerance && iter > 20) break
        current_probs <- new_probs
        
        weights <- (current_probs/100)^fuzziness
        total_weight <- sum(weights)
        if (total_weight > 0) {
          new_center <- colSums(test_points * weights) / total_weight
          center_diff <- new_center - current_center
          move_dist <- sqrt(sum(center_diff^2))
          if (move_dist > movement_limit) new_center <- current_center + center_diff * (movement_limit / move_dist)
          current_center <- 0.85 * current_center + 0.15 * new_center
        }
      }
      results[group_indices, source] <- current_probs
    }
  }
  
  source_probs <- results[, sources]
  results$Predicted_Source <- apply(source_probs, 1, function(x) { if (all(is.na(x)) || all(x == 0)) return(NA); sources[which.max(x)] })
  results$Max_Probability <- apply(source_probs, 1, function(x) { if (all(is.na(x)) || all(x == 0)) return(45); max(x, na.rm = TRUE) })
  results$Prob_Diff_AB <- abs(results$A - results$B)
  results$Prob_Diff_AC <- abs(results$A - results$C)
  results$Prob_Diff_BC <- abs(results$B - results$C)
  results$Second_Probability <- apply(source_probs, 1, function(x) {
    if (all(is.na(x)) || all(x == 0)) return(40)
    sorted_probs <- sort(x, decreasing = TRUE)
    if (length(sorted_probs) >= 2) sorted_probs[2] else sorted_probs[1]
  })
  results$Probability_Ratio <- results$Max_Probability / (results$Second_Probability + 1e-6)
  results$Confidence_Level <- ifelse(results$Max_Probability >= 82 & results$Probability_Ratio >= 2.0, "High",
                                     ifelse(results$Max_Probability >= 68 & results$Probability_Ratio >= 1.6, "Medium", "Low"))
  return(results)
}

# Post-processing function
enhanced_ABC_post_processing <- function(results, test_data, source_features, validation_info) {
  sources <- c("A", "B", "C")
  results$Optimized_Prediction <- results$Predicted_Source
  results$Is_True_AB_Overlap <- FALSE
  results$Final_Prediction <- results$Predicted_Source
  results$Confidence <- results$Confidence_Level
  
  high_confidence_overlap_indices <- which(test_data$Is_AB_Overlap_HighConfidence & test_data$AB_Overlap_Confidence >= 0.90 &
                                           test_data$Distance_Balance_AB >= 0.85 & test_data$Position_Centrality_AB >= 0.80)
  results$Is_True_AB_Overlap[high_confidence_overlap_indices] <- TRUE
  
  probability_balance_indices <- which(abs(results$A - results$B) < 15 & results$Max_Probability < 75 &
                                       results$Probability_Ratio < 1.8 & test_data$AB_Overlap_Confidence >= 0.75)
  results$Is_True_AB_Overlap[probability_balance_indices] <- TRUE
  
  overlap_indices <- which(results$Is_True_AB_Overlap)
  if (length(overlap_indices) > 0) {
    results$Optimized_Prediction[overlap_indices] <- "A_B_Overlap"
    results$Final_Prediction[overlap_indices] <- "A_B_Overlap"
    for (idx in overlap_indices) {
      overlap_confidence <- test_data$AB_Overlap_Confidence[idx]
      if (overlap_confidence >= 0.90) { results$Confidence[idx] <- "Very_High"; results$Max_Probability[idx] <- 95 }
      else if (overlap_confidence >= 0.75) { results$Confidence[idx] <- "High"; results$Max_Probability[idx] <- 85 }
      else { results$Confidence[idx] <- "Medium"; results$Max_Probability[idx] <- 75 }
    }
  }
  
  results$A_Misclassified_Enhanced <- FALSE; results$B_Misclassified_Enhanced <- FALSE; results$C_Misclassified_Enhanced <- FALSE
  results$Is_AB_Borderline_Enhanced <- FALSE; results$Is_AC_Borderline_Enhanced <- FALSE; results$Is_BC_Borderline_Enhanced <- FALSE
  
  a_misclassified <- results$Predicted_Source != "A" & !results$Is_True_AB_Overlap &
    (test_data$A_Protection_Enhanced > 0.78 | test_data$A_Enhanced_Score > pmax(test_data$B_Enhanced_Score, test_data$C_Enhanced_Score) * 1.3 |
     test_data$Enhanced_Zone %in% c("A_Strong_Core", "A_Moderate_Core") | test_data$Dist_ratio_AB < 0.65)
  
  b_misclassified <- results$Predicted_Source != "B" & !results$Is_True_AB_Overlap &
    (test_data$B_Protection_Enhanced > 0.78 | test_data$B_Enhanced_Score > pmax(test_data$A_Enhanced_Score, test_data$C_Enhanced_Score) * 1.3 |
     test_data$Enhanced_Zone %in% c("B_Strong_Core", "B_Moderate_Core") | test_data$Dist_ratio_AB > 1.5)
  
  c_misclassified <- results$Predicted_Source != "C" & !results$Is_True_AB_Overlap &
    (test_data$C_Protection_Enhanced > 0.75 | test_data$C_Enhanced_Score > pmax(test_data$A_Enhanced_Score, test_data$B_Enhanced_Score) * 1.35 |
     test_data$Enhanced_Zone %in% c("C_Strong_Core", "C_Moderate_Core") | test_data$Dist_to_C < 0.25)
  
  ab_borderline <- (results$Prob_Diff_AB < 20 | results$Probability_Ratio < 1.75) & results$Max_Probability < 80 &
    test_data$Enhanced_Zone %in% c("AB_Core_Overlap_HighConfidence", "AB_Extended_Overlap_MediumConfidence") & !results$Is_True_AB_Overlap
  ac_borderline <- (results$Prob_Diff_AC < 20 | results$Probability_Ratio < 1.75) & results$Max_Probability < 80 & test_data$Enhanced_Zone == "AC_Overlap_Zone"
  bc_borderline <- (results$Prob_Diff_BC < 20 | results$Probability_Ratio < 1.75) & results$Max_Probability < 80 & test_data$Enhanced_Zone == "BC_Overlap_Zone"
  
  if (length(a_misclassified) == nrow(results)) results$A_Misclassified_Enhanced <- a_misclassified
  if (length(b_misclassified) == nrow(results)) results$B_Misclassified_Enhanced <- b_misclassified
  if (length(c_misclassified) == nrow(results)) results$C_Misclassified_Enhanced <- c_misclassified
  if (length(ab_borderline) == nrow(results)) results$Is_AB_Borderline_Enhanced <- ab_borderline
  if (length(ac_borderline) == nrow(results)) results$Is_AC_Borderline_Enhanced <- ac_borderline
  if (length(bc_borderline) == nrow(results)) results$Is_BC_Borderline_Enhanced <- bc_borderline
  
  # A misclassification correction
  for (idx in which(results$A_Misclassified_Enhanced)) {
    correction_confidence <- 0
    if (test_data$A_Protection_Enhanced[idx] > 0.78) correction_confidence <- correction_confidence + 0.30
    if (test_data$A_Enhanced_Score[idx] > test_data$B_Enhanced_Score[idx] * 1.3) correction_confidence <- correction_confidence + 0.25
    if (test_data$A_Enhanced_Score[idx] > test_data$C_Enhanced_Score[idx] * 1.3) correction_confidence <- correction_confidence + 0.25
    if (test_data$Dist_ratio_AB[idx] < 0.65) correction_confidence <- correction_confidence + 0.20
    if (correction_confidence >= 0.75) results$Optimized_Prediction[idx] <- "A"
  }
  
  # B misclassification correction
  for (idx in which(results$B_Misclassified_Enhanced)) {
    correction_confidence <- 0
    if (test_data$B_Protection_Enhanced[idx] > 0.78) correction_confidence <- correction_confidence + 0.30
    if (test_data$B_Enhanced_Score[idx] > test_data$A_Enhanced_Score[idx] * 1.3) correction_confidence <- correction_confidence + 0.25
    if (test_data$B_Enhanced_Score[idx] > test_data$C_Enhanced_Score[idx] * 1.3) correction_confidence <- correction_confidence + 0.25
    if (test_data$Dist_ratio_AB[idx] > 1.5) correction_confidence <- correction_confidence + 0.20
    if (correction_confidence >= 0.75) results$Optimized_Prediction[idx] <- "B"
  }
  
  # C misclassification correction
  for (idx in which(results$C_Misclassified_Enhanced)) {
    correction_confidence <- 0
    if (test_data$C_Protection_Enhanced[idx] > 0.75) correction_confidence <- correction_confidence + 0.35
    if (test_data$C_Enhanced_Score[idx] > test_data$A_Enhanced_Score[idx] * 1.35) correction_confidence <- correction_confidence + 0.30
    if (test_data$C_Enhanced_Score[idx] > test_data$B_Enhanced_Score[idx] * 1.35) correction_confidence <- correction_confidence + 0.30
    if (test_data$Dist_to_C[idx] < 0.25) correction_confidence <- correction_confidence + 0.25
    if (correction_confidence >= 0.85) results$Optimized_Prediction[idx] <- "C"
  }
  
  # Borderline corrections
  for (idx in which(results$Is_AB_Borderline_Enhanced)) {
    current_pred <- results$Optimized_Prediction[idx]
    if (test_data$A_Enhanced_Score[idx] > test_data$B_Enhanced_Score[idx] * 1.25 && current_pred != "A") results$Optimized_Prediction[idx] <- "A"
    else if (test_data$B_Enhanced_Score[idx] > test_data$A_Enhanced_Score[idx] * 1.25 && current_pred != "B") results$Optimized_Prediction[idx] <- "B"
  }
  
  for (idx in which(results$Is_AC_Borderline_Enhanced)) {
    current_pred <- results$Optimized_Prediction[idx]
    if (test_data$A_Enhanced_Score[idx] > test_data$C_Enhanced_Score[idx] * 1.3 && current_pred != "A") results$Optimized_Prediction[idx] <- "A"
    else if (test_data$C_Enhanced_Score[idx] > test_data$A_Enhanced_Score[idx] * 1.3 && current_pred != "C") results$Optimized_Prediction[idx] <- "C"
  }
  
  for (idx in which(results$Is_BC_Borderline_Enhanced)) {
    current_pred <- results$Optimized_Prediction[idx]
    if (test_data$B_Enhanced_Score[idx] > test_data$C_Enhanced_Score[idx] * 1.3 && current_pred != "B") results$Optimized_Prediction[idx] <- "B"
    else if (test_data$C_Enhanced_Score[idx] > test_data$B_Enhanced_Score[idx] * 1.3 && current_pred != "C") results$Optimized_Prediction[idx] <- "C"
  }
  
  results$Final_Prediction <- results$Optimized_Prediction
  results$Confidence <- ifelse(results$Max_Probability >= 83 & results$Probability_Ratio >= 2.1, "High",
                               ifelse(results$Max_Probability >= 69 & results$Probability_Ratio >= 1.7, "Medium", "Low"))
  
  results$Sample_Type <- "Normal"
  results$Sample_Type[results$Is_True_AB_Overlap] <- "AB_True_Overlap"
  results$Sample_Type[results$A_Misclassified_Enhanced] <- "A_Misclassified"
  results$Sample_Type[results$B_Misclassified_Enhanced] <- "B_Misclassified"
  results$Sample_Type[results$C_Misclassified_Enhanced] <- "C_Misclassified"
  results$Sample_Type[results$Is_AB_Borderline_Enhanced] <- "AB_Borderline"
  results$Sample_Type[results$Is_AC_Borderline_Enhanced] <- "AC_Borderline"
  results$Sample_Type[results$Is_BC_Borderline_Enhanced] <- "BC_Borderline"
  
  return(results)
}

# Performance evaluation
performance_evaluation <- function(predictions, true_labels, test_data) {
  evaluation_data <- merge(predictions, true_labels, by = "Sample_ID")
  
  evaluation_data$is_correct <- case_when(
    evaluation_data$True_Source == "A_B_Overlap" & evaluation_data$Final_Prediction %in% c("A", "B") ~ TRUE,
    evaluation_data$True_Source != "A_B_Overlap" & evaluation_data$Final_Prediction == evaluation_data$True_Source ~ TRUE,
    TRUE ~ FALSE
  )
  
  overlap_data <- evaluation_data[evaluation_data$True_Source == "A_B_Overlap", ]
  overlap_accuracy <- if(nrow(overlap_data) > 0) mean(overlap_data$is_correct) else NA
  
  total_points <- nrow(evaluation_data)
  correct_points <- sum(evaluation_data$is_correct, na.rm = TRUE)
  accuracy_overall <- mean(evaluation_data$is_correct, na.rm = TRUE)
  
  unique_sources <- unique(evaluation_data$True_Source)
  class_accuracy <- numeric(length(unique_sources))
  names(class_accuracy) <- unique_sources
  for (source in unique_sources) {
    source_data <- evaluation_data[evaluation_data$True_Source == source, ]
    if (nrow(source_data) == 0) { class_accuracy[source] <- NA; next }
    if (source == "A_B_Overlap") class_accuracy[source] <- mean(source_data$is_correct, na.rm = TRUE)
    else class_accuracy[source] <- mean(source_data$Final_Prediction == source, na.rm = TRUE)
  }
  
  sample_type_accuracy <- list()
  sample_types <- c("AB_True_Overlap", "A_Misclassified", "B_Misclassified", "C_Misclassified", "AB_Borderline", "AC_Borderline", "BC_Borderline", "Normal")
  for (type in sample_types) {
    type_data <- evaluation_data[evaluation_data$Sample_Type == type, ]
    sample_type_accuracy[[type]] <- if (nrow(type_data) > 0) mean(type_data$is_correct, na.rm = TRUE) else NA
  }
  
  confusion_matrix <- table(True = evaluation_data$True_Source, Predicted = evaluation_data$Final_Prediction)
  
  return(list(accuracy_overall = accuracy_overall, class_accuracy = class_accuracy,
              overlap_performance = list(total_overlap_samples = nrow(overlap_data), overlap_accuracy = overlap_accuracy),
              sample_type_accuracy = sample_type_accuracy, confusion_matrix = confusion_matrix, evaluation_data = evaluation_data,
              correct_points_ratio = correct_points / total_points, total_points = total_points, correct_points = correct_points))
}

# Visualization function
create_two_enhanced_plots <- function(evaluation_data, performance_results, save_dir, reference_data, test_data_all, validation_info) {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  theme_custom <- theme_minimal() +
    theme(text = element_text(family = "Times New Roman", face = "bold", size = 14),
          axis.title = element_text(size = 16, face = "bold"), axis.text = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"), legend.text = element_text(size = 14, face = "bold"),
          plot.margin = margin(20, 20, 20, 20), plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12))
  
  ref_x_col <- ifelse("X" %in% colnames(reference_data), "X", ifelse("X_std" %in% colnames(reference_data), "X_std", colnames(reference_data)[grepl("^X|^x", colnames(reference_data))][1]))
  ref_y_col <- ifelse("Y" %in% colnames(reference_data), "Y", ifelse("Y_std" %in% colnames(reference_data), "Y_std", colnames(reference_data)[grepl("^Y|^y", colnames(reference_data))][1]))
  
  p_reference <- ggplot(reference_data, aes_string(x = ref_x_col, y = ref_y_col, color = "Known_Source")) +
    geom_point(alpha = 0.7, size = 1) + scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c")) +
    labs(x = "Feature X", y = "Feature Y", color = "Known Source") + theme_custom
  
  test_data_plot <- test_data_all
  validation_id_col <- ifelse("Sample_ID" %in% colnames(validation_info), "Sample_ID", colnames(validation_info)[grepl("ID|id|Id", colnames(validation_info))][1])
  if (!is.na(validation_id_col) && validation_id_col != "Sample_ID") validation_info$Sample_ID <- validation_info[[validation_id_col]]
  test_data_with_truth <- merge(test_data_plot[, c("Sample_ID", "Group", "X", "Y")], validation_info[, c("Sample_ID", "True_Source")], by = "Sample_ID", all.x = TRUE)
  
  p_test_groups <- ggplot(test_data_with_truth, aes(x = X, y = Y, color = Group)) +
    geom_point(alpha = 0.7, size = 2) + scale_color_manual(values = c("Group1" = "#E41A1C", "Group2" = "#377EB8", "Group3" = "#4DAF4A", "Group4" = "#984EA3")) +
    labs(x = "Feature X", y = "Feature Y", color = "Group") + theme_custom
  
  p_test_truth <- ggplot(test_data_with_truth, aes(x = X, y = Y, color = True_Source)) +
    geom_point(alpha = 0.7, size = 2) + scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c", "A_B_Overlap" = "#9467bd", "Unknown" = "#7f7f7f")) +
    labs(x = "Feature X", y = "Feature Y", color = "True Source") + theme_custom
  
  ref_plot_data <- reference_data; ref_plot_data$X <- ref_plot_data[[ref_x_col]]; ref_plot_data$Y <- ref_plot_data[[ref_y_col]]
  p_combined <- ggplot() +
    geom_point(data = ref_plot_data, aes(x = X, y = Y, color = Known_Source), alpha = 0.2, size = 1) +
    geom_point(data = test_data_with_truth, aes(x = X, y = Y, color = Group), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("A" = "#1f77b4", "B" = "#ff7f0e", "C" = "#2ca02c", "Group1" = "#E41A1C", "Group2" = "#377EB8", "Group3" = "#4DAF4A", "Group4" = "#984EA3")) +
    labs(x = "Feature X", y = "Feature Y", color = "Source/Group") + theme_custom
  
  p_prediction <- ggplot(evaluation_data, aes(x = X, y = Y, color = Final_Prediction)) +
    geom_point(size = 2.5, alpha = 0.8) + scale_color_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "C" = "#4DAF4A", "A_B_Overlap" = "#984EA3")) +
    labs(x = "X", y = "Y", color = "Predicted Source") + theme_custom
  
  p_correctness <- ggplot(evaluation_data, aes(x = X, y = Y, color = is_correct, shape = Final_Prediction)) +
    geom_point(size = 2.5, alpha = 0.8) + scale_color_manual(values = c("TRUE" = "#2E8B57", "FALSE" = "#DC143C"), labels = c("TRUE" = "True", "FALSE" = "False"), name = "Validity") +
    scale_shape_manual(values = c("A" = 16, "B" = 17, "C" = 15, "A_B_Overlap" = 18), name = "Source") +
    labs(subtitle = sprintf("Correct Sample: %d/%d (%.2f%%)", performance_results$correct_points, performance_results$total_points, performance_results$accuracy_overall * 100), x = "X", y = "Y") + theme_custom
  
  p_overlap <- ggplot(evaluation_data, aes(x = X, y = Y, color = Final_Prediction, shape = True_Source == "A_B_Overlap")) +
    geom_point(size = 2.5, alpha = 0.8) + scale_color_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "C" = "#4DAF4A", "A_B_Overlap" = "#984EA3")) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16), labels = c("TRUE" = "True Overlap", "FALSE" = "Non-overlap"), name = "True Type") +
    labs(subtitle = sprintf("Overlap Accuracy: %.2f%%, Total Overlap: %d", performance_results$overlap_performance$overlap_accuracy * 100, performance_results$overlap_performance$total_overlap_samples), x = "X", y = "Y") + theme_custom
  
  scatter_plots_path <- file.path(save_dir, "validation_scatter_plots")
  if (!dir.exists(scatter_plots_path)) dir.create(scatter_plots_path, recursive = TRUE)
  
  ggsave(file.path(scatter_plots_path, "plot1_reference_known_source.png"), p_reference, width = 10, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot2_test_groups.png"), p_test_groups, width = 10, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot3_test_true_source.png"), p_test_truth, width = 10, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot4_combined.png"), p_combined, width = 10, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot5_prediction.png"), p_prediction, width = 12, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot6_correctness.png"), p_correctness, width = 12, height = 8, dpi = 300)
  ggsave(file.path(scatter_plots_path, "plot7_overlap.png"), p_overlap, width = 12, height = 8, dpi = 300)
  
  return(list(reference_plot = p_reference, test_groups_plot = p_test_groups, test_truth_plot = p_test_truth, combined_plot = p_combined,
              prediction_plot = p_prediction, correctness_plot = p_correctness, overlap_plot = p_overlap))
}

# Group proportion analysis
group_proportion_analysis <- function(evaluation_data, save_dir) {
  predicted_proportions <- evaluation_data %>%
    group_by(Group, Final_Prediction) %>% summarise(Count = n(), .groups = 'drop') %>%
    group_by(Group) %>% mutate(Total_in_Group = sum(Count), Predicted_Proportion = Count / Total_in_Group * 100) %>% ungroup()
  
  true_proportions <- evaluation_data %>%
    group_by(Group, True_Source) %>% summarise(Count = n(), .groups = 'drop') %>%
    group_by(Group) %>% mutate(Total_in_Group = sum(Count), True_Proportion = Count / Total_in_Group * 100) %>% ungroup()
  
  proportion_comparison <- full_join(predicted_proportions, true_proportions, by = c("Group", "Final_Prediction" = "True_Source"), suffix = c("_pred", "_true")) %>%
    rename(Source = Final_Prediction, Predicted_Count = Count_pred, True_Count = Count_true, Predicted_Percentage = Predicted_Proportion, True_Percentage = True_Proportion)
  
  accuracy_by_group <- evaluation_data %>%
    group_by(Group) %>% summarise(Total_Samples = n(), Correct_Predictions = sum(is_correct), Group_Accuracy = Correct_Predictions / Total_Samples * 100, .groups = 'drop')
  
  write.csv(predicted_proportions, file.path(save_dir, "group_predicted_proportions.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(true_proportions, file.path(save_dir, "group_true_proportions.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(proportion_comparison, file.path(save_dir, "group_proportion_comparison.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(accuracy_by_group, file.path(save_dir, "group_accuracy.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  
  p_pred <- ggplot(predicted_proportions, aes(x = Group, y = Predicted_Proportion, fill = Final_Prediction)) +
    geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "C" = "#4DAF4A", "A_B_Overlap" = "#984EA3")) +
    labs(title = "Predicted Source Proportion by Group", x = "Group", y = "Proportion (%)", fill = "Predicted Source") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  p_true <- ggplot(true_proportions, aes(x = Group, y = True_Proportion, fill = True_Source)) +
    geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = c("A" = "#E41A1C", "B" = "#377EB8", "C" = "#4DAF4A", "A_B_Overlap" = "#984EA3")) +
    labs(title = "True Source Proportion by Group", x = "Group", y = "Proportion (%)", fill = "True Source") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  p_accuracy <- ggplot(accuracy_by_group, aes(x = Group, y = Group_Accuracy, fill = Group_Accuracy)) +
    geom_bar(stat = "identity") + geom_text(aes(label = sprintf("%.1f%%", Group_Accuracy)), vjust = -0.5, size = 4) +
    scale_fill_gradient(low = "#FF6B6B", high = "#4ECDC4", name = "Accuracy (%)") +
    labs(title = "Prediction Accuracy by Group", x = "Group", y = "Accuracy (%)") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(save_dir, "plot8_predicted_proportion.png"), p_pred, width = 10, height = 6, dpi = 300)
  ggsave(file.path(save_dir, "plot9_true_proportion.png"), p_true, width = 10, height = 6, dpi = 300)
  ggsave(file.path(save_dir, "plot10_group_accuracy.png"), p_accuracy, width = 10, height = 6, dpi = 300)
  
  return(list(predicted_proportions = predicted_proportions, true_proportions = true_proportions, proportion_comparison = proportion_comparison,
              accuracy_by_group = accuracy_by_group, plots = list(predicted_proportion_plot = p_pred, true_proportion_plot = p_true, accuracy_plot = p_accuracy)))
}

# Main execution function
multi_model_enhanced_ABC_version <- function(reference_data, test_data_all, validation_info) {
  sources <- c("A", "B", "C")
  
  tryCatch({
    preprocessed_data <- enhanced_ABC_feature_engineering(reference_data, test_data_all)
    fcm_results <- enhanced_ABC_fcm_algorithm(preprocessed_data$reference_data, preprocessed_data$test_data, sources, preprocessed_data$group_centers, preprocessed_data$source_features)
    final_results <- enhanced_ABC_post_processing(fcm_results, preprocessed_data$test_data, preprocessed_data$source_features, validation_info)
    performance_results <- performance_evaluation(final_results, validation_info, preprocessed_data$test_data)
    plots_enhanced <- create_two_enhanced_plots(performance_results$evaluation_data, performance_results, validation_dir, reference_data, test_data_all, validation_info)
    group_analysis_results <- group_proportion_analysis(performance_results$evaluation_data, validation_dir)
    
    final_output <- final_results[, c("Sample_ID", "Group", "X", "Y", "A", "B", "C", "Final_Prediction", "Max_Probability", "Confidence", "Sample_Type")]
    write.csv(final_output, file.path(validation_dir, "enhanced_ABC_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(performance_results$evaluation_data, file.path(validation_dir, "enhanced_ABC_evaluation.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    
    performance_summary <- data.frame(
      Metric = c("Overall_Accuracy", "A_Accuracy", "B_Accuracy", "C_Accuracy", "Overlap_Accuracy", "Total_Overlap_Samples", "Normal_Sample_Accuracy", "Total_Points", "Correct_Points"),
      Value = c(sprintf("%.2f%%", performance_results$accuracy_overall * 100), sprintf("%.2f%%", performance_results$class_accuracy["A"] * 100),
                sprintf("%.2f%%", performance_results$class_accuracy["B"] * 100), sprintf("%.2f%%", performance_results$class_accuracy["C"] * 100),
                sprintf("%.2f%%", performance_results$overlap_performance$overlap_accuracy * 100), performance_results$overlap_performance$total_overlap_samples,
                sprintf("%.2f%%", ifelse(!is.na(performance_results$sample_type_accuracy$Normal), performance_results$sample_type_accuracy$Normal * 100, 0)),
                performance_results$total_points, performance_results$correct_points)
    )
    write.csv(performance_summary, file.path(validation_dir, "enhanced_ABC_performance_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(performance_results$confusion_matrix, file.path(validation_dir, "enhanced_ABC_confusion_matrix.csv"), fileEncoding = "UTF-8")
    
    return(list(results = final_results, performance = performance_results, group_analysis = group_analysis_results, plots = plots_enhanced))
  }, error = function(e) NULL)
}

# Run
version12_enhanced_results <- multi_model_enhanced_ABC_version(reference_data, test_data_all, validation_info)