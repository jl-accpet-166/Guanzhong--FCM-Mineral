# Libraries
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(viridis)
library(ggrepel)
library(reshape2)
library(MASS)
library(mvtnorm)
library(FNN)
library(gridExtra)
library(cluster)
library(ggsci)

# Setup output directory
qinling_output_dir <- "./秦岭小区域溯源结果_增强版"
if (!dir.exists(qinling_output_dir)) dir.create(qinling_output_dir, recursive = TRUE)

# Load reference data
qinling_subregion_file <- "./data/qinling.xlsx"
qinling_subregion_data <- read_excel(qinling_subregion_file, sheet = 1)

region_map <- c("北秦岭" = "North Qinling", "东秦岭" = "East Qinling", 
                "南秦岭" = "South Qinling", "西秦岭" = "West Qinling", "小秦岭" = "Small Qinling")

qinling_clean <- data.frame(
  Sample_ID = if ("样品号" %in% names(qinling_subregion_data)) qinling_subregion_data$样品号 
              else if ("样品名称" %in% names(qinling_subregion_data)) qinling_subregion_data$样品名称 
              else paste0("Qinling_", 1:nrow(qinling_subregion_data)),
  Subregion = if ("地区" %in% names(qinling_subregion_data)) {
    sapply(qinling_subregion_data$地区, function(x) {
      if (x %in% names(region_map)) region_map[x]
      else { for (key in names(region_map)) if (grepl(key, x)) return(region_map[key]); "Unknown" }
    })
  } else rep("Unknown", nrow(qinling_subregion_data)),
  X = as.numeric(qinling_subregion_data$`207 Pb/206 Pb`),
  Y = as.numeric(qinling_subregion_data$`208Pb/206 Pb`),
  Mine_Name = if ("矿床名称" %in% names(qinling_subregion_data)) qinling_subregion_data$矿床名称 else NA
)
qinling_clean <- qinling_clean[!is.na(qinling_clean$X) & !is.na(qinling_clean$Y) & qinling_clean$Subregion != "Unknown", ]

# Load bronze samples
bronze_results_file <- "E:/李璟钰/博士/小论文/矿料/青铜器溯源FCM结果_平衡版/Bronze_Artefact_Ore_Source_Predictions_Balanced.csv"

if (file.exists(bronze_results_file)) {
  bronze_results <- read.csv(bronze_results_file, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
  colnames(bronze_results) <- gsub("^\\s+|\\s+$", "", gsub("\\.+", ".", gsub("\\s+", " ", trimws(colnames(bronze_results)))))
  
  final_pred_variants <- c("Final_Prediction", "Final Prediction", "Final.Prediction", "最终预测", "预测结果", "Predicted_Source", "预测区域", "Predicted Source", "Predicted.Region")
  final_pred_col <- final_pred_variants[final_pred_variants %in% colnames(bronze_results)][1]
  
  if (!is.na(final_pred_col)) {
    bronze_results$Final_Prediction <- as.character(bronze_results[[final_pred_col]])
    qinling_identifiers <- c("Qinling Region", "Qinling", "秦岭", "秦岭地区", "North Qinling", "South Qinling", 
                             "East Qinling", "West Qinling", "Small Qinling", "小秦岭", "Qinling region", "qinling region", "秦岭区域")
    qinling_samples <- bronze_results %>% filter(grepl(paste(qinling_identifiers, collapse = "|"), Final_Prediction, ignore.case = TRUE))
  }
  
  x_variants <- c("X", "207Pb/206Pb", "207Pb.206Pb", "Pb207_206", "207Pb/206 Pb")
  y_variants <- c("Y", "208Pb/206Pb", "208Pb.206Pb", "Pb208_206", "208Pb/206 Pb")
  id_variants <- c("Sample_ID", "Sample.ID", "SampleID", "样品号", "样品编号", "样品名称", "Sample Name", "Sample")
  period_variants <- c("period", "Period", "时期", "时代", "Dynasty", "Historical Period")
  
  for (v in x_variants) if (v %in% colnames(qinling_samples)) { qinling_samples$X <- qinling_samples[[v]]; break }
  for (v in y_variants) if (v %in% colnames(qinling_samples)) { qinling_samples$Y <- qinling_samples[[v]]; break }
  for (v in id_variants) if (v %in% colnames(qinling_samples)) { qinling_samples$Sample_ID <- qinling_samples[[v]]; break }
  for (v in period_variants) if (v %in% colnames(qinling_samples)) { qinling_samples$period <- qinling_samples[[v]]; break }
  
  if (!"X" %in% colnames(qinling_samples)) { nc <- which(sapply(qinling_samples, is.numeric)); if (length(nc) >= 1) qinling_samples$X <- qinling_samples[, nc[1]] }
  if (!"Y" %in% colnames(qinling_samples)) { nc <- which(sapply(qinling_samples, is.numeric)); if (length(nc) >= 2) qinling_samples$Y <- qinling_samples[, nc[2]] }
  if (!"Sample_ID" %in% colnames(qinling_samples)) qinling_samples$Sample_ID <- paste0("Bronze_", 1:nrow(qinling_samples))
  if (!"period" %in% colnames(qinling_samples)) qinling_samples$period <- "Unknown"
  if (!"Final_Prediction" %in% colnames(qinling_samples)) qinling_samples$Final_Prediction <- "Qinling Region"
  
  qinling_samples$X <- as.numeric(qinling_samples$X)
  qinling_samples$Y <- as.numeric(qinling_samples$Y)
  qinling_samples$Sample_ID <- as.character(qinling_samples$Sample_ID)
  qinling_samples$period <- as.character(qinling_samples$period)
  qinling_samples$Final_Prediction <- as.character(qinling_samples$Final_Prediction)
  qinling_samples <- qinling_samples[!is.na(qinling_samples$X) & !is.na(qinling_samples$Y), ]
  
  qinling_bronze_samples <- qinling_samples %>%
    mutate(Sample_ID = as.character(Sample_ID),
           period = ifelse(is.na(period) | period == "", "Unknown", period),
           period = factor(period, levels = c("Shang", "Western Zhou", "Spring and Autumn", "Warring States", "Unknown")),
           X = as.numeric(X), Y = as.numeric(Y))
} else {
  set.seed(123)
  n_samples <- 50
  qinling_bronze_samples <- data.frame(
    Sample_ID = paste0("Bronze_", 1:n_samples),
    period = factor(sample(c("Shang", "Western Zhou", "Spring and Autumn", "Warring States", "Unknown"), n_samples, replace = TRUE),
                    levels = c("Shang", "Western Zhou", "Spring and Autumn", "Warring States", "Unknown")),
    X = runif(n_samples, min(qinling_clean$X) - 0.01, max(qinling_clean$X) + 0.01),
    Y = runif(n_samples, min(qinling_clean$Y) - 0.01, max(qinling_clean$Y) + 0.01),
    Final_Prediction = "Qinling Region")
}

# Feature engineering function
enhanced_qinling_feature_engineering <- function(reference_data, test_data) {
  subregions <- unique(reference_data$Subregion)
  if (length(subregions) == 0) subregions <- c("North Qinling", "East Qinling", "South Qinling", "West Qinling", "Small Qinling")
  
  all_data <- rbind(reference_data[, c("X", "Y")] %>% mutate(Type = "Reference"), test_data[, c("X", "Y")] %>% mutate(Type = "Test"))
  median_vals <- apply(all_data[, c("X", "Y")], 2, median, na.rm = TRUE)
  iqr_vals <- apply(all_data[, c("X", "Y")], 2, IQR, na.rm = TRUE)
  iqr_vals[iqr_vals == 0] <- 1
  
  reference_data$X_std <- (reference_data$X - median_vals[1]) / iqr_vals[1]
  reference_data$Y_std <- (reference_data$Y - median_vals[2]) / iqr_vals[2]
  test_data$X_std <- (test_data$X - median_vals[1]) / iqr_vals[1]
  test_data$Y_std <- (test_data$Y - median_vals[2]) / iqr_vals[2]
  
  compute_distance_features <- function(data, ref_data, subs) {
    for (s in subs) { data[[paste0("Dist_to_", gsub(" ", "_", s))]] <- NA; data[[paste0("Dist_Ratio_", gsub(" ", "_", s))]] <- NA; data[[paste0("Angle_to_", gsub(" ", "_", s))]] <- NA }
    for (s in subs) {
      rd <- ref_data[ref_data$Subregion == s, ]
      if (nrow(rd) >= 3) {
        rc <- c(median(rd$X_std, na.rm = TRUE), median(rd$Y_std, na.rm = TRUE))
        data[[paste0("Dist_to_", gsub(" ", "_", s))]] <- sqrt((data$X_std - rc[1])^2 + (data$Y_std - rc[2])^2)
        data[[paste0("Angle_to_", gsub(" ", "_", s))]] <- atan2(data$Y_std - rc[2], data$X_std - rc[1])
      }
    }
    for (s in subs) {
      dc <- paste0("Dist_to_", gsub(" ", "_", s))
      if (dc %in% names(data)) {
        od <- sapply(setdiff(subs, s), function(os) { odc <- paste0("Dist_to_", gsub(" ", "_", os)); if (odc %in% names(data)) data[[odc]] else NA })
        mod <- mean(unlist(od), na.rm = TRUE)
        data[[paste0("Dist_Ratio_", gsub(" ", "_", s))]] <- data[[dc]] / (mod + 1e-6)
      }
    }
    data
  }
  
  reference_data <- compute_distance_features(reference_data, reference_data, subregions)
  test_data <- compute_distance_features(test_data, reference_data, subregions)
  
  boundaries <- list()
  for (s in subregions) {
    rd <- reference_data[reference_data$Subregion == s, ]
    if (nrow(rd) >= 3) {
      boundaries[[s]] <- list(x_min = min(rd$X_std, na.rm = TRUE), x_max = max(rd$X_std, na.rm = TRUE),
                              y_min = min(rd$Y_std, na.rm = TRUE), y_max = max(rd$Y_std, na.rm = TRUE),
                              center_x = mean(rd$X_std, na.rm = TRUE), center_y = mean(rd$Y_std, na.rm = TRUE))
    } else boundaries[[s]] <- list(x_min = NA, x_max = NA, y_min = NA, y_max = NA, center_x = NA, center_y = NA)
  }
  
  compute_protection <- function(data, bounds, subs) {
    for (s in subs) {
      b <- bounds[[s]]
      if (!is.na(b$x_min)) {
        xr <- b$x_max - b$x_min; yr <- b$y_max - b$y_min
        xp <- if (xr > 0) pmax(0, pmin(1, 1 - abs(2 * (data$X_std - b$center_x) / xr))) else 1
        yp <- if (yr > 0) pmax(0, pmin(1, 1 - abs(2 * (data$Y_std - b$center_y) / yr))) else 1
        prot <- sqrt(pmax(0, xp * yp))
        data[[paste0(gsub(" ", "_", s), "_Protection")]] <- prot
        dc <- paste0("Dist_to_", gsub(" ", "_", s))
        if (dc %in% names(data)) {
          md <- max(data[[dc]], na.rm = TRUE)
          nd <- if (md > 0) pmax(0, pmin(1, 1 - data[[dc]] / md)) else 1
          data[[paste0(gsub(" ", "_", s), "_Enhanced_Score")]] <- prot * nd * 100
        }
      } else { data[[paste0(gsub(" ", "_", s), "_Protection")]] <- 0; data[[paste0(gsub(" ", "_", s), "_Enhanced_Score")]] <- 0 }
    }
    data
  }
  
  reference_data <- compute_protection(reference_data, boundaries, subregions)
  test_data <- compute_protection(test_data, boundaries, subregions)
  
  groups <- unique(test_data$period)
  group_centers <- data.frame(Group = groups, Center_X = numeric(length(groups)), Center_Y = numeric(length(groups)))
  for (i in seq_along(groups)) {
    gd <- test_data[test_data$period == groups[i], ]
    if (nrow(gd) > 0) { group_centers$Center_X[i] <- median(gd$X_std, na.rm = TRUE); group_centers$Center_Y[i] <- median(gd$Y_std, na.rm = TRUE) }
  }
  
  source_features <- list()
  for (s in subregions) {
    sd <- reference_data[reference_data$Subregion == s, c("X_std", "Y_std")]
    n <- nrow(sd)
    if (n > 5) {
      sm <- colMeans(sd, na.rm = TRUE); sc <- cov(sd, use = "complete.obs")
      rc <- sm; rcov <- sc
      tryCatch({ mcd <- cov.mcd(sd, quantile.used = max(floor(0.75 * n), n)); rc <- mcd$center; rcov <- mcd$cov }, error = function(e) {})
      source_features[[s]] <- list(center = sm, robust_center = rc, cov = sc, robust_cov = rcov, n = n, boundaries = boundaries[[s]])
    } else if (n >= 2) {
      sm <- colMeans(sd, na.rm = TRUE); sc <- cov(sd, use = "complete.obs")
      source_features[[s]] <- list(center = sm, robust_center = sm, cov = sc, robust_cov = sc, n = n, boundaries = boundaries[[s]])
    } else if (n == 1) {
      source_features[[s]] <- list(center = as.numeric(sd), robust_center = as.numeric(sd), cov = diag(2) * 0.01, robust_cov = diag(2) * 0.01, n = 1, boundaries = boundaries[[s]])
    } else {
      source_features[[s]] <- list(center = c(0, 0), robust_center = c(0, 0), cov = diag(2) * 0.05, robust_cov = diag(2) * 0.05, n = 0, boundaries = list(x_min = NA, x_max = NA, y_min = NA, y_max = NA, center_x = 0, center_y = 0))
    }
  }
  
  list(reference_data = reference_data, test_data = test_data, subregions = subregions, group_centers = group_centers, source_features = source_features, boundaries = boundaries)
}

# FCM algorithm function
enhanced_qinling_fcm_algorithm <- function(reference_data, test_data, subregions, group_centers, source_features) {
  results <- data.frame(Sample_ID = test_data$Sample_ID, period = test_data$period, X = test_data$X, Y = test_data$Y, X_std = test_data$X_std, Y_std = test_data$Y_std)
  
  region_params <- list(
    "North Qinling" = list(mw = 0.60, ew = 0.40, rf = 0.30),
    "East Qinling" = list(mw = 0.58, ew = 0.42, rf = 0.35),
    "South Qinling" = list(mw = 0.55, ew = 0.45, rf = 0.40),
    "West Qinling" = list(mw = 0.62, ew = 0.38, rf = 0.28),
    "Small Qinling" = list(mw = 0.52, ew = 0.48, rf = 0.45))
  
  for (s in subregions) results[[s]] <- 50
  
  for (period in unique(test_data$period)) {
    ptd <- test_data[test_data$period == period, ]
    if (nrow(ptd) == 0) next
    pi <- which(test_data$period == period)
    tp <- as.matrix(ptd[, c("X_std", "Y_std")])
    gcr <- group_centers[group_centers$Group == period, ]
    gc <- if (nrow(gcr) > 0) c(gcr$Center_X, gcr$Center_Y) else c(0, 0)
    
    for (s in subregions) {
      params <- if (s %in% names(region_params)) region_params[[s]] else list(mw = 0.6, ew = 0.4, rf = 0.3)
      sc <- if (!is.null(source_features[[s]]$robust_center)) source_features[[s]]$robust_center else c(0, 0)
      gsd <- sqrt(sum((gc - sc)^2))
      sw <- exp(-gsd / 0.3)
      ac <- 0.8 * sc + 0.2 * gc * sw
      cc <- ac
      
      scr <- tryCatch({
        ev <- eigen(source_features[[s]]$robust_cov)$values
        me <- min(ev)
        if (me < 1e-6 || is.na(me)) diag(2) * (0.02 * params$rf + 0.01) else source_features[[s]]$robust_cov + diag(2) * 0.02 * params$rf
      }, error = function(e) diag(2) * 0.03)
      
      cp <- rep(50, nrow(tp))
      for (iter in 1:100) {
        md <- tryCatch({ m <- mahalanobis(tp, cc, scr); m[is.infinite(m) | m < 0 | is.na(m)] <- 1e-6; m }, error = function(e) sqrt(rowSums((tp - matrix(cc, nrow = nrow(tp), ncol = 2, byrow = TRUE))^2)))
        ed <- sqrt(rowSums((tp - matrix(cc, nrow = nrow(tp), ncol = 2, byrow = TRUE))^2))
        ed[is.infinite(ed) | ed < 0 | is.na(ed)] <- 1e-6
        cd <- params$mw * md + params$ew * ed
        if (all(is.na(cd))) cd <- rep(1, length(cd))
        sf <- max(cd, na.rm = TRUE) / 1.5
        if (sf <= 0 || is.na(sf) || is.infinite(sf)) sf <- 1
        np <- 100 * exp(-cd / sf)
        if (mean(abs(np - cp), na.rm = TRUE) < 1e-6 && iter > 10) break
        cp <- np
        w <- (cp/100)^1.5
        tw <- sum(w, na.rm = TRUE)
        if (tw > 0 && !is.na(tw)) {
          nc <- colSums(tp * w, na.rm = TRUE) / tw
          cdf <- nc - cc
          mvd <- sqrt(sum(cdf^2))
          if (mvd > 0.001) nc <- cc + cdf * (0.001 / mvd)
          cc <- 0.9 * cc + 0.1 * nc
        }
      }
      results[pi, s] <- cp
    }
  }
  
  sp <- results[, subregions, drop = FALSE]
  results$Predicted_Subregion <- apply(sp, 1, function(x) { if (all(is.na(x)) || all(x == 0)) NA else subregions[which.max(x)] })
  results$Max_Probability <- apply(sp, 1, function(x) { if (all(is.na(x)) || all(x == 0)) 50 else max(x, na.rm = TRUE) })
  results$Second_Probability <- apply(sp, 1, function(x) { if (all(is.na(x)) || all(x == 0)) 40 else { s <- sort(x, decreasing = TRUE); if (length(s) >= 2) s[2] else s[1] } })
  results$Probability_Ratio <- results$Max_Probability / (results$Second_Probability + 1e-6)
  results$Confidence_Level <- ifelse(results$Max_Probability >= 80 & results$Probability_Ratio >= 2.0, "High", ifelse(results$Max_Probability >= 65 & results$Probability_Ratio >= 1.6, "Medium", "Low"))
  results
}

# Post-processing function
enhanced_qinling_post_processing <- function(results, test_data, source_features) {
  subregions <- unique(test_data$Subregion)
  if (length(subregions) == 0) subregions <- c("North Qinling", "East Qinling", "South Qinling", "West Qinling", "Small Qinling")
  
  results$Final_Prediction <- results$Predicted_Subregion
  results$Optimized_Prediction <- results$Predicted_Subregion
  results$Is_Overlap <- FALSE
  results$Overlap_Type <- "None"
  results$Overlap_Confidence <- 0
  
  for (i in 1:(length(subregions)-1)) {
    for (j in (i+1):length(subregions)) {
      s1 <- subregions[i]; s2 <- subregions[j]
      pd <- abs(results[[s1]] - results[[s2]])
      oc <- pd < 15 & results$Max_Probability < 75 & results$Probability_Ratio < 1.8 & pmax(results[[s1]], results[[s2]]) > 55
      oconf <- ifelse(pd < 10 & results$Probability_Ratio < 1.5 & results$Max_Probability > 60, 0.85, ifelse(pd < 15 & results$Probability_Ratio < 1.8 & results$Max_Probability > 55, 0.70, 0.50))
      oi <- which(oc)
      if (length(oi) > 0) {
        results$Is_Overlap[oi] <- TRUE
        results$Overlap_Type[oi] <- ifelse(results$Overlap_Type[oi] == "None", paste(s1, s2, "Overlap", sep = "_"), paste(results$Overlap_Type[oi], paste(s1, s2, sep = "_"), sep = "+"))
        results$Overlap_Confidence[oi] <- pmax(results$Overlap_Confidence[oi], oconf[oi])
      }
    }
  }
  
  for (s in subregions) results[[paste0(s, "_Misclassified")]] <- FALSE
  results$Is_Borderline <- results$Probability_Ratio < 1.8 | results$Max_Probability < 70
  results$Sample_Type <- "Normal"
  results$Sample_Type[results$Is_Borderline] <- "Borderline"
  results$Sample_Type[results$Is_Overlap] <- "Overlap"
  results$Sample_Type[results$Max_Probability < 50] <- "Low_Confidence"
  
  for (s in subregions) {
    pc <- paste0(gsub(" ", "_", s), "_Protection"); scc <- paste0(gsub(" ", "_", s), "_Enhanced_Score")
    if (pc %in% names(test_data) && scc %in% names(test_data)) {
      esc <- paste0(gsub(" ", "_", setdiff(subregions, s)), "_Enhanced_Score")
      esc <- esc[esc %in% names(test_data)]
      if (length(esc) > 0) {
        mc <- results$Predicted_Subregion != s & !results$Is_Overlap & test_data[[pc]] > 0.75 & test_data[[scc]] > apply(test_data[, esc, drop = FALSE], 1, max, na.rm = TRUE) * 1.3
        tryCatch({ if (sum(mc, na.rm = TRUE) > 0) results[[paste0(s, "_Misclassified")]] <- mc }, error = function(e) {})
      }
    }
  }
  
  for (s in subregions) {
    mcc <- paste0(s, "_Misclassified")
    if (mcc %in% names(results)) {
      mi <- which(results[[mcc]])
      if (length(mi) > 0) {
        pc <- paste0(gsub(" ", "_", s), "_Protection"); scc <- paste0(gsub(" ", "_", s), "_Enhanced_Score"); drc <- paste0("Dist_Ratio_", gsub(" ", "_", s))
        for (idx in mi) {
          conf <- 0
          if (pc %in% names(test_data) && test_data[[pc]][idx] > 0.78) conf <- conf + 0.30
          if (scc %in% names(test_data)) {
            osc <- sapply(setdiff(subregions, s), function(os) { oscc <- paste0(gsub(" ", "_", os), "_Enhanced_Score"); if (oscc %in% names(test_data)) test_data[[oscc]][idx] else NA })
            mos <- max(osc, na.rm = TRUE)
            if (test_data[[scc]][idx] > mos * 1.3) conf <- conf + 0.25
          }
          if (drc %in% names(test_data) && !is.na(test_data[[drc]][idx]) && test_data[[drc]][idx] < 0.7) conf <- conf + 0.20
          if (conf >= 0.75) results$Optimized_Prediction[idx] <- s
        }
      }
    }
  }
  
  bi <- which(results$Is_Borderline & !results$Is_Overlap)
  for (idx in bi) {
    cp <- results$Optimized_Prediction[idx]
    for (s in subregions) {
      if (s == cp) next
      scc <- paste0(gsub(" ", "_", s), "_Enhanced_Score"); cscc <- paste0(gsub(" ", "_", cp), "_Enhanced_Score")
      if (scc %in% names(test_data) && cscc %in% names(test_data) && test_data[[scc]][idx] > test_data[[cscc]][idx] * 1.25) { results$Optimized_Prediction[idx] <- s; break }
    }
  }
  
  results$Final_Prediction <- results$Optimized_Prediction
  oi <- which(results$Is_Overlap)
  if (length(oi) > 0) for (idx in oi) if (results$Overlap_Confidence[idx] >= 0.85) results$Final_Prediction[idx] <- paste("Mixed", results$Overlap_Type[idx], sep = "_")
  results$Confidence_Level <- ifelse(results$Max_Probability >= 80 & results$Probability_Ratio >= 2.0, "High", ifelse(results$Max_Probability >= 65 & results$Probability_Ratio >= 1.6, "Medium", "Low"))
  results
}

# Visualization function
create_enhanced_qinling_visualizations <- function(results, reference_data, output_dir) {
  plot_list <- list()
  subregion_colors <- c("North Qinling" = "#3498db", "East Qinling" = "#e67e22", "South Qinling" = "#2ecc71", "West Qinling" = "#e74c3c", "Small Qinling" = "#9b59b6", "Mixed_North_East_Overlap" = "#1abc9c", "Mixed_North_South_Overlap" = "#f1c40f", "Mixed_East_South_Overlap" = "#e67e22")
  period_shapes <- c("Shang" = 16, "Western Zhou" = 18, "Spring and Autumn" = 15, "Warring States" = 17, "Unknown" = 4)
  
  p1 <- ggplot(results, aes(x = X, y = Y)) +
    geom_point(aes(color = Final_Prediction, shape = period), size = 4, alpha = 0.85, stroke = 1.2) +
    scale_color_manual(values = subregion_colors, name = "Predicted Subregion") +
    scale_shape_manual(values = period_shapes, name = "Historical Period") +
    labs(x = expression(paste(""^207, "Pb/", ""^206, "Pb")), y = expression(paste(""^208, "Pb/", ""^206, "Pb")), title = "Bronze Artefacts Subregion Prediction") +
    theme_classic() +
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.text.x = element_text(size = 14, face = "bold", color = "black"), axis.text.y = element_text(size = 14, face = "bold", color = "black"), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.position = "right", panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2), plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"))
  ggsave(file.path(output_dir, "Bronze_Artefacts_Scatter_Enhanced.png"), p1, width = 15, height = 9, dpi = 300, bg = "white")
  plot_list[["bronze_scatter_enhanced"]] <- p1
  
  if ("Confidence_Level" %in% names(results)) {
    p2 <- ggplot(results, aes(x = X, y = Y)) +
      geom_point(aes(color = Confidence_Level, shape = period), size = 5, alpha = 0.8) +
      scale_color_manual(values = c("High" = "#2ecc71", "Medium" = "#f39c12", "Low" = "#e74c3c"), name = "Confidence Level") +
      scale_shape_manual(values = period_shapes, name = "Historical Period") +
      labs(x = expression(paste(""^207, "Pb/", ""^206, "Pb")), y = expression(paste(""^208, "Pb/", ""^206, "Pb")), title = "Prediction Confidence Distribution") +
      theme_classic() +
      theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"), axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12, face = "bold"), panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2), plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt"))
    ggsave(file.path(output_dir, "Confidence_Distribution.png"), p2, width = 14, height = 8, dpi = 300, bg = "white")
    plot_list[["confidence_distribution"]] <- p2
  }
  plot_list
}

# Main execution function
enhanced_qinling_subregion_traceability <- function(reference_data, test_data) {
  tryCatch({
    preprocessed_data <- enhanced_qinling_feature_engineering(reference_data, test_data)
    fcm_results <- enhanced_qinling_fcm_algorithm(preprocessed_data$reference_data, preprocessed_data$test_data, preprocessed_data$subregions, preprocessed_data$group_centers, preprocessed_data$source_features)
    final_results <- enhanced_qinling_post_processing(fcm_results, preprocessed_data$test_data, preprocessed_data$source_features)
    plots <- create_enhanced_qinling_visualizations(final_results, reference_data, qinling_output_dir)
    
    write.csv(final_results, file.path(qinling_output_dir, "Qinling_Subregion_Traceability_Enhanced_Results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    
    subregion_cols <- c("North Qinling", "East Qinling", "South Qinling", "West Qinling", "Small Qinling")
    actual_cols <- subregion_cols[subregion_cols %in% colnames(final_results)]
    selected_cols <- c("Sample_ID", "period", "X", "Y", actual_cols, "Final_Prediction", "Max_Probability", "Sample_Type", "Is_Overlap")
    if ("Overlap_Type" %in% colnames(final_results)) selected_cols <- c(selected_cols, "Overlap_Type")
    selected_cols <- selected_cols[selected_cols %in% colnames(final_results)]
    write.csv(final_results[, selected_cols, drop = FALSE], file.path(qinling_output_dir, "Qinling_Subregion_Traceability_Results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    
    list(results = final_results, plots = plots, preprocessed_data = preprocessed_data)
  }, error = function(e) NULL)
}

# Execute analysis
test_data_for_enhanced <- data.frame(
  Sample_ID = as.character(qinling_bronze_samples$Sample_ID),
  period = as.character(qinling_bronze_samples$period),
  X = as.numeric(qinling_bronze_samples$X),
  Y = as.numeric(qinling_bronze_samples$Y))
test_data_for_enhanced$period <- factor(ifelse(test_data_for_enhanced$period %in% c("Shang", "Western Zhou", "Spring and Autumn", "Warring States", "Unknown"), test_data_for_enhanced$period, "Unknown"), levels = c("Shang", "Western Zhou", "Spring and Autumn", "Warring States", "Unknown"))
test_data_for_enhanced <- test_data_for_enhanced[!is.na(test_data_for_enhanced$X) & !is.na(test_data_for_enhanced$Y), ]

if (nrow(test_data_for_enhanced) > 0) enhanced_analysis_results <- enhanced_qinling_subregion_traceability(qinling_clean, test_data_for_enhanced)