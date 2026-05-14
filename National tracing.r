# Install required packages
if (!require("dplyr")) install.packages("dplyr")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("readxl")) install.packages("readxl")
if (!require("tidyr")) install.packages("tidyr")
if (!require("viridis")) install.packages("viridis")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("reshape2")) install.packages("reshape2")
if (!require("MASS")) install.packages("MASS")
if (!require("mvtnorm")) install.packages("mvtnorm")
if (!require("FNN")) install.packages("FNN")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cluster")) install.packages("cluster")
if (!require("colorspace")) install.packages("colorspace")

# Load libraries
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
library(colorspace)

# Create output directory
output_dir <- "./bronze_fcm_results_balanced"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Directory created: %s\n", output_dir))
}

# 1. Read data
cat("Reading data...\n")
bronze_data <- read_excel("./data/shang_zhou_bronze.xlsx", sheet = "关中地区文献中铅同位素")
mine_data <- read_excel("./data/National Mineral data.xlsx", sheet = "Sheet1")
qinling_data <- read_excel("./data/qinlingore1.xlsx", sheet = "Sheet1")

# Read Tongling, SouthChina and SouthwestChina data
data_group_dir <- "data_groups"
if (dir.exists(data_group_dir)) {
  cat(sprintf("Found data group directory: %s\n", data_group_dir))
  
  tongling_file <- file.path(data_group_dir, "Tongling.xlsx")
  southchina_file <- file.path(data_group_dir, "SouthChina.xlsx")
  southwest_file <- file.path(data_group_dir, "SouthwestChina.xlsx")
  
  if (file.exists(tongling_file)) {
    cat("Reading Tongling data...\n")
    tongling_data <- read_excel(tongling_file)
  } else {
    cat("Warning: Tongling.xlsx file does not exist\n")
    tongling_data <- NULL
  }
  
  if (file.exists(southchina_file)) {
    cat("Reading SouthChina data...\n")
    southchina_data <- read_excel(southchina_file)
  } else {
    cat("Warning: SouthChina.xlsx file does not exist\n")
    southchina_data <- NULL
  }
  
  if (file.exists(southwest_file)) {
    cat("Reading SouthwestChina data...\n")
    southwest_data <- read_excel(southwest_file, sheet = 1)
  } else {
    cat("Warning: SouthwestChina.xlsx file does not exist\n")
    southwest_data <- NULL
  }
} else {
  cat(sprintf("Warning: Data group directory does not exist: %s\n", data_group_dir))
  tongling_data <- NULL
  southchina_data <- NULL
  southwest_data <- NULL
}

# 2. Data preprocessing
# Qinling data processing
qinling_columns <- c("地区", "矿床名称", "样品号", "测试矿物", "207 Pb/206 Pb", "208Pb/206 Pb")
qinling_data_clean <- qinling_data[, qinling_columns]
colnames(qinling_data_clean) <- c("Region", "Mine_Name", "Sample_ID", "Mineral", "X", "Y")
qinling_data_clean$data_type <- "ore"
qinling_data_clean$Province <- "Shaanxi/Gansu"
qinling_data_clean$important_mine <- "Qinling Region"
qinling_data_clean$source <- "Qinling Region"
qinling_data_clean <- qinling_data_clean[!is.na(qinling_data_clean$X) & !is.na(qinling_data_clean$Y), ]

# Tongling, Southwest, South China data processing
mine_columns <- c("Sample No", "Site_Cn", "Admin_Province", "207Pb/206Pb", "208Pb/206Pb")
mine_data_clean <- mine_data[, mine_columns]
colnames(mine_data_clean) <- c("Sample_ID", "Mine_Name", "Province", "X", "Y")
mine_data_clean$data_type <- "ore"
mine_data_clean$important_mine <- ifelse(
  grepl("秦岭", mine_data_clean$Mine_Name), "Qinling Region",
  ifelse(
    grepl("铜陵", mine_data_clean$Mine_Name), "Tongling",
    ifelse(
      mine_data_clean$Province %in% c("广东", "福建", "Guangdong", "Fujian"), "South China",
      ifelse(
        mine_data_clean$Province %in% c("广西", "Guanxi"), "Southwest China",
        "Other Regions"
      )
    )
  )
)
mine_data_clean$source <- mine_data_clean$important_mine
mine_data_clean <- mine_data_clean[!is.na(mine_data_clean$X) & !is.na(mine_data_clean$Y), ]

# Process Tongling data
if (!is.null(tongling_data)) {
  tongling_cols <- names(tongling_data)
  find_column <- function(patterns, col_names) {
    for (pattern in patterns) {
      matched <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
      if (length(matched) > 0) return(matched[1])
    }
    return(NULL)
  }
  
  tongling_pb207_206 <- find_column(c("207Pb/206Pb", "207Pb.*206Pb", "Pb207_206", "X"), tongling_cols)
  tongling_pb208_206 <- find_column(c("208Pb/206Pb", "208Pb.*206Pb", "Pb208_206", "Y"), tongling_cols)
  tongling_sample_id <- find_column(c("Sample", "样品", "编号", "ID"), tongling_cols)
  
  if (!is.null(tongling_pb207_206) && !is.null(tongling_pb208_206)) {
    tongling_data_clean <- data.frame(
      Sample_ID = if (!is.null(tongling_sample_id)) tongling_data[[tongling_sample_id]] else paste0("Tongling_", 1:nrow(tongling_data)),
      X = as.numeric(tongling_data[[tongling_pb207_206]]),
      Y = as.numeric(tongling_data[[tongling_pb208_206]]),
      important_mine = "Tongling"
    )
    tongling_data_clean <- tongling_data_clean[!is.na(tongling_data_clean$X) & !is.na(tongling_data_clean$Y), ]
  } else {
    tongling_data_clean <- data.frame()
  }
} else {
  tongling_data_clean <- data.frame()
}

# Process SouthChina data
if (!is.null(southchina_data)) {
  southchina_cols <- names(southchina_data)
  find_column <- function(patterns, col_names) {
    for (pattern in patterns) {
      matched <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
      if (length(matched) > 0) return(matched[1])
    }
    return(NULL)
  }
  
  southchina_pb207_206 <- find_column(c("207Pb/206Pb", "207Pb.*206Pb", "Pb207_206", "X"), southchina_cols)
  southchina_pb208_206 <- find_column(c("208Pb/206Pb", "208Pb.*206Pb", "Pb208_206", "Y"), southchina_cols)
  southchina_sample_id <- find_column(c("Sample", "样品", "编号", "ID"), southchina_cols)
  
  if (!is.null(southchina_pb207_206) && !is.null(southchina_pb208_206)) {
    southchina_data_clean <- data.frame(
      Sample_ID = if (!is.null(southchina_sample_id)) southchina_data[[southchina_sample_id]] else paste0("SouthChina_", 1:nrow(southchina_data)),
      X = as.numeric(southchina_data[[southchina_pb207_206]]),
      Y = as.numeric(southchina_data[[southchina_pb208_206]]),
      important_mine = "South China"
    )
    southchina_data_clean <- southchina_data_clean[!is.na(southchina_data_clean$X) & !is.na(southchina_data_clean$Y), ]
  } else {
    southchina_data_clean <- data.frame()
  }
} else {
  southchina_data_clean <- data.frame()
}

# Process SouthwestChina data
if (!is.null(southwest_data)) {
  southwest_cols <- names(southwest_data)
  find_column <- function(patterns, col_names) {
    for (pattern in patterns) {
      matched <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
      if (length(matched) > 0) return(matched[1])
    }
    return(NULL)
  }
  
  pb207_206_col <- find_column(c("207Pb/206Pb", "207Pb.*206Pb", "Pb207.*Pb206", "207/206"), southwest_cols)
  pb208_206_col <- find_column(c("208Pb/206Pb", "208Pb.*206Pb", "Pb208.*Pb206", "208/206"), southwest_cols)
  pb206_204_col <- find_column(c("206Pb/204Pb", "206Pb.*204Pb", "206/204"), southwest_cols)
  pb207_204_col <- find_column(c("207Pb/204Pb", "207Pb.*204Pb", "207/204"), southwest_cols)
  pb208_204_col <- find_column(c("208Pb/204Pb", "208Pb.*204Pb", "208/204"), southwest_cols)
  southwest_sample_id <- find_column(c("Sample", "样品", "编号", "ID", "Sample No", "Sample.No", "ID"), southwest_cols)
  
  southwest_data_clean <- data.frame(
    Sample_ID = if (!is.null(southwest_sample_id)) {
      as.character(southwest_data[[southwest_sample_id]])
    } else {
      paste0("Southwest_", 1:nrow(southwest_data))
    },
    important_mine = "Southwest China"
  )
  
  southwest_data_clean$X <- NA
  southwest_data_clean$Y <- NA
  
  if (!is.null(pb207_206_col) && !is.null(pb208_206_col)) {
    valid_rows <- !is.na(southwest_data[[pb207_206_col]]) & !is.na(southwest_data[[pb208_206_col]])
    if (sum(valid_rows) > 0) {
      southwest_data_clean$X[valid_rows] <- as.numeric(southwest_data[[pb207_206_col]][valid_rows])
      southwest_data_clean$Y[valid_rows] <- as.numeric(southwest_data[[pb208_206_col]][valid_rows])
    }
  }
  
  if (!is.null(pb206_204_col) && !is.null(pb207_204_col) && !is.null(pb208_204_col)) {
    need_calc <- is.na(southwest_data_clean$X) | is.na(southwest_data_clean$Y)
    calc_rows <- need_calc & 
      !is.na(southwest_data[[pb206_204_col]]) & 
      !is.na(southwest_data[[pb207_204_col]]) & 
      !is.na(southwest_data[[pb208_204_col]])
    
    if (sum(calc_rows) > 0) {
      southwest_data_clean$X[calc_rows] <- as.numeric(southwest_data[[pb207_204_col]][calc_rows]) / 
        as.numeric(southwest_data[[pb206_204_col]][calc_rows])
      southwest_data_clean$Y[calc_rows] <- as.numeric(southwest_data[[pb208_204_col]][calc_rows]) / 
        as.numeric(southwest_data[[pb206_204_col]][calc_rows])
    }
  }
  
  southwest_data_clean <- southwest_data_clean[!is.na(southwest_data_clean$X) & !is.na(southwest_data_clean$Y), ]
} else {
  southwest_data_clean <- data.frame()
}

# Merge all mine data as reference data
cat("Merging reference data...\n")
mine_subset <- mine_data_clean[, c("Sample_ID", "X", "Y", "important_mine")]
qinling_subset <- qinling_data_clean[, c("Sample_ID", "X", "Y", "important_mine")]
reference_data <- rbind(mine_subset, qinling_subset)

if (nrow(tongling_data_clean) > 0) {
  reference_data <- rbind(reference_data, tongling_data_clean)
}
if (nrow(southchina_data_clean) > 0) {
  reference_data <- rbind(reference_data, southchina_data_clean)
}
if (nrow(southwest_data_clean) > 0) {
  reference_data <- rbind(reference_data, southwest_data_clean)
}

colnames(reference_data) <- c("Sample_ID", "X", "Y", "Known_Source")
reference_data <- reference_data[reference_data$Known_Source %in% 
                                   c("Qinling Region", "South China", "Tongling", "Southwest China"), ]
reference_data$Group <- "Reference"

# Process bronze data
bronze_columns <- c("Dynesty", "地区", "样品编号", "207Pb/206Pb", "208Pb/206Pb")
bronze_data_clean <- bronze_data[, bronze_columns]
colnames(bronze_data_clean) <- c("Dynasty", "Region", "Sample_ID", "X", "Y")
bronze_data_clean$data_type <- "bronze"
bronze_data_clean$period <- ifelse(
  grepl("商", bronze_data_clean$Dynasty), "Shang",
  ifelse(
    grepl("西周", bronze_data_clean$Dynasty), "Western Zhou",
    ifelse(
      grepl("春秋", bronze_data_clean$Dynasty), "Spring and Autumn",
      ifelse(
        grepl("战国", bronze_data_clean$Dynasty), "Warring States",
        "Other"
      )
    )
  )
)
bronze_data_clean$source <- "Bronze"
bronze_data_clean$Group <- bronze_data_clean$period
bronze_data_clean <- bronze_data_clean[!is.na(bronze_data_clean$X) & !is.na(bronze_data_clean$Y) & 
                                         bronze_data_clean$period != "Other", ]

# Output data statistics
cat(sprintf("Reference data (mines) total samples: %d\n", nrow(reference_data)))
cat(sprintf("Test data (bronzes) total samples: %d\n", nrow(bronze_data_clean)))
cat(sprintf("Shang period samples: %d\n", sum(bronze_data_clean$period == "Shang")))
cat(sprintf("Western Zhou period samples: %d\n", sum(bronze_data_clean$period == "Western Zhou")))
cat(sprintf("Spring and Autumn period samples: %d\n", sum(bronze_data_clean$period == "Spring and Autumn")))
cat(sprintf("Warring States period samples: %d\n", sum(bronze_data_clean$period == "Warring States")))

# 3. Balanced feature engineering function (moderate score adjustment, angle features retained)
balanced_feature_engineering <- function(reference_data, test_data_all) {
  cat("Starting balanced feature engineering...\n")
  
  # Robust standardisation
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
  
  # Compute group centres
  groups <- unique(test_data_all$Group)
  group_centers <- data.frame(
    Group = groups,
    Center_X = numeric(length(groups)),
    Center_Y = numeric(length(groups))
  )
  
  for (i in 1:length(groups)) {
    group_data <- test_data_all[test_data_all$Group == groups[i], ]
    group_centers$Center_X[i] <- median(group_data$X_std)
    group_centers$Center_Y[i] <- median(group_data$Y_std)
  }
  
  # MCD robust estimation for source feature extraction
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
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
        data_range = list(
          x = c(min(source_data$X_std) - boundary_extension, max(source_data$X_std) + boundary_extension),
          y = c(min(source_data$Y_std) - boundary_extension, max(source_data$Y_std) + boundary_extension)
        ), mcd_used = TRUE
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
        data_range = list(
          x = c(min(source_data$X_std) - boundary_extension, max(source_data$X_std) + boundary_extension),
          y = c(min(source_data$Y_std) - boundary_extension, max(source_data$Y_std) + boundary_extension)
        ), mcd_used = FALSE
      )
    } else {
      cov_factor <- 0.025
      boundary_extension <- 0.3
      
      source_features[[source]] <- list(
        center = as.numeric(source_data), robust_center = as.numeric(source_data),
        cov = diag(2) * cov_factor, robust_cov = diag(2) * cov_factor, n = 1,
        data_range = list(
          x = c(as.numeric(source_data[1]) - boundary_extension, as.numeric(source_data[1]) + boundary_extension),
          y = c(as.numeric(source_data[2]) - boundary_extension, as.numeric(source_data[2]) + boundary_extension)
        ), mcd_used = FALSE
      )
    }
  }
  
  # Compute basic distance features
  test_data_all <- test_data_all %>%
    mutate(
      Dist_to_Qinling = sqrt((X_std - source_features[["Qinling Region"]]$robust_center[1])^2 + 
                               (Y_std - source_features[["Qinling Region"]]$robust_center[2])^2),
      Dist_to_SouthChina = sqrt((X_std - source_features[["South China"]]$robust_center[1])^2 + 
                                  (Y_std - source_features[["South China"]]$robust_center[2])^2),
      Dist_to_Tongling = sqrt((X_std - source_features[["Tongling"]]$robust_center[1])^2 + 
                                (Y_std - source_features[["Tongling"]]$robust_center[2])^2),
      Dist_to_Southwest = sqrt((X_std - source_features[["Southwest China"]]$robust_center[1])^2 + 
                                 (Y_std - source_features[["Southwest China"]]$robust_center[2])^2),
      Dist_ratio_Qinling_SouthChina = Dist_to_Qinling / (Dist_to_SouthChina + 1e-6),
      Dist_ratio_Qinling_Tongling = Dist_to_Qinling / (Dist_to_Tongling + 1e-6),
      Dist_ratio_Qinling_Southwest = Dist_to_Qinling / (Dist_to_Southwest + 1e-6),
      Dist_ratio_SouthChina_Tongling = Dist_to_SouthChina / (Dist_to_Tongling + 1e-6),
      Dist_ratio_SouthChina_Southwest = Dist_to_SouthChina / (Dist_to_Southwest + 1e-6),
      Dist_ratio_Tongling_Southwest = Dist_to_Tongling / (Dist_to_Southwest + 1e-6)
    )
  
  # Angle features
  test_data_all <- test_data_all %>%
    mutate(
      Vector_Qinling_X = X_std - source_features[["Qinling Region"]]$robust_center[1],
      Vector_Qinling_Y = Y_std - source_features[["Qinling Region"]]$robust_center[2],
      Vector_SouthChina_X = X_std - source_features[["South China"]]$robust_center[1],
      Vector_SouthChina_Y = Y_std - source_features[["South China"]]$robust_center[2],
      Angle_Qinling_SouthChina = atan2(Vector_Qinling_Y, Vector_Qinling_X) - 
        atan2(Vector_SouthChina_Y, Vector_SouthChina_X),
      Angle_Qinling_SouthChina_norm = (Angle_Qinling_SouthChina + pi) / (2 * pi),
      Vector_Tongling_X = X_std - source_features[["Tongling"]]$robust_center[1],
      Vector_Tongling_Y = Y_std - source_features[["Tongling"]]$robust_center[2],
      Angle_Qinling_Tongling = atan2(Vector_Qinling_Y, Vector_Qinling_X) - 
        atan2(Vector_Tongling_Y, Vector_Tongling_X),
      Angle_Qinling_Tongling_norm = (Angle_Qinling_Tongling + pi) / (2 * pi),
      Angle_SouthChina_Tongling = atan2(Vector_SouthChina_Y, Vector_SouthChina_X) - 
        atan2(Vector_Tongling_Y, Vector_Tongling_X),
      Angle_SouthChina_Tongling_norm = (Angle_SouthChina_Tongling + pi) / (2 * pi)
    )
  
  # Balanced zone classification
  test_data_all <- test_data_all %>%
    mutate(
      Balanced_Zone = case_when(
        X_std < -0.10 & Y_std > 0.22 ~ "Qinling_Core",
        X_std < -0.03 & Y_std > 0.12 ~ "Qinling_Inner",
        X_std < 0.04 & Y_std > 0.05 ~ "Qinling_Outer",
        X_std > 0.25 & Y_std < -0.15 ~ "SouthChina_Core",
        X_std > 0.18 & Y_std < -0.08 ~ "SouthChina_Inner",
        X_std > 0.12 & Y_std < 0.02 ~ "SouthChina_Outer",
        (X_std > 0.36 & Y_std > 0.36) ~ "Tongling_Core",
        (X_std > 0.28 & Y_std > 0.28) ~ "Tongling_Inner",
        (X_std < -0.25 & Y_std < -0.25) ~ "Southwest_Core",
        (X_std < -0.18 & Y_std < -0.18) ~ "Southwest_Inner",
        Dist_ratio_Qinling_SouthChina >= 0.87 & Dist_ratio_Qinling_SouthChina <= 1.16 ~ "Qinling_SouthChina_Overlap",
        Dist_ratio_Qinling_Tongling >= 0.87 & Dist_ratio_Qinling_Tongling <= 1.16 ~ "Qinling_Tongling_Overlap",
        Dist_ratio_Qinling_Southwest >= 0.87 & Dist_ratio_Qinling_Southwest <= 1.16 ~ "Qinling_Southwest_Overlap",
        TRUE ~ "Neutral_Zone"
      )
    )
  
  # Balanced scoring system
  test_data_all <- test_data_all %>%
    mutate(
      Qinling_Score = (
        (1/(Dist_to_Qinling + 0.04)) * 0.45 +  
          exp(-Dist_ratio_Qinling_SouthChina) * 0.18 +
          (1 - abs(Angle_Qinling_SouthChina_norm - 0.5)) * 0.08 +  
          case_when(
            Balanced_Zone == "Qinling_Core" ~ 0.29,
            Balanced_Zone == "Qinling_Inner" ~ 0.22,
            Balanced_Zone == "Qinling_Outer" ~ 0.12,
            TRUE ~ 0
          )
      ),
      SouthChina_Score = (
        (1/(Dist_to_SouthChina + 0.04)) * 0.45 +
          exp(Dist_ratio_Qinling_SouthChina - 1) * 0.18 +
          abs(Angle_Qinling_SouthChina_norm - 0.5) * 0.08 +
          case_when(
            Balanced_Zone == "SouthChina_Core" ~ 0.29,
            Balanced_Zone == "SouthChina_Inner" ~ 0.22,
            Balanced_Zone == "SouthChina_Outer" ~ 0.12,
            TRUE ~ 0
          )
      ),
      Tongling_Score = (
        (1/(Dist_to_Tongling + 0.05)) * 0.48 +  
          exp(-pmin(Dist_ratio_Qinling_Tongling, Dist_ratio_Qinling_Southwest)) * 0.18 +
          (1 - abs(Angle_Qinling_Tongling_norm - 0.5)) * 0.05 +  
          case_when(
            Balanced_Zone == "Tongling_Core" ~ 0.24,
            Balanced_Zone == "Tongling_Inner" ~ 0.17,
            TRUE ~ 0
          )
      ),
      Southwest_Score = (
        (1/(Dist_to_Southwest + 0.08)) * 0.35 +  
          exp(-pmin(Dist_ratio_Qinling_Southwest, Dist_ratio_Qinling_Tongling)) * 0.12 +  
          (1 - abs(Angle_SouthChina_Tongling_norm - 0.5)) * 0.03 +  
          case_when(
            Balanced_Zone == "Southwest_Core" ~ 0.18,  
            Balanced_Zone == "Southwest_Inner" ~ 0.12,  
            TRUE ~ 0
          )
      )
    )
  
  # Determine dominant source
  test_data_all <- test_data_all %>%
    mutate(
      Dominant_Source = case_when(
        Qinling_Score > SouthChina_Score & Qinling_Score > Tongling_Score & Qinling_Score > Southwest_Score ~ "Qinling Region",
        SouthChina_Score > Qinling_Score & SouthChina_Score > Tongling_Score & SouthChina_Score > Southwest_Score ~ "South China",
        Tongling_Score > Qinling_Score & Tongling_Score > SouthChina_Score & Tongling_Score > Southwest_Score ~ "Tongling",
        Southwest_Score > Qinling_Score & Southwest_Score > SouthChina_Score & Southwest_Score > Tongling_Score ~ "Southwest China",
        TRUE ~ "Tie"
      )
    )
  
  cat("Balanced feature engineering completed!\n")
  
  return(list(
    reference_data = reference_data, test_data = test_data_all, group_centers = group_centers,
    source_features = source_features, norm_params = list(median = median_vals, iqr = iqr_vals)
  ))
}

# 4. Balanced FCM algorithm
balanced_fcm_algorithm <- function(reference_data, test_data, sources, group_centers, source_features) {
  cat("Executing balanced FCM algorithm...\n")
  
  results <- data.frame(
    Sample_ID = test_data$Sample_ID, Group = test_data$Group, period = test_data$period,
    X = test_data$X, Y = test_data$Y, X_std = test_data$X_std, Y_std = test_data$Y_std
  )
  
  # FCM parameters
  max_iterations <- 120
  tolerance <- 1e-7
  fuzziness <- 1.4
  movement_limit <- 0.0012
  base_regularization <- 0.016
  
  # Balanced source-specific parameters
  source_specific_params <- list(
    "Qinling Region" = list(mahalanobis_weight = 0.58, euclidean_weight = 0.42, angle_weight = 0.08, regularization_factor = 0.35),
    "South China" = list(mahalanobis_weight = 0.52, euclidean_weight = 0.48, angle_weight = 0.08, regularization_factor = 0.35),
    "Tongling" = list(mahalanobis_weight = 0.60, euclidean_weight = 0.40, angle_weight = 0.04, regularization_factor = 0.30),
    "Southwest China" = list(mahalanobis_weight = 0.50, euclidean_weight = 0.50, angle_weight = 0.03, regularization_factor = 0.25)
  )
  
  # Create probability columns for each source
  for (source in sources) {
    results[[source]] <- 0
  }
  
  # Process data by period group
  groups <- unique(test_data$Group)
  
  for (group in groups) {
    group_test_data <- test_data[test_data$Group == group, ]
    if (nrow(group_test_data) == 0) next
    
    group_indices <- which(test_data$Group == group)
    test_points <- as.matrix(group_test_data[, c("X_std", "Y_std")])
    
    # Compute probability for each source
    for (source in sources) {
      params <- source_specific_params[[source]]
      
      # Source centre strategy
      source_center <- source_features[[source]]$robust_center
      
      # Group centre adjustment
      group_center_row <- group_centers[group_centers$Group == group, ]
      group_center <- c(group_center_row$Center_X, group_center_row$Center_Y)
      group_source_dist <- sqrt(sum((group_center - source_center)^2))
      similarity_weight <- exp(-group_source_dist / 0.30)
      
      adjusted_center <- 0.80 * source_center + (1 - 0.80) * group_center * similarity_weight
      current_center <- adjusted_center
      
      # Enhanced covariance processing
      tryCatch({
        eigen_vals <- eigen(source_features[[source]]$robust_cov)$values
        min_eigen <- min(eigen_vals)
        reg_factor <- base_regularization * params$regularization_factor
        if (min_eigen < 1e-6) {
          source_cov_reg <- source_features[[source]]$robust_cov + diag(2) * (reg_factor + 0.010)
        } else {
          source_cov_reg <- source_features[[source]]$robust_cov + diag(2) * reg_factor
        }
      }, error = function(e) {
        source_cov_reg <- diag(2) * 0.032
      })
      
      current_probs <- rep(0.45, nrow(test_points))
      
      # FCM iteration
      for (iter in 1:max_iterations) {
        # Distance computation
        tryCatch({
          mahalanobis_dist <- mahalanobis(test_points, current_center, source_cov_reg)
          mahalanobis_dist[is.infinite(mahalanobis_dist) | mahalanobis_dist < 0 | is.na(mahalanobis_dist)] <- 1e-6
        }, error = function(e) {
          mahalanobis_dist <- sqrt(rowSums((test_points - matrix(current_center, nrow = nrow(test_points), ncol = 2, byrow = TRUE))^2))
        })
        
        euclidean_dist <- sqrt(rowSums((test_points - matrix(current_center, nrow = nrow(test_points), ncol = 2, byrow = TRUE))^2))
        euclidean_dist[is.infinite(euclidean_dist) | euclidean_dist < 0 | is.na(euclidean_dist)] <- 1e-6
        
        # Angle distance feature (balanced weight)
        if (source == "Qinling Region") {
          angle_dist <- abs(test_data$Angle_Qinling_SouthChina_norm[group_indices] - 0.5) * 2
        } else if (source == "South China") {
          angle_dist <- (1 - abs(test_data$Angle_Qinling_SouthChina_norm[group_indices] - 0.5)) * 2
        } else if (source == "Tongling") {
          angle_dist <- abs(test_data$Angle_Qinling_Tongling_norm[group_indices] - 0.5) * 2
        } else {
          angle_dist <- abs(test_data$Angle_SouthChina_Tongling_norm[group_indices] - 0.5) * 2
        }
        
        # Balanced distance combination
        combined_dist <- params$mahalanobis_weight * mahalanobis_dist + 
          params$euclidean_weight * euclidean_dist +
          params$angle_weight * angle_dist * 8
        
        if (all(is.na(combined_dist))) {
          combined_dist <- rep(1, length(combined_dist))
        }
        
        # Compute new probabilities with numerical stability
        scale_factor <- max(combined_dist, na.rm = TRUE) / 1.65
        if (scale_factor <= 0 || !is.finite(scale_factor)) scale_factor <- 1
        combined_dist[!is.finite(combined_dist)] <- 1
        new_probs <- 100 * exp(-combined_dist / scale_factor)
        new_probs[!is.finite(new_probs)] <- 0
        new_probs[new_probs < 0] <- 0
        new_probs[new_probs > 100] <- 100
        
        # Convergence check
        prob_change <- mean(abs(new_probs - current_probs))
        if (prob_change < tolerance && iter > 20) {
          break
        }
        
        current_probs <- new_probs
        
        # Update centre with safe weight computation
        prob_scaled <- current_probs/100
        prob_scaled[prob_scaled < 0] <- 0
        prob_scaled[prob_scaled > 1] <- 1
        weights <- prob_scaled^fuzziness
        weights[!is.finite(weights)] <- 0.001
        total_weight <- sum(weights)
        
        if (total_weight > 0 && !is.na(total_weight)) {
          new_center <- colSums(test_points * weights) / total_weight
          center_diff <- new_center - current_center
          move_dist <- sqrt(sum(center_diff^2))
          if (!is.na(move_dist) && !is.infinite(move_dist) && move_dist > movement_limit) {
            new_center <- current_center + center_diff * (movement_limit / move_dist)
          }
          if (!any(is.na(new_center)) && !any(is.infinite(new_center))) {
            current_center <- 0.85 * current_center + 0.15 * new_center
          }
        }
      }
      
      # Save results
      results[group_indices, source] <- current_probs
    }
  }
  
  # Basic prediction and metrics
  source_probs <- results[, sources]
  results$Predicted_Source <- apply(source_probs, 1, function(x) {
    if (all(is.na(x)) || all(x == 0)) return(NA)
    sources[which.max(x)]
  })
  results$Max_Probability <- apply(source_probs, 1, function(x) {
    if (all(is.na(x)) || all(x == 0)) return(45)
    max(x, na.rm = TRUE)
  })
  
  # Probability differences
  results$Prob_Diff_Qinling_SouthChina <- abs(results$`Qinling Region` - results$`South China`)
  results$Prob_Diff_Qinling_Tongling <- abs(results$`Qinling Region` - results$Tongling)
  results$Prob_Diff_Qinling_Southwest <- abs(results$`Qinling Region` - results$`Southwest China`)
  results$Prob_Diff_SouthChina_Tongling <- abs(results$`South China` - results$Tongling)
  results$Prob_Diff_SouthChina_Southwest <- abs(results$`South China` - results$`Southwest China`)
  results$Prob_Diff_Tongling_Southwest <- abs(results$Tongling - results$`Southwest China`)
  
  # Second probability
  results$Second_Probability <- apply(source_probs, 1, function(x) {
    if (all(is.na(x)) || all(x == 0)) return(40)
    sorted_probs <- sort(x, decreasing = TRUE)
    if (length(sorted_probs) >= 2) sorted_probs[2] else sorted_probs[1]
  })
  
  results$Probability_Ratio <- results$Max_Probability / (results$Second_Probability + 1e-6)
  
  cat("Balanced FCM algorithm completed!\n")
  
  return(results)
}

# 5. Balanced post-processing function
balanced_post_processing <- function(results, test_data) {
  cat("Executing balanced post-processing...\n")
  
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
  results$Optimized_Prediction <- results$Predicted_Source
  
  # Compute distance ratios
  test_data$Dist_ratio_Qinling_SouthChina <- ifelse(is.finite(test_data$Dist_to_Qinling) & is.finite(test_data$Dist_to_SouthChina), test_data$Dist_to_Qinling / (test_data$Dist_to_SouthChina + 1e-6), 1)
  test_data$Dist_ratio_Qinling_Tongling <- ifelse(is.finite(test_data$Dist_to_Qinling) & is.finite(test_data$Dist_to_Tongling), test_data$Dist_to_Qinling / (test_data$Dist_to_Tongling + 1e-6), 1)
  test_data$Dist_ratio_Qinling_Southwest <- ifelse(is.finite(test_data$Dist_to_Qinling) & is.finite(test_data$Dist_to_Southwest), test_data$Dist_to_Qinling / (test_data$Dist_to_Southwest + 1e-6), 1)
  test_data$Dist_ratio_SouthChina_Tongling <- ifelse(is.finite(test_data$Dist_to_SouthChina) & is.finite(test_data$Dist_to_Tongling), test_data$Dist_to_SouthChina / (test_data$Dist_to_Tongling + 1e-6), 1)
  test_data$Dist_ratio_SouthChina_Southwest <- ifelse(is.finite(test_data$Dist_to_SouthChina) & is.finite(test_data$Dist_to_Southwest), test_data$Dist_to_SouthChina / (test_data$Dist_to_Southwest + 1e-6), 1)
  test_data$Dist_ratio_Tongling_Southwest <- ifelse(is.finite(test_data$Dist_to_Tongling) & is.finite(test_data$Dist_to_Southwest), test_data$Dist_to_Tongling / (test_data$Dist_to_Southwest + 1e-6), 1)
  
  # Identify overlap samples
  overlap_categories <- list()
  
  # 1. Qinling-South China overlap
  qinling_southchina_overlap <- which(
    test_data$Dist_ratio_Qinling_SouthChina >= 0.87 & test_data$Dist_ratio_Qinling_SouthChina <= 1.16 &
      results$Prob_Diff_Qinling_SouthChina < 22 & results$Max_Probability < 78)
  if (length(qinling_southchina_overlap) > 0) {
    results$Optimized_Prediction[qinling_southchina_overlap] <- "Qinling_SouthChina_Overlap"
    overlap_categories[["Qinling-SouthChina Overlap"]] <- length(qinling_southchina_overlap)
  }
  
  # 2. Qinling-Tongling overlap
  qinling_tongling_overlap <- which(
    test_data$Dist_ratio_Qinling_Tongling >= 0.87 & test_data$Dist_ratio_Qinling_Tongling <= 1.16 &
      results$Prob_Diff_Qinling_Tongling < 22 & results$Max_Probability < 78)
  if (length(qinling_tongling_overlap) > 0) {
    results$Optimized_Prediction[qinling_tongling_overlap] <- "Qinling_Tongling_Overlap"
    overlap_categories[["Qinling-Tongling Overlap"]] <- length(qinling_tongling_overlap)
  }
  
  # 3. Qinling-Southwest overlap
  qinling_southwest_overlap <- which(
    test_data$Dist_ratio_Qinling_Southwest >= 0.87 & test_data$Dist_ratio_Qinling_Southwest <= 1.16 &
      results$Prob_Diff_Qinling_Southwest < 22 & results$Max_Probability < 78)
  if (length(qinling_southwest_overlap) > 0) {
    results$Optimized_Prediction[qinling_southwest_overlap] <- "Qinling_Southwest_Overlap"
    overlap_categories[["Qinling-Southwest Overlap"]] <- length(qinling_southwest_overlap)
  }
  
  # 4. South China-Tongling overlap
  southchina_tongling_overlap <- which(
    test_data$Dist_ratio_SouthChina_Tongling >= 0.87 & test_data$Dist_ratio_SouthChina_Tongling <= 1.16 &
      results$Prob_Diff_SouthChina_Tongling < 22 & results$Max_Probability < 78)
  if (length(southchina_tongling_overlap) > 0) {
    results$Optimized_Prediction[southchina_tongling_overlap] <- "SouthChina_Tongling_Overlap"
    overlap_categories[["SouthChina-Tongling Overlap"]] <- length(southchina_tongling_overlap)
  }
  
  # 5. South China-Southwest overlap
  southchina_southwest_overlap <- which(
    test_data$Dist_ratio_SouthChina_Southwest >= 0.87 & test_data$Dist_ratio_SouthChina_Southwest <= 1.16 &
      results$Prob_Diff_SouthChina_Southwest < 22 & results$Max_Probability < 78)
  if (length(southchina_southwest_overlap) > 0) {
    results$Optimized_Prediction[southchina_southwest_overlap] <- "SouthChina_Southwest_Overlap"
    overlap_categories[["SouthChina-Southwest Overlap"]] <- length(southchina_southwest_overlap)
  }
  
  # 6. Tongling-Southwest overlap
  tongling_southwest_overlap <- which(
    test_data$Dist_ratio_Tongling_Southwest >= 0.87 & test_data$Dist_ratio_Tongling_Southwest <= 1.16 &
      results$Prob_Diff_Tongling_Southwest < 22 & results$Max_Probability < 78)
  if (length(tongling_southwest_overlap) > 0) {
    results$Optimized_Prediction[tongling_southwest_overlap] <- "Tongling_Southwest_Overlap"
    overlap_categories[["Tongling-Southwest Overlap"]] <- length(tongling_southwest_overlap)
  }
  
  results$Final_Prediction <- results$Optimized_Prediction
  
  # Sample type identification
  results$Sample_Type <- "Normal"
  overlap_indices <- unique(c(
    qinling_southchina_overlap, qinling_tongling_overlap, qinling_southwest_overlap,
    southchina_tongling_overlap, southchina_southwest_overlap, tongling_southwest_overlap
  ))
  results$Sample_Type[overlap_indices] <- "Overlap_Region"
  
  # Output overlap statistics
  cat("\nOverlap sample statistics:\n")
  for (category in names(overlap_categories)) {
    cat(sprintf("  %s: %d samples\n", category, overlap_categories[[category]]))
  }
  cat(sprintf("Total overlap samples: %d\n", length(overlap_indices)))
  
  return(results)
}

# 6. Enhanced visualisation function
enhanced_visualizations <- function(results, output_dir) {
  cat("Generating enhanced visualisations...\n")
  
  plot_data <- results
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
  
  # Extract most probable source for each sample
  plot_data$Max_Prob_Source <- apply(plot_data[, sources], 1, function(row) {
    row[!is.finite(row)] <- 0
    if (all(row == 0)) return(NA)
    max_prob_idx <- which.max(row)
    sources[max_prob_idx]
  })
  
  # Mark overlap samples
  plot_data$Is_Overlap <- grepl("Overlap", plot_data$Final_Prediction)
  
  # Style configuration
  source_colors <- c(
    "Qinling Region" = "#5B2C6F", "South China" = "#2980B9",
    "Tongling" = "#F1C40F", "Southwest China" = "#27AE60"
  )
  
  period_shapes <- c(
    "Shang" = 16, "Spring and Autumn" = 15,
    "Warring States" = 17, "Western Zhou" = 18
  )
  
  # Plot 1: Main scatter plot (all samples, coloured by most probable source)
  p1 <- ggplot(plot_data, aes(x = X, y = Y, colour = Max_Prob_Source, shape = period)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_colour_manual(values = source_colors, name = "Most Probable Ore Source", limits = names(source_colors)) +
    scale_shape_manual(values = period_shapes, name = "Period", limits = names(period_shapes)) +
    labs(x = "207Pb/206Pb", y = "208Pb/206Pb",
         title = "Lead Isotope Ratios of Bronze Artefacts (Coloured by Most Probable Ore Source)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
          legend.position = "right")
  
  ggsave(file.path(output_dir, "Lead_Isotope_Distribution_Bronze.png"), p1, width = 14, height = 8, dpi = 300)
  cat("Generated: Lead_Isotope_Distribution_Bronze.png\n")
  
  # Plot 2: Overlap samples scatter plot
  if (sum(plot_data$Is_Overlap) > 0) {
    p2 <- ggplot(plot_data, aes(x = X, y = Y, colour = ifelse(Is_Overlap, "Overlap", "Normal"), shape = period)) +
      geom_point(size = 4) +
      scale_colour_manual(values = c("Overlap" = "#E74C3C", "Normal" = "#BDC3C7"),
                          name = "Sample Type", labels = c("Overlap Samples", "Normal Samples")) +
      scale_shape_manual(values = period_shapes, name = "Period", limits = names(period_shapes)) +
      labs(x = "207Pb/206Pb", y = "208Pb/206Pb", title = "Distribution of Overlap Samples") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12),
            legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11)) +
      annotate("text", x = min(plot_data$X), y = max(plot_data$Y),
               label = sprintf("Overlap Samples: %d", sum(plot_data$Is_Overlap)),
               hjust = 0, vjust = 1, size = 5, fontface = "bold", colour = "#E74C3C")
    
    ggsave(file.path(output_dir, "Overlap_Samples_Distribution.png"), p2, width = 14, height = 8, dpi = 300)
    cat("Generated: Overlap_Samples_Distribution.png\n")
  }
  
  # Plot 3: Source probability bar chart by period
  bar_data_all <- data.frame(
    period = plot_data$period,
    `Qinling Region` = plot_data$`Qinling Region`,
    `South China` = plot_data$`South China`,
    `Tongling` = plot_data$Tongling,
    `Southwest China` = plot_data$`Southwest China`
  ) %>%
    reshape2::melt(id.vars = "period", variable.name = "Mine", value.name = "Probability")
  
  bar_summary_all <- aggregate(Probability ~ period + Mine, bar_data_all, mean) %>%
    group_by(period) %>%
    arrange(desc(Probability), .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Mine = factor(Mine, levels = sources))
  
  p3 <- ggplot(bar_summary_all, aes(x = Mine, y = Probability, fill = Mine)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = source_colors) +
    facet_wrap(~ period, nrow = 2, ncol = 2) +
    labs(x = "Ore Source", y = "Average Probability (%)",
         title = "Average Ore Source Probabilities by Period") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          legend.position = "none",
          strip.text = element_text(size = 13, face = "bold"),
          panel.spacing = unit(1, "cm")) +
    geom_text(aes(label = sprintf("%.1f", Probability)), vjust = -0.3, size = 4, fontface = "bold") +
    scale_y_continuous(limits = c(0, 100))
  
  ggsave(file.path(output_dir, "Average_Ore_Source_Probabilities_by_Period.png"), p3, width = 14, height = 10, dpi = 300)
  cat("Generated: Average_Ore_Source_Probabilities_by_Period.png\n")
  
  if (sum(plot_data$Is_Overlap) > 0) {
    # Plot 4a: Overlap sample count by period
    overlap_count <- plot_data %>%
      filter(Is_Overlap) %>%
      group_by(period) %>%
      summarise(Count = n(), .groups = "drop")
    
    p4a <- ggplot(overlap_count, aes(x = period, y = Count, fill = period)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "Period", y = "Number of Overlap Samples",
           title = "Number of Overlap Samples by Period") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12),
            legend.position = "none") +
      geom_text(aes(label = Count), vjust = -0.3, size = 5, fontface = "bold")
    
    ggsave(file.path(output_dir, "Overlap_Sample_Count_by_Period.png"), p4a, width = 10, height = 8, dpi = 300)
    cat("Generated: Overlap_Sample_Count_by_Period.png\n")
    
    # Plot 4b: Overlap sample type distribution
    overlap_types <- plot_data %>%
      filter(Is_Overlap) %>%
      mutate(Overlap_Type = case_when(
        grepl("Qinling_SouthChina", Final_Prediction) ~ "Qinling-South China",
        grepl("Qinling_Tongling", Final_Prediction) ~ "Qinling-Tongling",
        grepl("Qinling_Southwest", Final_Prediction) ~ "Qinling-Southwest",
        grepl("SouthChina_Tongling", Final_Prediction) ~ "South China-Tongling",
        grepl("SouthChina_Southwest", Final_Prediction) ~ "South China-Southwest",
        grepl("Tongling_Southwest", Final_Prediction) ~ "Tongling-Southwest",
        TRUE ~ "Other"
      )) %>%
      group_by(Overlap_Type) %>%
      summarise(Count = n(), .groups = "drop") %>%
      arrange(desc(Count)) %>%
      mutate(Overlap_Type = factor(Overlap_Type, levels = Overlap_Type))
    
    p4b <- ggplot(overlap_types, aes(x = Overlap_Type, y = Count, fill = Overlap_Type)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_brewer(palette = "Set3") +
      labs(x = "Overlap Type", y = "Number of Samples",
           title = "Distribution of Overlap Sample Types") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            legend.position = "none") +
      geom_text(aes(label = Count), vjust = -0.3, size = 4, fontface = "bold") +
      ylim(0, max(overlap_types$Count) * 1.2)
    
    ggsave(file.path(output_dir, "Overlap_Sample_Types_Distribution.png"), p4b, width = 10, height = 8, dpi = 300)
    cat("Generated: Overlap_Sample_Types_Distribution.png\n")
  }
  
  return(list(
    main_scatter = p1,
    overlap_scatter = if (sum(plot_data$Is_Overlap) > 0) p2 else NULL,
    probability_bar = p3,
    overlap_count_chart = if (sum(plot_data$Is_Overlap) > 0) p4a else NULL,
    overlap_types_chart = if (sum(plot_data$Is_Overlap) > 0) p4b else NULL
  ))
}

# 7. Group analysis function
balanced_group_analysis <- function(results, output_dir) {
  cat("Starting group analysis...\n")
  
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
  
  # Replace all predictions by max probability source (overlap samples merged)
  results$Final_Prediction <- apply(results[, sources], 1, function(row) {
    sources[which.max(row)]
  })
  
  # Compute proportions for the four base sources
  period_source_prop <- results %>%
    group_by(period, Final_Prediction) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(period) %>%
    mutate(Total = sum(Count), Proportion = Count / Total * 100) %>%
    ungroup() %>%
    mutate(
      period = factor(period, levels = c("Shang", "Western Zhou", "Spring and Autumn", "Warring States")),
      Final_Prediction = factor(Final_Prediction, levels = sources)
    ) %>%
    arrange(period, desc(Proportion))
  
  # Save merged data
  write.csv(period_source_prop, file.path(output_dir, "Ore_Source_Proportion_by_Period.csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
  
  # Plot
  source_colors <- c(
    "Qinling Region" = "#5B2C6F", "South China" = "#2980B9",
    "Tongling" = "#F1C40F", "Southwest China" = "#27AE60"
  )
  
  p <- ggplot(period_source_prop, aes(x = period, y = Proportion, fill = Final_Prediction)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = source_colors, name = "Predicted Ore Source") +
    labs(title = "Predicted Ore Source Distribution by Period (Overlap Samples Merged to Most Probable Source)",
         x = "Historical Period", y = "Proportion (%)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11)) +
    geom_text(aes(label = sprintf("%.1f%%", Proportion)),
              position = position_stack(vjust = 0.5), size = 4, fontface = "bold")
  
  ggsave(file.path(output_dir, "Predicted_Ore_Source_Distribution_by_Period.png"), p, width = 12, height = 8, dpi = 300)
  
  return(list(data = period_source_prop, plot = p))
}

# 8. Main execution function
balanced_bronze_analysis <- function(reference_data, test_data_all) {
  cat("Starting balanced bronze artefact lead isotope provenance analysis...\n")
  
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
  
  tryCatch({
    # Step 1: Feature engineering
    cat("Step 1: Balanced feature engineering...\n")
    preprocessed_data <- balanced_feature_engineering(reference_data, test_data_all)
    
    # Step 2: FCM algorithm
    cat("Step 2: Balanced FCM algorithm...\n")
    fcm_results <- balanced_fcm_algorithm(
      preprocessed_data$reference_data, preprocessed_data$test_data,
      sources, preprocessed_data$group_centers, preprocessed_data$source_features
    )
    
    # Step 3: Post-processing
    cat("Step 3: Balanced post-processing...\n")
    final_results <- balanced_post_processing(fcm_results, preprocessed_data$test_data)
    
    # Key modification: ensure Final_Prediction is based on first probability only
    final_results$Final_Prediction <- apply(final_results[, sources], 1, function(row) {
      if (all(is.na(row)) || all(row == 0)) return(NA)
      max_prob_idx <- which.max(row)
      sources[max_prob_idx]
    })
    
    # Step 4: Visualisation
    cat("Step 4: Visualisation...\n")
    plots <- enhanced_visualizations(final_results, output_dir)
    
    # Step 5: Group analysis
    cat("Step 5: Group analysis...\n")
    group_analysis <- balanced_group_analysis(final_results, output_dir)
    
    # Step 6: Save results
    cat("Step 6: Saving results...\n")
    final_output <- final_results[, c("Sample_ID", "period", "X", "Y", 
                                      "Qinling Region", "South China", "Tongling", "Southwest China",
                                      "Final_Prediction", "Max_Probability", "Probability_Ratio", "Sample_Type")]
    
    write.csv(final_output, file.path(output_dir, "Bronze_Artefact_Ore_Source_Predictions_Balanced.csv"), 
              row.names = FALSE, fileEncoding = "UTF-8")
    
    # Save enhanced features
    enhanced_features <- preprocessed_data$test_data[, c(
      "Sample_ID", "period", "X_std", "Y_std",
      "Dist_to_Qinling", "Dist_to_SouthChina", "Dist_to_Tongling", "Dist_to_Southwest",
      "Angle_Qinling_SouthChina_norm", "Angle_Qinling_Tongling_norm", "Angle_SouthChina_Tongling_norm",
      "Qinling_Score", "SouthChina_Score", "Tongling_Score", "Southwest_Score",
      "Balanced_Zone", "Dominant_Source"
    )]
    
    write.csv(enhanced_features, file.path(output_dir, "Bronze_Artefact_Enhanced_Features.csv"), 
              row.names = FALSE, fileEncoding = "UTF-8")
    
    # Performance summary
    performance_summary <- data.frame(
      Period = unique(final_results$period),
      Total_Samples = as.numeric(table(final_results$period)),
      Qinling_Ratio = as.numeric(by(final_results, final_results$period, 
                                    function(x) mean(x$Final_Prediction == "Qinling Region", na.rm = TRUE) * 100)),
      SouthChina_Ratio = as.numeric(by(final_results, final_results$period, 
                                       function(x) mean(x$Final_Prediction == "South China", na.rm = TRUE) * 100)),
      Tongling_Ratio = as.numeric(by(final_results, final_results$period, 
                                     function(x) mean(x$Final_Prediction == "Tongling", na.rm = TRUE) * 100)),
      Southwest_Ratio = as.numeric(by(final_results, final_results$period, 
                                      function(x) mean(x$Final_Prediction == "Southwest China", na.rm = TRUE) * 100)),
      Avg_Probability = as.numeric(by(final_results, final_results$period, 
                                      function(x) mean(x$Max_Probability, na.rm = TRUE))),
      Overlap_Ratio = as.numeric(by(final_results, final_results$period,
                                    function(x) mean(x$Sample_Type == "Overlap_Region", na.rm = TRUE) * 100))
    )
    
    write.csv(performance_summary, file.path(output_dir, "Performance_Summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    
    # Output detailed statistics
    cat("\n=== Balanced Analysis Complete ===\n")
    cat(sprintf("Total analysed samples: %d\n", nrow(final_results)))
    cat(sprintf("Overlap samples: %d (%.1f%%)\n", 
                sum(final_results$Sample_Type == "Overlap_Region"),
                sum(final_results$Sample_Type == "Overlap_Region") / nrow(final_results) * 100))
    cat(sprintf("Average maximum probability: %.1f%%\n", mean(final_results$Max_Probability, na.rm = TRUE)))
    
    cat("\nOre source distribution by period:\n")
    for (period in unique(final_results$period)) {
      period_data <- final_results[final_results$period == period, ]
      total <- nrow(period_data)
      cat(sprintf("  %s period (%d samples):\n", period, total))
      for (source in sources) {
        count <- sum(period_data$Final_Prediction == source)
        if (count > 0) {
          cat(sprintf("    %s: %d (%.1f%%)\n", source, count, count/total*100))
        }
      }
    }
    
    return(list(
      results = final_results, enhanced_features = enhanced_features,
      plots = plots, group_analysis = group_analysis, performance = performance_summary
    ))
    
  }, error = function(e) {
    cat(sprintf("\nError during analysis: %s\n", e$message))
    return(NULL)
  })
}

# 9. Execute balanced analysis
cat("Executing balanced bronze artefact lead isotope provenance analysis...\n")
bronze_analysis_balanced <- balanced_bronze_analysis(reference_data, bronze_data_clean)

# 10. Check results
if (!is.null(bronze_analysis_balanced)) {
  cat("\n=== Balanced analysis completed successfully! ===\n")
  cat(sprintf("Results saved to: %s\n", output_dir))
  
  cat("\nAll analysis files generated:\n")
  cat("1. Main analysis results:\n")
  cat("   - Bronze_Artefact_Ore_Source_Predictions_Balanced.csv\n")
  cat("   - Bronze_Artefact_Enhanced_Features.csv\n")
  cat("   - Performance_Summary.csv\n")
  cat("   - Ore_Source_Proportion_by_Period.csv\n")
  cat("2. Visualisation charts:\n")
  cat("   - Lead_Isotope_Distribution_Bronze.png\n")
  cat("   - Average_Ore_Source_Probabilities_by_Period.png\n")
  cat("   - Predicted_Ore_Source_Distribution_by_Period.png\n")
  if (sum(bronze_analysis_balanced$results$Sample_Type == "Overlap_Region") > 0) {
    cat("   - Overlap_Samples_Distribution.png\n")
    cat("   - Overlap_Sample_Count_by_Period.png\n")
    cat("   - Overlap_Sample_Types_Distribution.png\n")
  }
  
  # Generate final report summary
  cat("\n", paste0(rep("=", 80), collapse = ""), "\n")
  cat("         FCM Bronze Artefact Ore Source Provenance Analysis - Final Report\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  
  total_samples <- nrow(bronze_analysis_balanced$results)
  overlap_samples <- sum(bronze_analysis_balanced$results$Sample_Type == "Overlap_Region")
  avg_max_prob <- mean(bronze_analysis_balanced$results$Max_Probability, na.rm = TRUE)
  
  cat("\nOverall Statistics:\n")
  cat(paste0(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("Total samples: %d\n", total_samples))
  cat(sprintf("Overlap samples: %d (%.1f%%)\n", overlap_samples, overlap_samples/total_samples*100))
  cat(sprintf("Average maximum probability: %.1f%%\n", avg_max_prob))
  cat(paste0(rep("-", 80), collapse = ""), "\n")
  
  # Main source by period
  cat("\nMain ore source by period (overlap samples merged to most probable source):\n")
  sources <- c("Qinling Region", "South China", "Tongling", "Southwest China")
  for (period in unique(bronze_analysis_balanced$results$period)) {
    period_data <- bronze_analysis_balanced$results[bronze_analysis_balanced$results$period == period, ]
    period_data$Merged_Source <- apply(period_data[, sources], 1, function(row) {
      sources[which.max(row)]
    })
    total <- nrow(period_data)
    main_source <- names(which.max(table(period_data$Merged_Source)))
    main_ratio <- max(table(period_data$Merged_Source)) / total * 100
    cat(sprintf("  %s period: main source is %s (%.1f%%)\n", period, main_source, main_ratio))
  }
  
  # Output file location
  cat("\nResults saved to:\n")
  cat(sprintf("Main directory: %s\n", output_dir))
  
  # Save final report
  final_report_path <- file.path(output_dir, "Final_Analysis_Report.txt")
  sink(final_report_path)
  cat("FCM Bronze Artefact Ore Source Provenance Analysis Report\n")
  cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("I. Data Overview\n")
  cat(sprintf("  Total analysed samples: %d\n", total_samples))
  cat(sprintf("  Reference mine data samples: %d\n", nrow(reference_data)))
  cat(sprintf("  Bronze artefact samples: %d\n", nrow(bronze_data_clean)))
  cat("\nII. Sample Distribution by Period\n")
  period_counts <- table(bronze_analysis_balanced$results$period)
  for (period in names(period_counts)) {
    cat(sprintf("  %s period: %d samples\n", period, period_counts[[period]]))
  }
  cat("\nIII. Overlap Sample Analysis\n")
  cat(sprintf("  Total overlap samples: %d (%.1f%%)\n", overlap_samples, overlap_samples/total_samples*100))
  cat("\nIV. Generated File List\n")
  cat("  1. Data files:\n")
  cat("     - Bronze_Artefact_Ore_Source_Predictions_Balanced.csv\n")
  cat("     - Ore_Source_Proportion_by_Period.csv\n")
  cat("     - Performance_Summary.csv\n")
  cat("  2. Visualisation charts:\n")
  cat("     - Lead_Isotope_Distribution_Bronze.png\n")
  cat("     - Average_Ore_Source_Probabilities_by_Period.png\n")
  cat("     - Predicted_Ore_Source_Distribution_by_Period.png\n")
  if (overlap_samples > 0) {
    cat("     - Overlap_Samples_Distribution.png\n")
    cat("     - Overlap_Sample_Count_by_Period.png\n")
    cat("     - Overlap_Sample_Types_Distribution.png\n")
  }
  sink()
  
  cat(sprintf("\nFinal report saved to: %s\n", final_report_path))
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat("Analysis pipeline completed!\n")
  
} else {
  cat("\nAnalysis failed, please check data format and file paths\n")
}