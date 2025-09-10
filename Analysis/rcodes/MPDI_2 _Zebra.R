###R code for real data analysis, including surrogate variable analysis (SVA) and differential expression (DE) analysis,. 

## codes for FSR variable selection results

# Set directory path
base_dir <- "./SimulationOutsvazebra"

# Find directories and files
model_dirs <- list.dirs(base_dir, recursive = FALSE)[grepl("ModelSize_\\d+_nGene_2000_B_100_alpha0_0.05_ideal_(TRUE|FALSE)$", list.dirs(base_dir, recursive = FALSE))]
all_files <- unlist(lapply(model_dirs, function(dir) list.files(dir, pattern = "^nrep_\\d+\\.rds$", full.names = TRUE)))

# Extract info and process data
file_info <- data.frame(
  file_path = all_files,
  model_size = as.numeric(gsub(".*ModelSize_(\\d+).*", "\\1", dirname(all_files))),
  scenario = ifelse(grepl("TRUE", dirname(all_files)), "TRUE", "FALSE")
)

# Process all groups
results <- list()
for(scenario in c("TRUE", "FALSE")) {
  group_files <- file_info$file_path[file_info$model_size == 1 & file_info$scenario == scenario]
  all_reps <- lapply(group_files, function(file) {
    x <- readRDS(file)
    c(FSP_FSR = x$FSROut["FSR.FSP.OWN.RE"], FSP_SVAall = x$FSROut.SVAallFSR["FSR.FSP.OWN.RE"],
      R_FSR = x$FSROut["FSR.S.OWN.RE"], R_SVAall = x$FSROut.SVAallFSR["FSR.S.OWN.RE"])
  })
  results[[scenario]] <- colMeans(do.call(rbind, all_reps), na.rm = TRUE)
}

# Create plot data
plot_data_reshaped <- data.frame()
methods <- c("FSR_sva", "SVAall_FSR")
method_keys <- c("FSP_FSR", "FSP_SVAall")
r_keys <- c("R_FSR", "R_SVAall")

for(i in 1:2) {
  for(scenario in c("Scenario 1", "Scenario 2")) {
    scenario_key <- ifelse(scenario == "Scenario 1", "TRUE", "FALSE")
    result <- results[[scenario_key]]
    plot_data_reshaped <- rbind(plot_data_reshaped,
                                data.frame(Method = methods[i], Scenario = scenario,
                                           FSP = result[paste0(method_keys[i], ".FSR.FSP.OWN.RE")],
                                           R = result[paste0(r_keys[i], ".FSR.S.OWN.RE")]))
  }
}
# Rename columns in the dataframe
table_data <- plot_data_reshaped
colnames(table_data)[colnames(table_data) == "FSP"] <- "Empirical FSR"
colnames(table_data)[colnames(table_data) == "R"] <- "Average R"

library(ggplot2)
library(tidyr)
library(dplyr)

# reshape the data to long format for ggplot
plot_data_long <- table_data %>%
  pivot_longer(
    cols = c("Empirical FSR", "Average R"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  # Correct order  of the methods
  mutate(Metric = factor(Metric, levels = c("Empirical FSR", "Average R")))

# create plot
ggplot(plot_data_long, aes(x = Method, y = Value)) +
  geom_bar(stat = "identity", fill = "gray40", width = 0.7) +
  facet_grid(
    Metric ~ Scenario,
    scales = "free_y",
    space = "fixed",
    labeller = labeller(
      Scenario = c(
        "Scenario 1" = "Scenario 1: No hidden covariates",
        "Scenario 2" = "Scenario 2: Hidden covariates"
      )
    )
  ) +
  labs(
    title = expression(paste("FSR performance Comparison for ", k[R] == 1)),
    x = "Methods",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

## codes for DE analysis results 
# Set directory path
base_dir <- "./SimulationOutsvazebra"
# List and load .rds files
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", full.names = TRUE)
invisible(lapply(rds_files, function(file) {
  obj_name <- tools::file_path_sans_ext(basename(file))
  assign(obj_name, readRDS(file), envir = .GlobalEnv)
}))
# Calculate means for both scenarios
invisible({
  # Scenario 1 (TRUE)
  means_list_T <- lapply(1, function(i) {
    obj_name <- paste0("ModelSize_", i, "_nGene_2000_B_100_m_3_alpha0_0.05_ideal_TRUE_nSim_1_100")
    x <- colMeans(get(obj_name)$nSimSelCov)
    x[grepl("FDP|NTP|PAUC", names(x))]
  })
  
  # Scenario 2 (FALSE) - FIXED: Now properly stores results
  means_list_F <- lapply(1, function(i) {
    obj_name <- paste0("ModelSize_", i, "_nGene_2000_B_100_m_3_alpha0_0.05_ideal_FALSE_nSim_1_100")
    x <- colMeans(get(obj_name)$nSimSelCov)
    x[grepl("FDP|NTP|PAUC", names(x))]
  })
})
# Extract metrics for both scenarios
extract_metrics <- function(means_list, metric) {
  sapply(means_list, function(x) x[paste0("DEA.", metric)])
}

# Create data with method-metric grouping
plot_data <- data.frame()
methods <- c("FSR_sva", "SVAall_FSR", "SVA0", "Oracle")  # Changed order here
method_suffixes <- c("FSRsva", "SVAall_FSR", "SVA0", "Oracle")

for(i in 1:length(methods)) {
  for(metric in c("FDP", "NTP", "PAUC")) {
    # Add Scenario 1 row
    plot_data <- rbind(plot_data,
                       data.frame(
                         Method = methods[i],
                         Scenario = "Scenario 1",
                         Metric = ifelse(metric == "FDP", "FDR", metric),
                         Value = extract_metrics(means_list_T, paste0(metric, ".", method_suffixes[i]))
                       ))
    
    # Add Scenario 2 row
    plot_data <- rbind(plot_data,
                       data.frame(
                         Method = methods[i],
                         Scenario = "Scenario 2", 
                         Metric = ifelse(metric == "FDP", "FDR", metric),
                         Value = extract_metrics(means_list_F, paste0(metric, ".", method_suffixes[i]))
                       ))
  }
}

# ADD THIS ONE LINE - Convert Method to factor with desired order
plot_data$Method <- factor(plot_data$Method, levels = c("FSR_sva", "SVAall_FSR", "SVA0", "Oracle"))

library(ggplot2)
library(tidyr)
library(dplyr)

# First, ensure correct order of metrics
plot_data_ordered <- plot_data %>%
  mutate(Metric = factor(Metric, levels = c("FDR", "NTP", "PAUC")))

# Create the plot (NO CHANGES HERE)
ggplot(plot_data_ordered, aes(x = Method, y = Value)) +
  geom_bar(stat = "identity", fill = "gray40", width = 0.7) +
  facet_grid(
    Metric ~ Scenario,
    scales = "free_y",
    space = "fixed",
    labeller = labeller(
      Scenario = c(
        "Scenario 1" = "Scenario 1: No hidden covariates",
        "Scenario 2" = "Scenario 2: Hidden covariates"
      )
    )
  ) +
  labs(
    title = expression(paste("Differential expression analysis results for ", k[R] == 1)),
    x = "Methods",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )


## Rcodes for surrogate variable analysis results
library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)

# 1. Load and prepare data
combined_data <- readRDS("./SimulationOutsvazebra/combined_surrogate_frequency.rds")

long_data <- pivot_longer(
  combined_data,
  cols = c(Scenario_1, Scenario_2),
  names_to = "Scenario",
  values_to = "Frequency"
) %>%
  mutate(
    Method_Scenario = paste0(Method, " (", gsub("Scenario_", "S", Scenario), ")")
  ) %>%
  # Filter out SVA0 (S1) plot
  filter(!(Method == "SVA0" & Scenario == "Scenario_1"))

# 2. Create plot
ggplot(long_data, aes(x = factor(Num_Surrogates), y = Frequency)) +
  geom_bar(stat = "identity", fill = "gray40", width = 0.7) +
  facet_grid(~ Method_Scenario) +
  labs(
    title = expression(paste("Bar Plot for Both Scenarios for ", k[R] == 1, " relevant covariates")),
    x = "Number of surrogate variables estimated",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5)
  )


## codes for creating Tables for the average $R^2$ value along with their standard deviation across 100 replication for the three methods.

results='asis'
library(dplyr)
library(tidyr)
library(kableExtra)

# 1. Define covariates
covariate_map <- list(
  `1` = "Batch" 
)

# 2. Main analysis function 
calculate_r_squared <- function(model_sizes = 1) {  
  base_dir <- "./SimulationOutsvazebra"
  all_results <- list()
  
  for (size in model_sizes) {
    cov_cols <- covariate_map[[as.character(size)]]
    if (is.null(cov_cols)) next
    
    model_dir <- file.path(base_dir,
                           sprintf("ModelSize_%s_nGene_2000_B_100_alpha0_0.05_ideal_FALSE", size))
    
    if (!dir.exists(model_dir)) next
    
    rep_files <- list.files(model_dir, pattern = "nrep_\\d+\\.rds$", full.names = TRUE)
    if (length(rep_files) == 0) next
    
    method_results <- list()
    
    for (method in c("FSRsva", "SVA0", "SVAall_FSR")) {
      cov_name <- switch(method,
                         "FSRsva" = "FSRsvacov",
                         "SVA0" = "svacov0",
                         "SVAall_FSR" = "svacovall")
      
      r_squared_values <- sapply(rep_files, function(rep_file) {
        rep_data <- readRDS(rep_file)
        if (!cov_name %in% names(rep_data)) return(NA)
        
        mean(apply(rep_data[[cov_name]], 2, function(x) {
          fit <- lm(x ~ ., data = rep_data$VarCov0[, cov_cols, drop = FALSE])
          summary(fit)$r.squared
        }), na.rm = TRUE)
      })
      
      method_results[[method]] <- data.frame(
        ModelSize = size,
        Method = ifelse(method == "FSRsva", "FSR_sva", method),
        R_squared = r_squared_values,
        stringsAsFactors = FALSE
      )
    }
    
    all_results[[size]] <- bind_rows(method_results)
  }
  
  bind_rows(all_results)
}

# 3. Create summary table with your exact formatting style
create_summary_table <- function(results) {
  # Process results
  summary_data <- results %>%
    group_by(ModelSize, Method) %>%
    summarise(
      Mean_R2 = mean(R_squared, na.rm = TRUE),
      SD_R2 = sd(R_squared, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Display = sprintf("%.3f (%.3f)", Mean_R2, SD_R2)
    ) %>%
    select(ModelSize, Method, Display) %>%
    pivot_wider(
      names_from = Method,
      values_from = Display
    ) %>%
    arrange(ModelSize) %>%
    rename(
      `$k_R$` = ModelSize,
      `FSR\\_sva` = FSR_sva,
      `SVAall\\_FSR` = SVAall_FSR
    )
  
  # Apply your exact table formatting
  kable(summary_data,
        caption = "Average $R^2$ values with corresponding standard errors across 100 replications.",
        format = "latex",
        linesep = "\\hline",
        escape = FALSE,
        align = c("l", "c", "c", "c")) %>%
    kable_styling(
      bootstrap_options = c("striped", "scale_down", "HOLD_position"),
      font_size = 10,
      full_width = FALSE
    ) %>%
    column_spec(1, width = "2cm", border_left = TRUE) %>%
    column_spec(2, width = "3cm") %>%
    column_spec(3, width = "3cm") %>%
    column_spec(4, width = "3cm", border_right = TRUE)
}

# Run the analysis
results <- calculate_r_squared()
final_table <- create_summary_table(results)
final_table




