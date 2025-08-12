###R code for real data analysis, including surrogate variable analysis (SVA) and differential expression (DE) analysis,. 

library(csrnaseq)
FixCov <- csrnaseq::FixCov
VarCov <- csrnaseq::VarCov
counts <- csrnaseq::counts


filtered_id <- apply(counts, 1, function(x) length(x[x>=5]) >= 2)
SimCnt_filtered <-counts[filtered_id,]

lib.size <- apply(SimCnt_filtered, 2, quantile, .75)
y <- t(log2(t(SimCnt_filtered + 0.5)/(lib.size + 1) * 1e+06))

SimCntOut<-FSRAnalysisBS(
counts,
FixCov,
VarCov,
option = "OWN",
B = 100,
m = 3,
lambdamax = 10,
alpha0 = 0.05,
ncores =1,
print.progress= FALSE,
saveall = TRUE
)

BestCovOut <- VarCov[names(SimCntOut$BestCovOut$BestRE)]

# FSR and then sva--------
if(ncol(BestCovOut)>0){
  mod0 <- model.matrix(~., data = BestCovOut)
  mod <- model.matrix(~., data = cbind(FixCov, BestCovOut))
}else{
  mod0 <- model.matrix(~1, data = BestCovOut)
  mod <- model.matrix(~., data = FixCov)
}
FSRsvaout <- sva::sva(dat = y, mod = mod, mod0 = mod0)
FSRsvacov <- data.frame(FSRsvaout$sv)

if(ncol(FSRsvacov) >=1){names(FSRsvacov) <- paste0("sva", 1:ncol(FSRsvacov))

svacov_on_FSR <- cbind(BestCovOut, FSRsvacov)} else{
  svacov_on_FSR<-BestCovOut
  
}


# sva on nothing-----
mod1 <- model.matrix(~., data = FixCov)
svaout0 <- sva::sva(dat = y, mod = mod1)
svacov0 <- data.frame(svaout0$sv)
if(ncol(svacov0) >=1){
  names(svacov0) <- paste0("sva", 1:ncol(svacov0))
  FixCovOut_sva_FSR <- cbind(FixCov,svacov0)
  
}

# sva on everything then FSR
mod0 <- model.matrix(~., data =  VarCov)
mod <- model.matrix(~., data = cbind(FixCov, VarCov))
# svaoutall <- sva::svaseq(dat = SimCnt_filtered, mod = mod, mod0 = mod0)
svaoutall <- sva::sva(dat = y, mod = mod, mod0 = mod0)
svacovall <- data.frame(svaoutall$sv)
if(ncol(svacovall) >=1){
  names(svacovall) <- paste0("sva", 1:ncol(svacovall))
  AllCovOut_sva_FSR <- cbind(FixCov, VarCov, svacovall)
  
}
cbind(VarCov,svacovall)

SVAall_FSRvar<-csrnaseq::FSRAnalysisBS(
  counts,
  FixCov,
  VarCov=cbind(VarCov,svacovall),
  option = "OWN",
  B = 100,
  m = 3,
  lambdamax = 10,
  alpha0 = 0.05,
  ncores = 1,
  print.progress = FALSE,
  saveall = TRUE
)

Allsvaallcov<-cbind(VarCov,svacovall)
BestCovOutFSR_all <- Allsvaallcov[names(SVAall_FSRvar$BestCovOut$BestRE)]


Allcovsva <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,BestCovOutFSR_all ))
deallsva <- sum(csrnaseq:::jabes.q(Allcovsva$pvs[,1])<= 0.05)


# Create the list with all objects (now including the actual values)
results <- list(
  "FSR-SVA" = FSRsvacov,
  "SVA0" = svacov0,
  "SVA-All" = svacovall,
  "svacov-on-FSR" = svacov_on_FSR,
  "FixCovOut-sva-FSR" = FixCovOut_sva_FSR,
  "BestCovOutFSR_all" = BestCovOutFSR_all
)

# Save to RDS file
save_path <- "../analysis/fsr_sva_real_results.rds"
saveRDS(results, file = save_path)


# DE including only Line
lineonly <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = csrnaseq::FixCov)
deline <- sum(csrnaseq:::jabes.q(lineonly$pvs[,1])<= 0.05)
#238

# DE Genes including sva and line
linesva <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,svacov0))
delinesva <- sum(csrnaseq:::jabes.q(linesva$pvs[,1])<= 0.05)
#774

# DE including only the bestcov
Bestcovonly <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov ,BestCovOut))
debestcov <- sum(csrnaseq:::jabes.q(Bestcovonly$pvs[,1])<= 0.05)
#505

# DE including only the bestcov and sva
Bestcovsav <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov ,BestCovOut,FSRsvacov))
debestcovsva <- sum(csrnaseq:::jabes.q(Bestcovsav$pvs[,1])<= 0.05)
#751

# DE genes including all variables
Allcov <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,VarCov))
deallcov <- sum(csrnaseq:::jabes.q(Allcov$pvs[,1])<= 0.05)

# DE genes including all variables and estimated surrogate variables
Allcovsva <- csrnaseq:::VoomPv(counts = csrnaseq::counts, AllCov = cbind(csrnaseq::FixCov,BestCovOutFSR_all ))
deallsva <- sum(csrnaseq:::jabes.q(Allcovsva$pvs[,1])<= 0.05)
#293



### Rcodes to produce Table 2 in the manuscript that contains 8 sets of truly relevant covariates.

library(dplyr)
SelCov <- plyr::laply(list(1, 2, 3, 4, 5, 6, 7, 8), function(i) { # i <- 6
  # Use the full local path to read RDS files
  simRFI <- readRDS(paste0("./../analysis/RealDataOutBS/ModelSize_", i, ".rds"))
  
  lowercase_name <- c("baso", "eosi", "lymp", "mono", "neut")
  best_re_name <- names(simRFI$BestCovOut$BestER)
  best_re_name[best_re_name %in% lowercase_name] <- stringr::str_to_title(best_re_name[best_re_name %in% lowercase_name])
  
  c(i, paste(rev(best_re_name), collapse = ", "))
})

colnames(SelCov) <- c("Number of relevant covariates $k_R$", "Relevant covariates")

SelCov %>% 
  kableExtra::kable(
    caption = "Simulation scenarios corresponding to eigth sets of truly relevant covariates.", 
    align = c("c", "r", "r"), 
    escape = FALSE, 
    format = "latex"
  ) %>% 
  kableExtra::kable_styling(font_size = 11) %>%
  kableExtra::column_spec(1, border_left = TRUE) %>%
  kableExtra::column_spec(2, border_right = TRUE)


### The rcode to create the line plot for the empirical FSR and average R for two different scenarios for different $k_R$ over 100 replications.

# Set directory path
base_dir <- "../analysis/SimulationOutsvaruv"
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
invisible({
  for(size in 1:8) {
    for(scenario in c("TRUE", "FALSE")) {
      group_files <- file_info$file_path[file_info$model_size == size & file_info$scenario == scenario]
      all_reps <- lapply(group_files, function(file) {
        x <- readRDS(file)
        c(FSP_FSR = x$FSROut["FSR.FSP.OWN.RE"], FSP_SVAall = x$FSROut.SVAallFSR["FSR.FSP.OWN.RE"],
          R_FSR = x$FSROut["FSR.S.OWN.RE"], R_SVAall = x$FSROut.SVAallFSR["FSR.S.OWN.RE"])
      })
      results[[paste0("Size_", size, "_", scenario)]] <- colMeans(do.call(rbind, all_reps), na.rm = TRUE)
    }
  }
})

# Create plot data
plot_data_reshaped <- data.frame()
invisible({
  for(size in 1:8) {
    for(scenario in c("Scenario 1", "Scenario 2")) {
      key <- paste0("Size_", size, "_", ifelse(scenario == "Scenario 1", "TRUE", "FALSE"))
      result <- results[[key]]
      plot_data_reshaped <- rbind(plot_data_reshaped,
                                  data.frame(ModelSize = size, Scenario = scenario, Method = "FSR_sva", 
                                             FSP = result["FSP_FSR.FSR.FSP.OWN.RE"], R = result["R_FSR.FSR.S.OWN.RE"]),
                                  data.frame(ModelSize = size, Scenario = scenario, Method = "SVAall_FSR", 
                                             FSP = result["FSP_SVAall.FSR.FSP.OWN.RE"], R = result["R_SVAall.FSR.S.OWN.RE"]))
    }
  }
})

# Plot
library(tidyr)
library(ggplot2)
plot_data_long <- pivot_longer(plot_data_reshaped, cols = c(FSP, R), names_to = "Metric", values_to = "Value")

ggplot(plot_data_long, aes(x = ModelSize, y = Value, color = Method, shape = Method)) +
  geom_line() + 
  geom_point(size = 3) +
  geom_hline(
    data = data.frame(Metric = "FSP", yint = 0.05),
    aes(yintercept = yint),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  facet_grid(
    Metric ~ Scenario, 
    scales = "free_y",
    labeller = labeller(
      Scenario = c(
        "Scenario 1" = "Scenario 1: No hidden covariates", 
        "Scenario 2" = "Scenario 2: Hidden covariates"
      ),
      Metric = c("FSP" = "Empirical FSR", "R" = "Average R")
    )
  ) +
  scale_x_continuous(breaks = 1:8) + 
  labs(x = "Number of Relevant Covariates", y = NULL) +
  scale_color_manual(values = c("FSR_sva" = "#FFB366", "SVAall_FSR" = "#FFD700")) +
  scale_shape_manual(values = c(17, 8)) + 
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray90"), 
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"), 
    strip.text = element_text(size = 10),
    legend.position = "right"
  )


### The rcode to produce the figure that presents the performance of differential expression analysis of the four methods

# Set directory path
base_dir <- "../analysis/SimulationOutsvaruv"
# List and load .rds files
rds_files <- list.files(path = base_dir, pattern = "\\.rds$", full.names = TRUE)
invisible(lapply(rds_files, function(file) {
  obj_name <- tools::file_path_sans_ext(basename(file))
  assign(obj_name, readRDS(file), envir = .GlobalEnv)
}))
# Calculate means for both scenarios
invisible({
  # Scenario 1 (TRUE)
  means_list_T <- lapply(1:8, function(i) {
    obj_name <- paste0("ModelSize_", i, "_nGene_2000_B_100_m_3_alpha0_0.05_ideal_TRUE_nSim_1_100")
    x <- colMeans(get(obj_name)$nSimSelCov)
    x[grepl("FDP|NTP|PAUC", names(x))]
  })
  
  # Scenario 2 (FALSE) - FIXED: Now properly stores results
  means_list_F <- lapply(1:8, function(i) {
    obj_name <- paste0("ModelSize_", i, "_nGene_2000_B_100_m_3_alpha0_0.05_ideal_FALSE_nSim_1_100")
    x <- colMeans(get(obj_name)$nSimSelCov)
    x[grepl("FDP|NTP|PAUC", names(x))]
  })
})
# Extract metrics for both scenarios
extract_metrics <- function(means_list, metric) {
  sapply(means_list, function(x) x[paste0("DEA.", metric)])
}
# Scenario 1
metrics_1 <- lapply(c("FDP", "NTP", "PAUC"), function(m) {
  data.frame(
    ModelSize = 1:8,
    Scenario = "Scenario 1",
    Method = rep(c("FSR_sva", "SVA0", "SVAall_FSR", "Oracle"), each = 8),
    Metric = ifelse(m == "FDP", "FDR", m),
    Value = c(extract_metrics(means_list_T, paste0(m, ".FSRsva")),
              extract_metrics(means_list_T, paste0(m, ".SVA0")),
              extract_metrics(means_list_T, paste0(m, ".SVAall_FSR")),
              extract_metrics(means_list_T, paste0(m, ".Oracle")))
  )
})
# Scenario 2
metrics_2 <- lapply(c("FDP", "NTP", "PAUC"), function(m) {
  data.frame(
    ModelSize = 1:8,
    Scenario = "Scenario 2",
    Method = rep(c("FSR_sva", "SVA0", "SVAall_FSR", "Oracle"), each = 8),
    Metric = ifelse(m == "FDP", "FDR", m),
    Value = c(extract_metrics(means_list_F, paste0(m, ".FSRsva")),
              extract_metrics(means_list_F, paste0(m, ".SVA0")),
              extract_metrics(means_list_F, paste0(m, ".SVAall_FSR")),
              extract_metrics(means_list_F, paste0(m, ".Oracle")))
  )
})
# Combine all data
plot_data <- do.call(rbind, c(metrics_1, metrics_2))
# Create plot
library(ggplot2)
ggplot(plot_data, aes(x = ModelSize, y = Value, color = Method, shape = Method)) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(
    data = subset(plot_data, Metric == "FDR"),  # Only for FDR panels
    aes(yintercept = 0.05), 
    linetype = "dashed", 
    color = "black",
    linewidth = 0.5
  ) +
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
  scale_x_continuous(breaks = 1:8) +
  labs(x = "Number of Relevant Covariates", y = NULL) +
  scale_color_manual(
    values = c(
      "FSR_sva" = "#FFB366",    
      "SVA0" = "#66B2FF",
      "SVAall_FSR" = "#FFD700",
      "Oracle" = "#808080"
    ),
    breaks = c("FSR_sva", "SVAall_FSR", "SVA0", "Oracle")
  ) +
  scale_shape_manual(
    values = c(17, 4, 18, 8),
    breaks = c("FSR_sva", "SVAall_FSR", "SVA0", "Oracle")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10),
    legend.position = "right"
  )


### Rcode to produce table2
library(knitr)
library(kableExtra)

# Create the data frame with formatted column names
table_data <- data.frame(
  "k\\textsubscript{R}" = 1:8,
  "Relevant Covariates" = c(
    "Mono",
    "Concb, Mono",
    "Neut, Concb, Mono",
    "Block, Neut, Concb, Mono",
    "RINa, Block, Neut, Concb, Mono",
    "Baso, RINa, Block, Neut, Concb, Mono",
    "Lymp, Baso, RINa, Block, Neut, Concb, Mono",
    "RFI, Lymp, Baso, RINa, Block, Neut, Concb, Mono"
  ),
  "Selected by FSR (S2)" = c(
    "RINa",
    "RINa, Concb",
    "RINa, Concb",
    "Block, Concb",
    "RINa, Block, Concb",
    "RINa, Block, Concb",
    "RINa, Block, Concb",
    "RINa, Block, Concb"
  ),
  check.names = FALSE  # Preserves LaTeX formatting
)

# Generate the table
kable(table_data, 
      format = "latex",
      booktabs = FALSE,
      align = c("c", "l", "l"),
      escape = FALSE,
      caption = "Eight simulation scenarios with FSR selection results (S2)",
      linesep = "\\hline") %>%
  kable_styling(
    font_size = 10,
    full_width = FALSE,
    latex_options = c("HOLD_position", "repeat_header")
  ) %>%
  column_spec(1, width = "1.5cm", border_left = TRUE) %>%
  column_spec(2, width = "7.5cm") %>%
  column_spec(3, width = "3.5cm", border_right = TRUE)



#### Rcode to create the combined_modelsize_data.rds 
library(dplyr)
library(knitr)
library(purrr)
library(kableExtra)
library(tidyr)

# Function to process frequency tables
generate_combined_frequency_tables <- function(model_sizes) {
  base_dir <- "../analysis"
  
  all_results <- map_dfr(c(FALSE, TRUE), function(scenario) {
    scenario_label <- ifelse(scenario, "Scenario 1", "Scenario 2")
    
    map_dfr(model_sizes, function(model_size) {
      sim_file <- sprintf("4-Simulation-sva_out/4-Simulation-sva_4163060_%s_%s.out", 
                          model_size, toupper(scenario))
      sim_results <- tryCatch(
        readLines(file.path(base_dir, sim_file)),
        error = function(e) {
          message("Error reading file for model size ", model_size, " scenario ", scenario, ": ", e$message)
          return(NULL)
        }
      )
      
      if (is.null(sim_results)) return(NULL)
      
      surrogate_lines <- grep("significant surrogate variables", sim_results, value = TRUE)
      num_surrogates <- ifelse(
        grepl("No significant surrogate variables", surrogate_lines),
        0,
        as.numeric(gsub(".*Number of significant surrogate variables is: ", "", surrogate_lines))
      )
      
      data.frame(
        Replication = rep(1:(length(num_surrogates)/3), each = 3),
        Method = rep(c("FSR_sva", "SVA0", "SVAall_FSR"), length.out = length(num_surrogates)),
        Num_Surrogates = num_surrogates,
        ModelSize = model_size,
        Scenario = scenario_label
      )
    }) %>% 
      filter(!is.na(Num_Surrogates))
  })
  
  if (nrow(all_results) > 0) {
    all_results %>%
      group_by(Scenario, ModelSize, Method, Num_Surrogates) %>%
      summarise(Frequency = n(), .groups = "drop") %>%
      arrange(Scenario, ModelSize, Method, Num_Surrogates)
  } else {
    return(NULL)
  }
}

create_formatted_tables <- function(data) {
  if (is.null(data) || nrow(data) == 0) {
    cat("No frequency data available for the specified model sizes.")
    return()
  }
  
  # Generate all tables first
  tables <- data %>%
    split(.$ModelSize) %>%
    map(~ {
      current_model_size <- unique(.x$ModelSize)
      
      scenario_tables <- .x %>%
        split(.$Scenario) %>%
        map(~ {
          .x %>%
            select(Method, Num_Surrogates, Frequency) %>%
            arrange(Method, Num_Surrogates)
        })  # <-- Closing parenthesis for inner map(~ {...})
      
      combined <- scenario_tables[[1]] %>%
        full_join(scenario_tables[[2]], 
                  by = c("Method", "Num_Surrogates"),
                  suffix = c("_S1", "_S2")) %>%
        replace_na(list(Frequency_S2 = 0, Frequency_S1 = 0)) %>%
        arrange(Method, Num_Surrogates)
      
      # Save the combined dataset as RDS file
      base_dir <- "base_dir <- "../analysis""
      save_file <- file.path(base_dir, paste0("combined_frequency_modelsize_", current_model_size, ".rds"))
      saveRDS(combined, save_file)
    })  
}

library(ggplot2)
library(gridExtra)
library(grid)

# Set the base path
base_path <- "../analysis/"

# Create a list of all 8 dataframes
df_list <- list()

# Loop through each model size (1 to 8) and load the dataframes
for (i in 1:8) {
  file_name <- paste0("combined_frequency_modelsize_", i, ".rds")
  file_path <- paste0(base_path, file_name)
  
  # Load the dataframe and add model size identifier
  df_list[[i]] <- readRDS(file_path) %>% 
    dplyr::mutate(Model_Size = paste0(i))
}

# Combine all dataframes into one
combined_df <- dplyr::bind_rows(df_list)

# Save the combined dataframe
saveRDS(combined_df, file = paste0(base_path, "combined_modelsize_data.rds"))

# Print the structure to verify
str(combined_df)





### Rcode to create the Bar plot for the frequency for the number of surrogate variables estimated via three methods FSR_sva, SVA0 and SVAall_FSR

library(tidyr)
library(ggplot2)
library(dplyr)  
library(magrittr)  
# 1. Load and prepare data
combined_data <- readRDS("../analysis/combined_modelsize_data.rds")

long_data <- pivot_longer(
  combined_data,
  cols = c(Frequency_S1, Frequency_S2),
  names_to = "Scenario",
  values_to = "Frequency"
) %>%
  mutate(
    Model_Size = gsub("modelsize_", "", Model_Size),
    Method_Scenario = case_when(
      Method == "SVA0" & Scenario == "Frequency_S1" ~ NA_character_,
      TRUE ~ paste0(Method, " (", gsub("Frequency_", "", Scenario), ")")
    ),
    # Create new label column with proper formatting
    Model_Size_Label = paste0("k[R]==", Model_Size)
  ) %>%
  filter(!is.na(Method_Scenario))

# 2. Create plot
ggplot(long_data, aes(x = factor(Num_Surrogates), y = Frequency)) +
  geom_bar(stat = "identity", fill = "gray40", width = 0.7) +
  facet_grid(
    Model_Size_Label ~ Method_Scenario,
    labeller = label_parsed  # This interprets the expressions
  ) +
  labs(
    title = expression(paste("Bar Plot for Both Scenarios for 8 sets of ", k[R], " relevant covariates")),
    x = "Number of surrogate variables estimated",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.5, "lines"),
    strip.text = element_text(size = 6),
    axis.text.x = element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(hjust = 0.5)
  )


### Rcode to produce the Rsquared table 

library(dplyr)
library(tidyr)
library(kableExtra)

# 1. Define covariates
covariate_map <- list(
  `1` = "mono",
  `2` = "mono",
  `3` = c("neut", "mono"),
  `4` = c("neut", "mono"),
  `5` = c("neut", "mono"),
  `6` = c("baso", "neut", "mono"),
  `7` = c("lymp", "baso", "neut", "mono"), 
  `8` = c("RFI", "lymp", "baso", "neut", "mono")
)

# 2. Main analysis function (unchanged)
calculate_r_squared <- function(model_sizes = 1:8) {
  base_dir <- "./../analysis"
  all_results <- list()
  
  for (size in model_sizes) {
    cov_cols <- covariate_map[[as.character(size)]]
    if (is.null(cov_cols)) next
    
    model_dir <- file.path(base_dir,
                           sprintf("SimulationOutsvaruv/ModelSize_%s_nGene_2000_B_100_alpha0_0.05_ideal_FALSE", size))
    
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
        caption = "Average $R^2$ values with standard deviations across 100 replicates ",
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






