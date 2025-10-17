##### Libraries #####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(castr)
library(missForest)
library(doParallel)
library(purrr)
library(strucchange)

# Call the functions 

source("code_analyses/functions_pre_clean_plankton.R")
source("code_analyses/functions_to_regularize.R")
source("code_analyses/functions_BP.R")


#### Step 0 :Load the data tables #####
residuals_df     <- read_tsv("data/residuals_df.tsv")
# Load the residuals of the gls regression --> time series without season and no trend

##### Step 1 : Perform BP analysis #####

#### Parameters ####
depth_sel <- 10                  # selected depth
varcol    <- "resid_gls"         # what I analyze
                  

#### Filter by depth ####
resid <- residuals_df %>%
  dplyr::filter(depth == depth_sel) %>%
  dplyr::mutate(target_date = as.Date(target_date)) %>%
  dplyr::arrange(name, target_date)

series <- sort(unique(resid$name))

#### For each variable  ####
for (nm in series) {
  df_i <- resid %>%
    dplyr::filter(name == nm) %>%
    dplyr::filter(is.finite(.data[[varcol]])) %>%
    dplyr::arrange(target_date)
  
  if (nrow(df_i) < min_n) next
  
  bp_i <- calculate_nb_BP_try(df_i, var = varcol, date_col = "target_date")
  p_i  <- plot_BP(df_i, var = varcol, break_dates = bp_i$break_dates, date_col = "target_date")
  
  message(sprintf("Plot %s (n=%d, breaks=%d)", nm, nrow(df_i), length(bp_i$break_dates)))
  print(p_i)  
}

#### With a tibble ####
bp_list <- purrr::map(series, function(nm) {
  df_i <- resid %>%
    dplyr::filter(name == nm, is.finite(.data[[varcol]])) %>%
    dplyr::arrange(target_date)
  if (nrow(df_i) < min_n) {
    return(tibble(name = nm, n = nrow(df_i), n_breaks = 0L, break_dates = list(as.Date(character()))))
  }
  out <- calculate_nb_BP_try(df_i, var = varcol, date_col = "target_date")
  tibble(name = nm, n = nrow(df_i), n_breaks = length(out$break_dates), break_dates = list(out$break_dates))
})

bp_summary <- dplyr::bind_rows(bp_list)
bp_summary
