# Header
# Objectives:
#   1) Choose optimal Embedding dimension for each column variable
#   2) CCM for each pair of variables (columns, target) at Tp = 0, -1, -2
#   3) Test for convergence + seasonal surrogates
#   4) Identify causal links and look at rho depending on Tp
#
# Interpretation:
#   - columns xmap target => target "causes" columns



##### load the libraries #####
library(doParallel)
library(doRNG)
library(foreach)
library(dplyr)
library(tibble)
library(rEDM)
library(tidyverse)
library(igraph)
library(ggplot2)
library(readr)
library(ggplot2)
library(furrr)

# Load function
source("code_CCM/function_ccm_plot.R")

##### Load the data table#####

ccm_data_scaled <- read_tsv(file="code_CCM/ccm_data_scaled.tsv") # this data set contains the environmental and plankton variables that are scaled

head(ccm_data_scaled)

#Load the wind datasets as well

wind <- read_tsv(file="code_CCM/wind_scaled.tsv")

head(wind)

#Load the couranto data set
currents <- read_tsv(file="code_CCM/currents_scaled.tsv")
head(currents)

# Join the three : the three data sets are already scaled and have the values at the same dates, with a bimonthly resolution

ccm_data_scaled_2 <- cbind(ccm_data_scaled, wind, currents)
head(ccm_data_scaled_2)
tail(ccm_data_scaled_2)

##### Step 1: Find the best embedding dimension #####

#  Create an empty data frame to store the results
embedding_dim_tot <- data.frame(variables=colnames(ccm_data_scaled_2),
                                E.rho=NA,E_best=NA)
embed_stats_tot<- data.frame()


n_obs <- nrow(ccm_data_scaled_2)
lib_range  <- c(1, n_obs)
pred_range <- c(1, n_obs)

# For each time series variable
for (j in colnames(ccm_data_scaled_2)) {
  E.rho <- EmbedDimension(
    dataFrame = ccm_data_scaled_2,
    columns = j,
    target = j,
    showPlot = TRUE,
    noTime = TRUE,
    maxE = 10,
    lib = lib_range,
    pred = pred_range
  )
  
  # Find best E
  best_E <- E.rho$E[which.max(E.rho$rho)][1] # 1 or last?
  
  # Add to the dataframe
  embed_stats_tot <- embed_stats_tot <- rbind(embed_stats_tot,
                                              data.frame(variables = j, E_best = best_E))
  
}


embed_stats_tot

##### Step 2 :  Get the rho at max lib size for each Tp #####

# Load the time sequence data to be able to generate shuffled time series
ref<- read_tsv(file="code_CCM/ref.tsv") # this contains the dates 
env_window <- ccm_data_scaled_2
env_window$date <- ref$target_date

# Generate it by shuffling per week for each variable
env_test_shuffle <- env_window %>% mutate(week=week(date))%>%select(-date)%>%group_by(week) %>%
  mutate(across(-any_of("week"), ~ sample(.x, dplyr::n()), .names = "{.col}_shuffled"))

# Separate data in two parts 
df_real <- env_window
df_shuf <- env_test_shuffle 


# Function to perform CCM and all the tests to assume convergence 


# causal_tests_shuffle():
#  For a couple of variable (columns, target) at a given Tp:
#     - uses E_best(columns) from embed_stats_tot
#     - does a CCM  (columns xmap target) on 4 libSizes
#     - tests for :
#         * convergence: rho(50%) > rho(20%) + strict monotony increase (rho) vs LibSize
#         * non-seasonality: real (lib max) > series "weekly shuffled" of columns
#     - output : a line with all the stats (p-values, rhos, flags)


causal_tests_shuffle <- function(target, columns, df_real, df_shuf, n = 100, Tp = 0, tau = -1) {
  # Shuffle columned 
  shuf_name <- paste0(columns, "_shuffled")
  if (!all(c(columns, target) %in% names(df_real))) stop("columns/target absent of df_real")
  if (!shuf_name %in% names(df_shuf)) stop(sprintf("Colonne absent of df_shuf : %s", shuf_name)) 
  
  # E adapted to columns
  Ei <- embed_stats_tot$E_best[embed_stats_tot$variables == columns][1]
  if (is.na(Ei)) stop(sprintf("No best E'%s' in embed_stats_tot", columns))
  
  # Lib sizes
  n_obs <- nrow(df_real)
  lib_sizes <- c(round(0.20*n_obs), round(0.30*n_obs), round(0.50*n_obs), round(0.90*n_obs))
  
  #CCM ref : columns xmap target
  ref_ccm <- rEDM::CCM(
    dataFrame   = df_real,
    columns     = columns,
    target      = target,
    E           = as.integer(Ei),
    Tp          = Tp,
    tau         = tau,
    libSizes    = lib_sizes,
    sample      = n,
    includeData = TRUE,
    showPlot    = FALSE,
    noTime      = TRUE
  )
  ref1 <- ref_ccm$CCM1_PredictStat  # columns:target
  ref2 <- ref_ccm$CCM2_PredictStat  # target:columns
  
  # Check that rho increases linearly with libsize 
  
  # 1) Convergence (50% > 20%)
  conv_test <- t.test(
    dplyr::filter(ref1, LibSize == lib_sizes[3])$rho,
    dplyr::filter(ref1, LibSize == lib_sizes[1])$rho,
    alternative = "greater"
  )
  
  # Monotonously increasing
  mono_tbl <- ref1 |>
    dplyr::group_by(LibSize) |>
    dplyr::summarise(mu = mean(rho, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(LibSize) 
  
  # Monotony strict
  mono_strict <- all(diff(mono_tbl$mu) > 0)
  
  
  # 2) Shuffle surrogate : columns_shuffled xmaps target
  
  shuf_ccm <- rEDM::CCM(
    dataFrame   = df_shuf,
    columns     = shuf_name,       # ex: "gelatinous_pred_shuffled xmaps temperature"
    target      = target,
    E           = as.integer(Ei),
    Tp          = Tp,
    tau         = tau,
    libSizes    = lib_sizes[4],
    sample      = n,
    includeData = TRUE,
    showPlot    = FALSE,
    noTime      = TRUE
  )
  shuf1 <- shuf_ccm$CCM1_PredictStat  # shuffled:target 
  
  # Observed > shuffled
  non_season <- t.test(
    x = dplyr::filter(ref1, LibSize == lib_sizes[4])$rho, # at maximum lib size
    y = shuf1$rho,
    alternative = "greater"
  )
  surr_greater_obs <- mean(shuf1$rho > non_season$estimate[1])
  
  # Get a tibble at the end with everything 
  
  dplyr::tibble(
    columns = columns, target = target,
    E_used = as.integer(Ei), Tp = Tp, tau = tau,
    # convergence
    converg_stat = unname(conv_test$statistic),
    converg_p    = conv_test$p.value,
    converg_lib_low  = unname(conv_test$estimate[2]),  # ~20%
    converg_lib_high = unname(conv_test$estimate[1]),  # ~80%
    mono_strict=mono_strict,
    # real vs shuffle
    non_season_stat       = unname(non_season$statistic),
    non_season_p          = non_season$p.value,
    obs_mean_lib_high     = unname(non_season$estimate[1]),
    shuffle_mean_lib_high = unname(non_season$estimate[2]),
    `% surr>obs` = surr_greater_obs
  )
}


# Define the parameters you want to test

# Pair of variables
vars_all <- names(df_real)
vars <- vars_all[!grepl("_shuffled$", vars_all) & !vars_all %in% c("date","week")]
tps <- c(0, -1, -2) # We test it for three Tp

# All possible combinations
combin <- tidyr::crossing(
  columns = vars,
  target  = vars,
  Tp      = tps
) %>%
  dplyr::filter(columns != target)

# Parallelize the code otherwise it is too long
future::plan(multisession, workers = 21)

# Do it for each couple of variables and for the three Tp values you want to test
res_all <- furrr::future_pmap_dfr(
  list(combin$target, combin$columns, combin$Tp),
  ~ causal_tests_shuffle(target = ..1, columns = ..2,
                         df_real = df_real, df_shuf = df_shuf,
                         n = 100, Tp = ..3, tau = -1),
  .options = furrr::furrr_options(seed = TRUE)
)

# Stop parallelizing
future::plan(sequential)

# Get the data table
res_all_tp <- res_all
res_all_tp %>% write_tsv(file="code_CCM/res_all_tp.tsv")


##### Step 3: Plot Tp vs rho based on that #####

# Get significant causal links

threshold <- 0.01

res_all_flag <- res_all %>%
  dplyr::mutate(
    signif = converg_p < threshold & # convergence has to be true
      non_season_p < threshold & # t.test between surrogates and real rhos at max lib size has to be significant
      `% surr>obs` < threshold &
      mono_strict # has to increase monotonically with lib size
  )

table(res_all_flag$signif) # If all those conditions are fulfilled, then signif= TRUE

# signif = TRUE if :
#   - convergence(libSize) is significant
#   - pred at max lib size > surrogates
#   => causal link (in the direction target -> columns, meaning columns xmaps target)


# Look at the links
links <- res_all_flag %>%
  dplyr::filter(signif) %>%              
  dplyr::arrange(converg_p)

links

# Look at the links direction
links_dir <- links %>%
  dplyr::transmute(
    cause  = target,      # target xmap columns => target cause columns
    effect = columns,
    weight = converg_lib_high,
    Tp     = Tp
  )

links_dir



##### Last step : perform the plot #####

# Look only at pairs of variables for which at least there is a significant link at at least one Tp value
pairs_sig <- res_all_flag %>%
  dplyr::group_by(columns, target) %>%
  dplyr::summarise(any_sig = any(signif), .groups = "drop") %>%
  dplyr::filter(any_sig)

res_all_sig_pairs <- res_all_flag %>%
  dplyr::semi_join(pairs_sig, by = c("columns", "target"))


# Focus on some variables 
vars_focus_targ <- c("S")
#vars_focus_col<- c(")

res_focus <- res_all_flag %>%
  dplyr::filter(#columns %in% vars_focus_col,
                target  %in% vars_focus_targ)

# Plot it 
ggplot(
  res_focus,
  aes(x = Tp, y = obs_mean_lib_high, group = 1)
) +
  geom_line(alpha = 0.4) +
  geom_point(
    aes(fill = signif),
    shape = 21,
    color = "black",
    size  = 2
  ) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white")) +
  facet_wrap(~ columns, scales = "free_y", ncol = 4) +
  theme_minimal()

##### Check it is consistent by looking at one plot more precisely#####
# The function make_ccm_plot comes from the code source functions_ccm_plot.R
# Try with one

lib_sizes <- c(100, 200, 300, 400,  500, 600, 700)
  
make_ccm_plot(
  columns = "curr_speed_mean_around",
  target  = "S",
  df_real = df_real,
  df_shuf = df_shuf,
  embed_stats_tot = embed_stats_tot,
  lib_sizes = lib_sizes,
  n = 100,
  Tp = -1,
  tau = -1
)
