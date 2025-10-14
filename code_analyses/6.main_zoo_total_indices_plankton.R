##### Libraries
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(vegan) 

# Load the functions
source("code_analyses/functions_gls.R") 
source("code_analyses/functions_to_regularize.R")

##### What was already regularized, at the family level
final_panel <- read_tsv("data/final_panel_rf_plankton.tsv")
clim_wk     <- read_tsv("data/clim_wk_smooth.tsv")


##### Step 1 : Prepare the community matrix ######
depth_one <- 10
fp <- final_panel %>% dplyr::filter(depth == depth_one)


# All the family-like columns
fam_cols <- setdiff(names(fp), c("target_date","depth"))

##### Come back to raw value 
A_long <- fp %>%
  pivot_longer(all_of(fam_cols), names_to = "family", values_to = "logv") %>%
  mutate(value = pmax(10^logv - 1, 0))  

# Community matrix
comm_mat <- A_long %>%
  select(target_date, family, value) %>%
  pivot_wider(names_from = family, values_from = value, values_fill = 0) %>%
  arrange(target_date)

# Only taxonomic groups
comm_m <- as.matrix(select(comm_mat, -target_date))

# Zoo total
zoo_total_raw <- rowSums(comm_m, na.rm = TRUE)
zoo_total_log <- log10(zoo_total_raw + 1)

# Indices
Srich       <- specnumber(comm_m)                    # Richness
Shannon     <- diversity(comm_m, index = "shannon")  # H'
Pielou      <- ifelse(Srich > 0, Shannon / log(Srich), NA_real_)   # J'
Simpson     <- diversity(comm_m, index = "simpson")  # 1 - D
InvSimpson  <- diversity(comm_m, index = "invsimpson")

# Table of indices
idx_long <- tibble(
  target_date   = comm_mat$target_date,
  depth         = depth_one,
  ZooTot_log    = zoo_total_log,
  ZooTot_raw    = zoo_total_raw,
  Srich         = Srich,
  Shannon       = Shannon,
  Pielou        = Pielou,
  Simpson       = Simpson,
  InvSimpson    = InvSimpson
)


##### Step 2: Remove seasonality #####
panel_idx <- idx_long %>%
  select(target_date, depth,
         ZooTot_log, Srich, Shannon, Pielou, Simpson, InvSimpson)

vars_idx <- c("ZooTot_log","Srich","Shannon","Pielou","Simpson","InvSimpson")

idx_long2 <- panel_idx %>%
  tidyr::pivot_longer(dplyr::all_of(vars_idx), names_to = "name", values_to = "value")


clim_idx <- weekly_clim_smooth(idx_long2)
#
dstl_idx <- build_dstl(
  final_panel = panel_idx,
  clim_wk     = clim_idx,                
  vars        = vars_idx,
  depths_keep = depth_one
)


##### Step 3 : Gls regression ######

stats_idx <- fit_trends(dstl_idx, depth_one = as.character(depth_one),
                        acf_threshold = 0.2, min_n = 5)

stats_idx %>%
  mutate(model_used = case_when(
    gls_cor.struct == "corAR1"  ~ "AR1",
    gls_cor.struct == "corARMA" ~ "AR2",
    TRUE                        ~ "GLS"
  )) %>%
  select(name, gls_slope, gls_p.value, model_used, mk_p.value, start, end) %>%
  print(n = Inf)

##### Step 4 : Plots ######
p_idx <- plot_trends(dstl_idx, stats_idx, vars_show = vars_idx)
print(p_idx)

