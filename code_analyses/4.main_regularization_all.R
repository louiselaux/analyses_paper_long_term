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

# Call the functions 

source("code_analyses/functions_pre_clean_plankton.R")
source("code_analyses/functions_to_regularize.R")

#### Step 0 :Load the data tables #####

# Environment hydro

std <- read_tsv("data/std_d.tsv") %>%
  filter(name %in% c("T","S","O","NO3","CHLA","SIOH4","MES"))%>%select(-NOP) 

# Plankton

tablo <- read_tsv(file="data/tablo.tsv")

###### Step 1: For plankton only, pre-clean and format to have it structured as the one of hydro #####

df_plankton <- extract_wanted_level_from_tablo(tablo, level_col="family_like", value_col="concentration")
df_plankton_filled <- fill_missing_dates_taxa(df_plankton)
head(df_plankton_filled)



##### Step 2 : Regularize #####

# Regularize plankton
res_plankton <- run_to_rf(df_long = df_plankton_filled)

# Regularize hydro
res_std <- run_to_rf(df_long=std)


# Inspection of how many points were imputed
res_plankton$df_long_rf %>%
  group_by(var) %>%
  summarise(
    n_total = n(),
    n_rf_imputed = sum(show_red),
    pct_rf_imputed = round(100 * n_rf_imputed / n_total, 2)
  ) %>%
  arrange(desc(pct_rf_imputed))

# How many by RF compared to Linear interpolation
interp_tags <- res_plankton$debug_interp %>%
  filter(size_gap > 0, size_gap <= 5, is.na(value), !is.na(value_final)) %>%
  transmute(target_date, depth, var = name, was_interpolated = TRUE)

res_plankton$df_long_rf %>%
  left_join(interp_tags, by = c("target_date","depth","var")) %>%
  mutate(was_interpolated = replace_na(was_interpolated, FALSE)) %>%
  group_by(var) %>%
  summarise(
    n_total = n(),
    n_interp = sum(was_interpolated),
    n_rf = sum(show_red),
    pct_interp = round(100 * n_interp / n_total, 2),
    pct_rf = round(100 * n_rf / n_total, 2)
  ) %>%
  arrange(desc(pct_rf))

##### Step 3 : Plot to confirm it worked as expected #####
g_hydro_T_9395 <- plot_rf_series(res_std$df_long_rf,
                                 vars  = "T",
                                 years = 1993:1995)
print(g_hydro_T_9395)

g_plk_multi <- plot_rf_series(res_plankton$df_long_rf,
                              vars   = c("Chaetognatha"))
print(g_plk_multi)

##### Step 4 : Save outputs #####
write_tsv(res_plankton$df_long_rf, "data/df_long_rf_plankton.tsv")
write_tsv(res_plankton$final_panel_rf, "data/final_panel_rf_plankton.tsv")

