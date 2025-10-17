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

# HPLC

hplc <- read_tsv("data/hplc_coeffs.tsv")

# CTD

final_point_B_ctd <- read_tsv("data/final_point_B_ctd.tsv")

###### Step 1: For plankton only, pre-clean and format to have it structured as the one of hydro #####

df_plankton <- extract_wanted_level_from_tablo(tablo, level_col="family_like", value_col="concentration")
df_plankton_filled <- fill_missing_dates_taxa(df_plankton)
head(df_plankton_filled)



##### Step 2 : Regularize #####

# Regularize plankton
res_plankton <- run_to_clim(df_long = df_plankton_filled)

# Regularize hydro
res_std <- run_to_rf(df_long=std)

# Regularize HPLC
res_hplc <- run_to_clim(
  df_long        = hplc_long_multi,
  start_date     = as.Date("2012-02-01"),
  by_days        = 14,
  tolerance_days = 3,
  max_gap_interp = 3
)


# Regularize CTD

ctd <- final_point_B_ctd %>%
       mutate(date = ymd(date)) %>%
       pivot_longer(cols = "stratif_index", names_to = "var", values_to = "value") %>%
       mutate(depth =10) %>%
       select(date, var, depth, value)

ctd <- ctd %>%
      dplyr::rename(name = var) %>%  
      dplyr::mutate(depth = 10)

# When to start the grid to make it match with plankton
anchor <- as.Date("1967-01-04")
grid_1995 <- seq(anchor, as.Date("1995-12-31"), by = "14 days")
grid_1995 <- grid_1995[year(grid_1995) == 1995]


res_ctd <- run_to_clim(
  df_long=ctd,
  start_date = as.Date("1995-01-18"),
  by_days=14,
  tolerance_days=3,
  max_gap_interp=5)

##### Step 3 : Plot to confirm it worked as expected #####
g_hydro_T_9395 <- plot_rf_series_2(res_std$df_long_rf,
                                 vars  = "T",
                                 years = 1993:1995)
print(g_hydro_T_9395)

g_plk_multi <- plot_rf_series_2(res_plankton$df_long_rf,
                              vars   = c("Chaetognatha"))
print(g_plk_multi)

##### Step 4 : Save outputs #####
write_tsv(res_plankton$df_long_rf, "data/df_long_rf_plankton.tsv")
write_tsv(res_plankton$final_panel_rf, "data/final_panel_rf_plankton.tsv")

##### Step 5 : A better plot #####
rf_plot_data <- res_std$df_long_rf %>%
  dplyr::left_join(
    res_std$debug_interp %>%
      dplyr::select(target_date, name, depth, value_final) %>%
      dplyr::rename(value_final_debug = value_final),
    by = c("target_date" = "target_date", "var" = "name", "depth" = "depth")
  ) %>%
  dplyr::mutate(
    is_rf = is.na(value_final_debug) 
  )


rf_plot_data %>%filter(year(target_date)%in%c(1990:1996))%>%
  dplyr::filter(var == "T", depth == 10) %>%
  ggplot(aes(x = target_date, y = final)) +
  geom_line() +
  geom_point(aes(color = is_rf), size = 1.3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                     labels = c("FALSE" = "observé / linéaire", "TRUE" = "RF imputé"),
                     name = "") +
  theme_bw() +
  labs(
    title = "",
    x = "Date", y = "Concentration",
    subtitle = ""
  )

# Plot HPLC
rf_plot_data <- res_hplc$df_long_rf %>%
  dplyr::left_join(
    res_hplc$debug_interp %>%
      dplyr::select(target_date, name, depth, value_final) %>%
      dplyr::rename(value_final_debug = value_final),
    by = c("target_date" = "target_date", "var" = "name", "depth" = "depth")
  ) %>%
  dplyr::mutate(
    is_rf = is.na(value_final_debug) 
  )

rf_plot_data %>%filter(year(target_date)%in%c(2012:2023))%>%
  dplyr::filter(var == "pico_chla", depth == 10) %>%
  ggplot(aes(x = target_date, y = final)) +
  geom_line() +
  geom_point(aes(color = is_rf), size = 1.3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                     labels = c("FALSE" = "observé / linéaire", "TRUE" = "RF imputé"),
                     name = "") +
  theme_bw() +
  labs(
    title = "",
    x = "Date", y = "Concentration",
    subtitle = ""
  )

# Small tests
res_plankton$debug_interp%>%mutate(year=year(target_date))%>%filter(is.na(value_final))%>%select(year,target_date)%>%distinct()%>%group_by(year)%>%summarize(count=n())

res_plankton$debug_interp%>%mutate(year=year(target_date))%>%filter(year%in%c(2011:2015))%>%filter(name=="Chaetognatha")

res_hplc$df_long_rf%>%mutate(year=year(target_date))%>%filter(year%in%c(2011:2023))%>%filter(var=="tchla")%>%select(target_date,var,season,obs,lin_value,final)


# Save them in interm files

# For the environment
# Ctd
final_panel_ctd <- res_ctd$final_panel_rf
clim_wk <- res_ctd$clim_wk_smooth
write_tsv(final_panel_ctd, file="data/interm_files/final_panel_ctd.tsv")
write_tsv(clim_wk, file="data/interm_files/clim_wk_ctd.tsv")


#Std 
final_panel_std <- res_std$final_panel_rf
clim_wk_std <- res_std$clim_wk_smooth
write_tsv(final_panel_std, file="data/interm_files/final_panel_std.tsv")
write_tsv(clim_wk_std, file="data/interm_files/clim_wk_std.tsv")
