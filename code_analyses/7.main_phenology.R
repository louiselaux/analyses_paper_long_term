##### Libraries #####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(zoo)
library(castr)
library(broom)

# Functions
source("code_analyses/functions_gls.R")       
source("code_analyses/functions_phenology.R")  


# Read what we have
dstl  <- read_tsv("data/dtsl_all.tsv")        
stats <- read_tsv("data/stats_all.tsv")      



##### Step 1 : Detrend the time series #####
dstl_det <- compute_detrended(dstl, stats)


##### Step 2 : Daily linear interpolation #####
interp_det <- interp_daily(dstl_det,
                           group_col = name,
                           date_col  = target_date,
                           value_col = detrended_raw)

##### Step 3 : Compute the number of days above a threshold #####
prep_det   <- prep_days_above(interp_det, k_smooth = 20,
                              thr_fun = function(v) quantile(v, 0.75, na.rm = TRUE))
days_det   <- count_days_above_by_year(prep_det)
glm_days   <- fit_days_glm(days_det)
print(glm_days %>% arrange(p.value))

p_days <- plot_days_models(days_det)
print(p_days)

##### Step 4: For plankton, COG #####
cog_tbl  <- compute_cog(interp_det, min_total = 0)
lm_cog   <- fit_cog_lm(cog_tbl)
print(lm_cog %>% arrange(p.value))

p_cog <- plot_cog_models(cog_tbl)
print(p_cog)

##### Step 5 : Compare with raw data #####
interp_raw <- interp_daily(dstl_det,
                           group_col = name,
                           date_col  = target_date,
                           value_col = value)

prep_raw <- prep_days_above(interp_raw, k_smooth = 20,
                            thr_fun = function(v) quantile(v, 0.75, na.rm = TRUE))
days_raw <- count_days_above_by_year(prep_raw)

# Merge to see
both_days <- dplyr::bind_rows(
  days_det %>% dplyr::mutate(approach = "Detrended (Q75)"),
  days_raw %>% dplyr::mutate(approach = "Raw (Q75)")
)


# Plot
ggplot(both_days, aes(year, days_above, color = approach)) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "glm",
    method.args = list(family = quasipoisson(link = "log")),
    se = TRUE
  ) +
  facet_wrap(~ name, scales = "free_y") +
  theme_bw() + labs(x = "Année", y = "Jours > seuil")
