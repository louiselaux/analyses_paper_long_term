##### Load the libraries #####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(trend)
library(broom)
library(nlme)

# Load the functions 
source("code_analyses/functions_gls.R")

##### Load the tables needed 
final_panel <- read_tsv("data/final_panel_rf.tsv")
clim_wk     <- read_tsv("data/clim_wk_smooth.tsv") 

final_panel <- res_plankton$final_panel_rf
clim_wk<-res_plankton$clim_wk_smooth

# Choose the variables and the depth
vars_keep   <- c("T","CHLA","NO3","S","O","SIOH4","MES")
depth_one   <- "10" 

vars_keep <- c("Acartiidae","Calanidae","Corycaeidae","Oithonidae")

##### Step 1 : deseasonalisation #####
dstl <- build_dstl(
  final_panel = final_panel %>% dplyr::mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk,
  vars        = intersect(vars_keep, names(final_panel)),
  depths_keep = as.numeric(depth_one)
)

# Save it
#readr::write_tsv(dstl, "data/dtsl_all.tsv")

##### Step 2 : Trends Gls/ MK #####
stats <- fit_trends(dstl, depth_one = depth_one, acf_threshold = 0.2, min_n = 5)

# Look at what model was retained
stats %>%
  dplyr::mutate(model_used = dplyr::case_when(
    gls_cor.struct == "corAR1"  ~ "AR1",
    gls_cor.struct == "corARMA" ~ "AR2",
    TRUE                        ~ "GLS"
  )) %>%
  dplyr::select(name, gls_slope, gls_p.value, model_used, mk_p.value, start, end) %>%
  print(n = Inf)

##### Step 3: Trend graphs #####
vars_show <- vars_keep

p_trends <- plot_trends(dstl, stats, vars_show = vars_show)
print(p_trends)

###### Step 4 : Residuals and detrended time series #####
res_out <- compute_residuals(dstl, depth_one = depth_one, acf_threshold = 0.2, min_n = 5)
stats_2       <- res_out$stats_2
residuals_df  <- res_out$residuals_df
dstl2         <- res_out$dstl2


# Little check, do those two ways give the same thing?
chk <- dstl2 %>%
  dplyr::mutate(diff = abs(detrended_A - detrended_B)) %>%
  dplyr::summarise(
    n_mismatch = sum(!is.na(diff) & diff > 1e-10),
    max_diff   = max(diff, na.rm = TRUE)
  )
print(chk)

# Save this out

readr::write_tsv(stats_2,      "data/stats_with_models.tsv")
readr::write_tsv(residuals_df, "data/residuals_depth10.tsv")
readr::write_tsv(dstl2,        "data/dstl_with_detrended.tsv")

