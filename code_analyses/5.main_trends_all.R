##### Load the libraries #####
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(trend)
library(broom)
library(nlme)
library(patchwork)

# Load the functions 
source("code_analyses/functions_gls.R")


############ Step 1: Compute trends For the environment ##### 
##### Example with the environment STD + CTD #####

# To have the environment merged together

# Load the data tables
final_panel_std <- read_tsv("output/final_panel_std.tsv")
clim_wk_std <- read_tsv("output/clim_wk_std.tsv")
long_std<- read_tsv("output/long_std.tsv")

#This was the first version, with B and only B points 
#final_panel_ctd <- read_tsv("data/interm_files/final_panel_ctd.tsv")
#clim_wk_ctd <- read_tsv("data/interm_files/clim_wk_ctd.tsv")

final_panel_ctd <- read_tsv("output/final_panel_ctd.tsv")
clim_wk_ctd <- read_tsv("output/clim_wk_ctd.tsv")
long_ctd<- read_tsv("output/long_ctd.tsv")

# For environment : get only the range of interest
cov_std <- long_std%>%
  dplyr::filter(!is.na(obs)) %>%                
  dplyr::group_by(var, depth) %>%
  dplyr::summarise(first_obs = min(target_date),
                   last_obs  = max(target_date),
                   .groups = "drop")


cov_ctd <- long_ctd %>%
  dplyr::filter(!is.na(obs)) %>%
  dplyr::group_by(var, depth) %>%
  dplyr::summarise(first_obs = min(target_date),
                   last_obs  = max(target_date),
                   .groups = "drop")



##### Step 1 : deseasonalisation #####

# To merge everything for the environment
vars_keep   <- c("T","CHLA","NO3","S","O","SIOH4","DEG_INDEX")

depth_one <- "10"

dstl_std <- build_dstl(
  final_panel = final_panel_std %>% dplyr::mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk_std,
  vars        = intersect(vars_keep, names(final_panel_std)),
  depths_keep = as.numeric(depth_one)
)

vars_keep_2   <- c("MLD","pycnocline","stratif_index","thermocline")

dstl_ctd <- build_dstl(
  final_panel = final_panel_ctd %>% dplyr::mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk_ctd,
  vars        = intersect(vars_keep_2, names(final_panel_ctd)),
  depths_keep = as.numeric(depth_one)
)

# Born it

dstl_ctd_bounded <- dstl_ctd %>%
  dplyr::rename(var = name) %>%
  dplyr::left_join(cov_ctd, by = c("var","depth")) %>%
  dplyr::filter(!is.na(first_obs),
                target_date >= first_obs,
                target_date <= last_obs) %>%
  dplyr::rename(name = var)


dstl_std_bounded <- dstl_std %>%
  dplyr::rename(var = name) %>%
  dplyr::left_join(cov_std, by = c("var","depth")) %>%
  dplyr::filter(!is.na(first_obs),
                target_date >= first_obs,
                target_date <= last_obs) %>%
  dplyr::rename(name = var)


# Save it
readr::write_tsv(dstl_std, "output/dstl_std_10m.tsv")
readr::write_tsv(dstl_ctd, "data/interm_files/dtsl_ctd.tsv")

##### Step 2 : Trends Gls/ MK #####

# To merge all the environment
stats_std <- fit_trends(dstl_std_bounded, depth_one = depth_one, acf_threshold = 0.2, min_n = 5)
stats_ctd <- fit_trends(dstl_ctd_bounded, depth_one = depth_one, acf_threshold = 0.2, min_n = 5)



##### Merge them all#####
dstl <- rbind(dstl_std_bounded,dstl_ctd_bounded)
stats <- rbind(stats_std, stats_ctd)

##### Step 3: Trend graphs #####
vars_show <- vars_keep

vars_show <- c("T","CHLA","NO3","S","O","DEG_index","stratif_index","MLD","pycnocline","thermocline","DEG_INDEX")

p_trends <- plot_trends(dstl, stats, vars_show = vars_show)
print(p_trends)

p_trends2 <- p_trends + facet_wrap(~ name, ncol = 3, nrow = 4, scales = "free_y")
print(p_trends2)

p_trends2 <- p_trends2 +
  theme(
    strip.background = element_rect(fill = "lightblue", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )
print(p_trends2)

# Change order of variables
custom_order <- c("T", "S", "CHLA", "NO3", "O", "SIOH4", 
                  "DEG_INDEX", "stratif_index", "MLD", "pycnocline", "thermocline")


dstl <- dstl %>%
  dplyr::mutate(
    name = factor(name, levels = custom_order)
  )

stats <- stats %>%
  dplyr::mutate(
    name = factor(name, levels = custom_order)
  )

vars_show <- custom_order

p_trends <- plot_trends(dstl, stats, vars_show = vars_show)
print(p_trends)

p_trends2 <- p_trends + facet_wrap(~ name, ncol = 3, nrow = 4, scales = "free_y")
p_trends2 <- p_trends2 +
  theme(
    strip.background = element_rect(fill = "lightblue", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )
print(p_trends2)


############ Step 2 : For the plankton #################################### #################################################################################
##### Need family and order level #####
res_plankton_family <- readRDS("output/res_plankton_family.rds")
#res_plankton_family <- readRDS("output/res_plankton_order.rds")
# res_plankton_family <- readRDS("output/res_plankton_functional.rds")
res_plankton_family <- res_plankton_family

res_plankton_order <- readRDS("output/res_plankton_order.rds")
#res_plankton_fg <- res_plankton

levels(as.factor(res_plankton_order$df_long_rf$var))
levels(as.factor(res_plankton_family$df_long_rf$var))

final_panel_family <- res_plankton_family$final_panel_rf
clim_wk_family     <- res_plankton_family$clim_wk_smooth
final_panel_order <- res_plankton_order$final_panel_rf
clim_wk_order     <- res_plankton_order$clim_wk_smooth



#final_panel<- res_plankton_fg$final_panel_rf
#clim_wk<- res_plankton_fg$clim_wk_smooth

# Choose the variables and the depth
vars_keep_order <- names(final_panel_order) %>%
  setdiff(c("target_date", "depth")) %>%
  (\(x) x[!str_detect(x, "^(multiple|tail|trunk|temp|seaweed|other)")] )()

vars_keep_family_cops <- c("Calanidae","Centropagidae","Candaciidae","Metridinidae","Temoridae","Acartiidae","Harpacticoida","Acartiidae","Corycaeidae","Oithonidae","Oncaeidae")

vars_keep_fg<- c("carnivorous_cops","gelatinous_filt","gelatinous_pred","herbivorous_cops","omnivorous_cops")

#vars_keep <- c("Acartiidae","Calanidae","Corycaeidae","Oithonidae","Abylidae")

#vars_keep <- c("tchla","micro_chla","nano_chla","pico_chla")

depth_one   <- "10" 

# Step 1 : Deseasonalisation
final_panel <- final_panel_family
clim_wk <- clim_wk_family

dstl <- build_dstl(
  final_panel = final_panel%>% dplyr::mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk,
  vars        = intersect(vars_keep_family_cops, names(final_panel)),
  depths_keep = as.numeric(depth_one)
)

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

#readr::write_tsv(dstl, "data/interm_files/dstl_plankton.tsv")
#readr::write_tsv(stats, "data/interm_files/stats_plankton.tsv")

##### Step 3: Trend graphs #####
vars_show <- vars_keep_family_cops

p_trends <- plot_trends(dstl, stats, vars_show = vars_show)
print(p_trends)

# Make it nicer order
p_trends_order <- p_trends + facet_wrap(~ name, ncol = 6, nrow = 6, scales = "free_y")
print(p_trends_order)

p_trends_order <- p_trends_order +
  theme(
    strip.background = element_rect(fill = "#d1ab75", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )
print(p_trends_order)

# Make it nicer family
p_trends_family <- p_trends + facet_wrap(~ name, ncol = 2, nrow = 6, scales = "free_y")
print(p_trends_family)

p_trends_family <- p_trends_family +
  theme(
    strip.background = element_rect(fill = "#d1ab75", color = NA),
    strip.text = element_text(color = "black", face = "bold")
  )
print(p_trends_family)

# Plot using patchwork with everything close to each other
(p_trends_order | p_trends_family )

(p_trends_order | p_trends_family) +
  plot_annotation(tag_levels = "A") +    
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0, 1))

pA <- p_trends_order  + ggtitle(NULL, subtitle = "Order level")
pB <- p_trends_family + ggtitle(NULL, subtitle = "Focus on Calanoida at the family level")

(pA | pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.subtitle = element_text(margin = margin(b = 6)))


#####Step 4 : Residuals and detrended time series plankton #####

# Order level
dstl_order <- build_dstl(
  final_panel = final_panel_order %>% mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk_order,
  vars        = intersect(vars_keep_order, names(final_panel_order)),
  depths_keep = as.numeric(depth_one)
)


res_out_plankton_order <- compute_residuals(dstl_order, depth_one= depth_one, acf_threshold = 0.2)

residuals_plankton_order <- res_out_plankton_order$residuals_df

# Save it
readr::write_tsv(residuals_plankton_order, "data/residuals_plankton_order.tsv")

# Family level 
dstl_family <- build_dstl(
  final_panel = final_panel_family %>% mutate(depth = as.numeric(depth)),
  clim_wk     = clim_wk_family,
  vars        = intersect(vars_keep_order, names(final_panel_family)),
  depths_keep = as.numeric(depth_one)
)


res_out_plankton_family <- compute_residuals(dstl_family, depth_one= depth_one, acf_threshold = 0.2)

residuals_plankton_family <- res_out_plankton_family$residuals_df

# Save it
readr::write_tsv(residuals_plankton_family, "data/residuals_plankton_family.tsv")
