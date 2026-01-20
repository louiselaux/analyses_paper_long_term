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
dstl  <- read_tsv("data/interm_files/dstl_std_bounded.tsv")        
stats <- read_tsv("data/interm_files/stats_std_bounded.tsv")      

dstl_ctd <- read_tsv("data/interm_files/dstl_ctd_bounded.tsv")
stats_ctd <- read_tsv("data/interm_files/stats_ctd_bounded.tsv")


dstl <- dstl
stats <- stats

dstl_ctd <- dstl_ctd %>% mutate(year=year(target_date))%>%filter(year!="2023")%>%select(-year)


##### 0) Combine STD and CTD #####
dstl_all  <- dplyr::bind_rows(dstl, dstl_ctd)
stats_all <- dplyr::bind_rows(stats, stats_ctd)

##### 1) Detrend #####
dstl_all_det <- compute_detrended(dstl_all, stats_all)

##### 2) Daily interpolation #####
interp_det_all <- interp_daily(
  dstl_all_det,
  group_col = name,
  date_col  = target_date,
  value_col = detrended_raw
)

interp_raw_all <- interp_daily(
  dstl_all_det,
  group_col = name,
  date_col  = target_date,
  value_col = value
)

##### 3) Days above Q75 #####
prep_det_all <- prep_days_above(
  interp_det_all, k_smooth = 20,
  thr_fun = function(v) quantile(v, 0.75, na.rm = TRUE)
)
days_det_all <- count_days_above_by_year(prep_det_all)

prep_raw_all <- prep_days_above(
  interp_raw_all, k_smooth = 20,
  thr_fun = function(v) quantile(v, 0.75, na.rm = TRUE)
)
days_raw_all <- count_days_above_by_year(prep_raw_all)

##### 4) Merge for plotting #####
both_days_all <- dplyr::bind_rows(
  days_det_all %>% dplyr::mutate(approach = "Detrended (Q75)"),
  days_raw_all %>% dplyr::mutate(approach = "Raw (Q75)")
)

##### Vars + order + labels #####
vars_plot <- c("T", "S", "CHLA", "stratif_index")

facet_order <- vars_plot
facet_labels <- c(
  T = "Temperature",
  S = "Salinity",
  CHLA = "Chlorophyll a",
  stratif_index = "Stratification index"
)

##### 5) Fit GLM on days_above for raw et detrended (to get the stars) #####
glm_det <- fit_days_glm(days_det_all) %>%
  dplyr::mutate(approach = "Detrended (Q75)")
glm_raw <- fit_days_glm(days_raw_all) %>%
  dplyr::mutate(approach = "Raw (Q75)")

glm_both <- dplyr::bind_rows(glm_det, glm_raw) %>%
  dplyr::filter(name %in% vars_plot, term == "year") %>%
  dplyr::mutate(
    stars = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

##### 6) Put the stars #####
ann <- both_days_all %>%
  dplyr::filter(name %in% vars_plot) %>%
  dplyr::group_by(name, approach) %>%
  dplyr::summarise(
    x = min(year, na.rm = TRUE) + 1,
    y = max(days_above, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(y = y * 0.95) %>%   # under the max
  dplyr::left_join(glm_both %>% dplyr::select(name, approach, stars), by = c("name","approach")) %>%
  dplyr::mutate(
    name = factor(name, levels = facet_order)
  )

##### 7) Plot #####
plot_data <- both_days_all %>%
  dplyr::filter(name %in% vars_plot) %>%
  dplyr::mutate(
    name = factor(name, levels = facet_order),
    approach = factor(approach, levels = c("Raw (Q75)", "Detrended (Q75)"))
  )

p <- ggplot(plot_data, aes(year, days_above, color = approach)) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "glm",
    method.args = list(family = quasipoisson(link = "log")),
    se = TRUE
  ) +
  facet_wrap(
    ~ name, scales = "free_y",
    labeller = labeller(name = facet_labels)
  ) +
  geom_text(
    data = ann,
    aes(x = x, y = y, label = stars, color = approach),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(x = "Year", y = "Days > Q75", color = "") +
  theme(
    strip.background = element_rect(fill = "lightblue", color = NA),
    strip.text = element_text(face = "bold")
  )

print(p)


##### Add the other facet #####
##### Load the libraries #####
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(trend)
library(FactoMineR)
library(factoextra)
library(lubridate)
library(vegan)


# Functions
source("code_analyses/functions_gls.R")       
source("code_analyses/functions_phenology.R")  


##### Load the data table at the family level, the finest taxonomical level, with time series already regularized ######

res_plankton_family <- readRDS("output/res_plankton_family.rds")

df <- res_plankton_family$final_panel_rf
head(df)

df_long <- df %>%
  filter(depth == 10) %>%         
  pivot_longer(
    cols = -c(target_date, depth),
    names_to = "name",
    values_to = "value"
  )


interp_raw <- interp_daily(df_long,
                           group_col = name,
                           date_col  = target_date,
                           value_col = value)

# COG 
cog_tbl_raw <- compute_cog(interp_raw, min_total = 0)


# Variation in time
lm_cog_raw <- fit_cog_lm(cog_tbl_raw)
lm_cog_raw %>% arrange(p.value)

# Plot
plot_cog_models(cog_tbl_raw)

# Raw data
interp_raw <- interp_daily(df_long,
                           group_col = name,
                           date_col  = target_date,
                           value_col = value)

# Linear scale
interp_raw_lin <- interp_raw %>%
  mutate(value = pmax(10^value - 1, 0))

# COG weighted by abundance
cog_tbl_raw <- compute_cog(interp_raw_lin, min_total = 0)

lm_cog_raw <- fit_cog_lm(cog_tbl_raw)
print(lm_cog_raw %>% arrange(p.value))

plot_cog_models(cog_tbl_raw)

##### Focus on some groups #####
groups_keep <- c(
  "Calanidae","Centropagidae","Temoridae","Acartiidae","Oithonidae",
  "Harpacticoida","Oncaeidae","Corycaeidae",
  "Cladocera","Chaetognatha",
  "Salpida","Doliolida","Oikopleuridae",
  "Thecofilosea","Hydrozoa","Siphonophorae","Pteropoda", "Thecofilosea",
  "Harpacticoida",
  "Pteropoda",
  "Siphonophorae",
  "Salpida",
  "Dinophyceae",
  "Annelida",
  "Cyclopoida",
  "Chaetognatha",
  "Doliolida",
  "Ctenopoda",
  "Copelata",
  "Actinopterygii",
  "Euphausiacea","Fritillaridae"
)

cog_sel <- cog_tbl_raw %>%
  dplyr::filter(name %in% groups_keep)


cog_stars <- lm_cog_raw %>%
  dplyr::filter(term == "year") %>%       
  dplyr::mutate(
    stars = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# stars
ann_cog <- cog_sel %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(
    x = min(year, na.rm = TRUE) + 1,
    y = max(cog,  na.rm = TRUE) * 0.96,
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    cog_stars %>% dplyr::select(name, stars),
    by = "name"
  )

cog_sel<- cog_sel %>%  filter(year!="2024")

p_cog <- ggplot(cog_sel, aes(x = year, y = cog)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ name, scales = "free_y") +
  geom_text(
    data = ann_cog,
    aes(x = x, y = y, label = stars),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  theme_bw() +
  labs(x = "Year", y = "COG (day-of-year)")

print(p_cog)

p_cog <- ggplot(cog_sel, aes(year, cog)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(strip.text = element_text(face="bold", size=9)) +
  labs(x="Year", y="COG (day-of-year)")

p_cog

# Figure
library(patchwork)

p_env <- p
p_plk <- p_cog

(p_env | p_plk) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.tag.position = c(0, 1)
  )

p_env + p_plk + plot_layout(widths = c(1.2, 1))


(p_env | p_cog) +
  plot_layout(widths = c(1, 2.3)) + 
  plot_annotation(tag_levels = "A")


# Adjustments
p <- p +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE)
  )

p <- p +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.2, "cm")
  )

p <- p +
  theme(
    aspect.ratio = 2,   # carré
    legend.position = "bottom left"
  )

p_cog <- p_cog + theme(legend.position = "none")



(p | p_cog) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 2.3), guides = "collect") &
  theme(legend.position = "bottom")


##### Last try
p <- p +
  theme(
    legend.position = "bottom",
    legend.box.just = "left",
    legend.justification = "left"
  )


ann2 <- ann %>% dplyr::filter(stars != "", !is.na(stars))

ann2 <- ann2 %>%
  mutate(
    x = ifelse(approach == "Raw (Q75)", x, x + 4),  # décale en x
    y = ifelse(approach == "Raw (Q75)", y, y * 0.92) # décale en y
  )


p <- p +
  geom_text(
    data = ann2,
    aes(x = x, y = y, label = stars, color = approach),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5.5,
    show.legend = FALSE
  )



p <- p +
  theme(
    aspect.ratio = 1.2,
    legend.position = "bottom",
    legend.box.just = "left",
    legend.justification = "left",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "lightblue", color = NA),
    strip.text = element_text(face = "bold")
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

# enlever légende sur B
p_cog <- p_cog + theme(legend.position = "none")

(p | p_cog) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1.1, 2.4), guides = "collect") &
  theme(legend.position = "bottom")

