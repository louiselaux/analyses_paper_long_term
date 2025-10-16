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

# Load the data set

hplc <- read_tsv("data/hplc_coeffs.tsv")

# Test to regularize
vars_hplc <- c("pico_chla","nano_chla","micro_chla","tchla")  

hplc_long_multi <- hplc %>%
  select(date, depth, all_of(vars_hplc)) %>%
  pivot_longer(all_of(vars_hplc), names_to = "name", values_to = "value") %>%
  mutate(
    date  = as.Date(date),
    depth = as.numeric(depth)
  )

hplc_reg_multi <- regularize_panel(
  df_long       = hplc_long_multi,
  start_date     = as.Date("2012-02-01"),
  by_days        = 14,
  tolerance_days = 3
)


hplc_long_multi <- hplc_long_multi%>%filter(depth !="40")
# Is this date of beginning aligned on the plankton one ?
as.integer(as.Date("2012-02-01") - origin) %% 14 == 0

# Inspect closest date to better choose the regularization
head(hplc_reg_multi)

# Plots to see the closest dates compared to target ones
xmin <- min(hplc_reg_multi$target_date, na.rm = TRUE)
xmax <- max(hplc_reg_multi$target_date, na.rm = TRUE)


tol_days <- 3

df_one <- hplc_reg_multi %>%
  filter(name == "tchla", depth == 10) %>%
  mutate(in_tol = date_diff <= tol_days)

targets <- df_one %>% distinct(target_date)
avail   <- hplc %>% filter(depth == 10) %>% distinct(date)

ggplot() +
  # 14 days 
  geom_linerange(data = targets,
                 aes(x = target_date, ymin = 0, ymax = 1),
                 color = "grey80", linewidth = 0.4) +
  # Sampling points
  geom_point(data = avail, aes(x = date, y = 0.5), size = 1) +
  # Target dates
  geom_segment(data = df_one,
               aes(x = target_date, xend = closest_date,
                   y = 0.5, yend = 0.5, color = in_tol),
               linewidth = 0.8, alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = "forestgreen", `FALSE` = "red"),
                     labels = c(`TRUE` = paste0("≤ ", tol_days, " j"),
                                `FALSE` = paste0("> ", tol_days, " j")),
                     name = "Écart à la cible") +
  scale_y_continuous(NULL, breaks = NULL) +
  theme_bw() +
  labs(x = "Date", y = NULL,
       title = "Availability of dates HPLC") +
  coord_cartesian(xlim = c(xmin, xmax))

# Other version
ggplot() +
  geom_rug(data = targets, aes(x = target_date), sides = "t", alpha = 0.4) +
  geom_rug(data = avail,   aes(x = date),        sides = "b", alpha = 0.8) +
  theme_bw() +
  labs(x = "Date", y = NULL,
       title = "Target vs available dates – tchla, 10 m") +
  coord_cartesian(xlim = c(xmin, xmax))


# What is the sampling frequency
 hplc %>% select(date)%>%distinct()%>%mutate(year=year(date))%>%group_by(year)%>%summarize(count=n())
 

 ##### See what dates are available 
 dbg <- res_hplc$debug_interp
 # Small graph to see
 res_hplc$debug_interp %>%
   filter(name %in% c("tchla"), depth == 10,
          year(target_date) == 2013) %>%
   ggplot() +
   geom_line(aes(x = target_date, y = value_final)) +                 
   geom_point(aes(x = target_date, y = value_interp), color = "lightgreen") +   
   geom_point(aes(x = closest_date, y = value_tol), color = "red") +  
   theme_bw() + labs(x = "Date", y = "Valeur (log ou brute)")

 
 # 
 res_hplc <- run_to_clim(
   df_long        = hplc_long_multi,
   start_date     = as.Date("2012-02-01"),
   by_days        = 14,
   tolerance_days = 3,
   max_gap_interp = 3
 )

 
 u<- res_hplc$debug_interp

 hplc_reg_multi <- hplc_reg_multi%>% filter(depth !="40")
 dt<- tag_gaps_2(hplc_reg_multi)
 head(dt)

 u<- interp_small_gaps_2(dt)
 head(u)
 
 u%>%mutate(year=year(target_date))%>%filter(year %in% c(2020:2023))%>%filter(name%in%c("tchla"))%>%ggplot() + geom_line(aes(x=target_date,y=value_final))+geom_point(aes(x=target_date,y=value_interp),color="lightgreen")+ geom_point(aes(x=closest_date,y=value_tol),color="red") + geom_point(aes(x=target_date, y=value_obs), color="blue")+facet_wrap(~depth)
 
##### Test with the same function as plankton
 vars_hplc <- c("pico_chla","nano_chla","micro_chla","tchla")
 
 hplc_long <- hplc %>%
  filter(depth != 40)
 
 res_hplc <- run_to_clim(
   df_long        = hplc_long,
   start_date     = as.Date("2012-02-01"), 
   by_days        = 14,
   tolerance_days = 3,
   max_gap_interp = 3
 )

 plot_rf_series_2(res_hplc$df_long_rf,    vars="nano_chla", depths=10)
 
