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


##### Tests on regularization #####

df_test <- regularize_panel(df_plankton_filled)

head(df_test)


# 
tag_gaps <- function(table_reg, tolerance_days_obs = 3, by_days = 14){
  table_reg %>%
    dplyr::arrange(name, target_date, depth) %>%
    dplyr::group_by(name, depth) %>%
    dplyr::mutate(
      # Real observations
      value_obs     = dplyr::if_else(date_diff <= tolerance_days_obs, value, NA_real_),
      is_na_obs     = is.na(value_obs),
      start_gap_obs = is_na_obs & !dplyr::lag(is_na_obs, default = FALSE),
      nb_gap_obs    = cumsum(start_gap_obs),
      gap_id_obs    = dplyr::if_else(is_na_obs, nb_gap_obs, NA_integer_),
      
      # Values for linear interpolation
      value_tol     = dplyr::if_else(date_diff <= 7, value, NA_real_),
      is_na_tol     = is.na(value_tol),
      start_gap_tol = is_na_tol & !dplyr::lag(is_na_tol, default = FALSE),
      nb_gap_tol    = cumsum(start_gap_tol),
      gap_id_tol    = dplyr::if_else(is_na_tol, nb_gap_tol, NA_integer_)
    ) %>%
    dplyr::group_by(name, depth, gap_id_obs) %>%
    dplyr::mutate(size_gap_obs = dplyr::if_else(is.na(gap_id_obs), 0L, dplyr::n())) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(name, depth, gap_id_tol) %>%
    dplyr::mutate(size_gap_tol = dplyr::if_else(is.na(gap_id_tol), 0L, dplyr::n())) %>%
    dplyr::ungroup()
}



dt<- tag_gaps(df_test)
head(dt)
dt%>%filter(size_gap_tol>5)


# Interpolation
interp_small_gaps <- function(ts_all, max_gap = 3){
  ts_all %>%
    group_by(name, depth) %>%
    mutate(
      raw_value    = value,
      value_interp = castr::interpolate(
        x    = closest_date[!is.na(raw_value)],
        y    = raw_value[!is.na(raw_value)],
        xout = target_date
      ),
      value_final  = dplyr::coalesce(value_obs, value_interp),
      value_final  = if_else(size_gap_tol >= max_gap, NA_real_, value_final)
    ) %>%
    ungroup()
}

u<- interp_small_gaps(dt)
head(u)


u%>%mutate(year=year(target_date))%>%filter(year==1981)%>%arrange(target_date)%>%filter(name=="Acartiidae")%>%print(n=50)

u%>%mutate(year=year(target_date))%>%filter(is.na(value_final))%>%select(year,target_date)%>%distinct()%>%group_by(year)%>%summarize(count=n())


# Small plot to investigate
u%>%mutate(year=year(target_date))%>%filter(year %in% c(2011))%>%filter(name%in%c("Acartiidae"))%>%ggplot() + geom_line(aes(x=target_date,y=value_final))+geom_point(aes(x=target_date,y=value_interp),color="lightgreen")+ geom_point(aes(x=closest_date,y=value_tol),color="red")

u%>%mutate(year=year(target_date))%>%filter(name%in%c("Chaetognatha"))%>%ggplot() + geom_line(aes(x=target_date,y=value_final))+geom_point(aes(x=target_date,y=value_interp),color="lightgreen")+ geom_point(aes(x=closest_date,y=value_tol),color="red")          
