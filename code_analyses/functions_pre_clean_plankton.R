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


##### Functions #####

# Define a function that will select at what level we are

extract_wanted_level_from_tablo <- function(
    tablo,
    level_col="order_like",
    value_col="concentration",
    depth=10) {
  df<- tablo %>% group_by(date=as.Date(object_date), name= .data[[level_col]], depth)%>%summarise(value= sum(.data[[value_col]], na.rm = TRUE),.groups = "drop") %>% dplyr::mutate(value = log10(value + 1))
  return(df)
}


# Function that will replace NA per 0 for plankton 
fill_missing_dates_taxa <- function(df,
                                    date_col = "date",
                                    name_col = "name",
                                    value_col = "value",
                                    depth_value = 10) {
  df <- df %>%
    dplyr::mutate(
      !!date_col := as.Date(.data[[date_col]]),
      !!name_col := as.character(.data[[name_col]])
    )
  
  all_dates <- sort(unique(df[[date_col]]))
  all_names <- sort(unique(df[[name_col]]))
  
  df_filled <- tidyr::crossing(
    !!date_col := all_dates,
    !!name_col := all_names
  ) %>%
    dplyr::left_join(df, by = c(date_col, name_col)) %>%
    dplyr::mutate(
      !!value_col := tidyr::replace_na(.data[[value_col]], 0),
      depth = depth_value
    )
  
  return(df_filled)
}


##### To deal with biovolumes
extract_biovol_from_tablo <- function(
    tablo,
    level_col = "family_like",   # family,order
    depth     = 10,
    shape     = c("ellipsoidal", "spherical")  # which one
) {
  shape <- match.arg(shape)
  
  tablo %>%
    # 1. Convert 
    dplyr::mutate(
      area_mm2  = area * (0.0106)^2,
      major_mm  = axis_major_length * 0.0106,
      minor_mm  = axis_minor_length * 0.0106,
      radius_mm = sqrt(area_mm2 / pi),
      vol_ind   = dplyr::case_when(
        shape == "spherical"  ~ (4/3) * pi * radius_mm^3,
        shape == "ellipsoidal" ~ (4/3) * pi * (major_mm / 2) * (minor_mm / 2)^2,
        TRUE ~ NA_real_
      )
    ) %>%
    # 2. aggregation date per taxon
    dplyr::group_by(
      date  = as.Date(object_date),
      name  = .data[[level_col]],
      depth = depth
    ) %>%
    dplyr::summarise(
      # nb objects with concentrations >0
      n_pos  = sum(concentration > 0, na.rm = TRUE),
      # biovolumes taken into account only if concentration >0 
      value  = dplyr::if_else(
        n_pos > 0,
        sum(vol_ind * concentration, na.rm = TRUE),
        NA_real_    # otherwise, concentration is NA
      ),
      .groups = "drop"
    )
}


