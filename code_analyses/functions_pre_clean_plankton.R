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
  df<- tablo %>% group_by(date=as.Date(object_date), name= .data[[level_col]], depth)%>%summarise(value= sum(.data[[value_col]], na.rm = TRUE),
                                                                                                  .groups = "drop")
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
