# Load the libraries
library(castr)
library(morphr)
library(stlplus)
library(readr)
library(dplyr)
library(tidyverse)
library(trend)
library(zoo)
library(missMDA)

#Load the files
std <- read_tsv("data/std_all_depths.tsv")

# Remove outliers
std <- std %>%
  mutate(value = ifelse(name == "S" & value < 37, NA, value))

# Remove this variable 
std <- std %>% filter(name != "pH")


# Pivot it here
std_pivot <- std %>%pivot_wider(names_from = "name", values_from="value" )
std_pivot<- as.data.frame(std_pivot)


# Save it
write_tsv(std_pivot, file= "data/std_pivot.tsv")
write_tsv(std, file= "data/std.tsv")
