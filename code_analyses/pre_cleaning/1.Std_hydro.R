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


# Add degradation index
std_pivot <- std_pivot %>% mutate(DEG_INDEX=PHEO/CHLA+1e-6)


# Log the chla variable, log10+1 to better deal with 0 that will be -inf otherwise
std_pivot <- std_pivot %>%mutate(
  CHLA= log10(CHLA+1)
)


std_d <- std_pivot %>%pivot_longer(cols=c("T","S","O","NO3","NO2","PO4","SIOH4","COP","MES","CHLA","Sigma","DEG_INDEX"), values_to="value", names_to="name")

# Save it
write_tsv(std_pivot, file= "data/std_pivot.tsv")
write_tsv(std, file= "data/std.tsv")
write_tsv(std_d, file="data/std_d.tsv")
