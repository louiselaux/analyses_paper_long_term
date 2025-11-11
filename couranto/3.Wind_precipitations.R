##### Load the libraries #####
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readr)

# Read the two files
vent <- read_delim("couranto/vent.csv.gz", delim= ";")
param <- read_delim("couranto/parametres.csv.gz", delim=";")

names(vent)
names(param)

vent%>%distinct(NOM_USUEL) %>%print(n=200)
param%>%
  distinct(NOM_USUEL) %>% print(n=200)

# Which one to take?

stations_col <- vent %>%
  filter(NOM_USUEL %in% c("NICE")) %>%
  distinct(NOM_USUEL, NUM_POSTE, LAT, LON) # They are the same apparently

stations_col

# Keep the only station that has wind and keep good QC 

data_meteo <- vent %>%
  filter(NUM_POSTE == "06088001") %>%   # Nice aéroport
  filter(
    QRR  %in% c(0, 1, 9),
    QFFM %in% c(0, 1, 9),
    QFXY %in% c(0, 1, 9),
    QDXY %in% c(0, 1, 9)  # direction              
  ) %>%
  mutate(
    date      = as.Date(as.character(AAAAMMJJ), format = "%Y%m%d"),
    rain_mm   = RR  / 10,
    wind_mean = FFM / 10,   # m/s
    #wind_gust = FXY / 10,   # m/s
    
    dir_deg = DXY %% 360,
    dir_rad = dir_deg * pi / 180,
    
    # composante zonale : >0 = vent d'ouest, <0 = vent d'est
    wind_zonal_mean = - wind_mean * sin(dir_rad),
    wind_merid_mean = - wind_mean * cos(dir_rad)
    #wind_zonal_gust = - wind_gust * sin(dir_rad)
  ) %>%
  filter(date > as.Date("1967-01-01")) %>%
  select(date, rain_mm, wind_mean,
         dir_deg, wind_zonal_mean, wind_merid_mean)


head(data_meteo)

# Little check
data_meteo %>%
  filter(dir_deg %in% c(80: 100) | dir_deg%in% c(260: 280)) %>%  # ~est et ~ouest
  select(date, wind_mean, dir_deg, wind_zonal_mean, wind_merid_mean) %>%
  head(10)


ggplot(data_meteo, aes(x = dir_deg, y = wind_merid_mean)) +
  geom_point(alpha = 0.1) +
  geom_smooth(se = FALSE) +
  theme_bw()

data_meteo <- data_meteo %>%filter(date>"1967-01-01")

hist(data_meteo$wind_mean)

data_meteo %>% mutate(year=year(date), month=month(date))%>%group_by(month)%>%summarise(wind=mean(rain_mm, na.rm=TRUE))%>%ggplot()+geom_path(aes(x=month,y=wind)) + theme_bw() 

# Save it here
data_meteo %>% as.data.frame()%>%write_tsv(file="couranto/data_meteo.tsv")


##### Link wind data with other point B data #####

dates_biweekly <- read_tsv(file="code_CCM/ref.tsv")%>%select(target_date)

head(dates_biweekly)

dates_vec <- dates_biweekly$target_date

##### Step 0 : Mean of wind around plankton dates #####
wind <- data_meteo

wind_biweekly_around <- tibble(
  target_date = dates_vec,
  
  wind_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      wind$wind_mean[
        wind$date >= (.x - 7) &   # 7 jours before
          wind$date <= (.x + 7)     # 7 jours after
      ],
      na.rm = TRUE
    )
  ),
  
  wind_zonal_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      wind$wind_zonal_mean[
        wind$date >= (.x - 7) &
          wind$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  ),
  wind_merid_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      wind$wind_merid_mean[
        wind$date >= (.x - 7) &
          wind$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  ),
  
  rain_sum_around = map_dbl(
    dates_vec,
    ~ sum(
      wind$rain_mm[
        wind$date >= (.x - 7) &
          wind$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  )
)

head(wind_biweekly_around)

# Small check
d_test <- as.Date("1992-01-15")

check_around <- wind %>%
  filter(
    date >= d_test - 7,
    date <= d_test + 7
  ) %>%
  summarise(
    n_days    = n(),
    mean_wind = mean(wind_mean, na.rm = TRUE),
    mean_zonal = mean(wind_zonal_mean, na.rm = TRUE),
    sum_rain  = sum(rain_mm,   na.rm = TRUE),
    first_date = min(date),
    last_date  = max(date)
  )

check_around

wind_biweekly_around %>%
  filter(target_date == d_test)


# Time series of zonal component to look at strong east winds
ggplot(data_meteo, aes(x = date, y = wind_zonal_mean)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  theme_bw() +
  labs(y = "Zonal component (m/s)",
       x = "Date",
       title = "Vent zonal : >0 ouest, <0 est")



data_meteo %>%
  mutate(strong_east = wind_zonal_mean <= -1) %>%
  ggplot(aes(x = date, y = wind_zonal_mean)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(color = "grey60") +
  geom_point(data = ~ filter(.x, strong_east),
             aes(x = date, y = wind_zonal_mean),
             color = "red", alpha = 0.8) +
  theme_bw() +
  labs(y = "Zonal component (m/s)",
       x = "Date",
       title = "East winds very strong in red")


# Wind rose
data_meteo %>%
  mutate(
    dir_class = cut(
      dir_deg,
      breaks = seq(0, 360, by = 45),   # 0, 45, ..., 360
      include.lowest = TRUE,
      right = FALSE                    # [0,45), [45,90), ..., [315,360]
    )
  ) %>%
  count(dir_class)



data_meteo %>%
  mutate(dir_class = cut(dir_deg,
                         breaks = seq(0, 360, by = 45),
                         include.lowest = TRUE,
                         right = FALSE)) %>%
  count(dir_class) %>%
  mutate(mid_angle = seq(22.5, 360, by = 45)) %>%
  ggplot(aes(x = mid_angle, y = n)) +
  geom_col(width = 40, fill = "skyblue3", color = "grey30") +
  coord_polar(start = -pi/8, direction = -1) +
  scale_x_continuous(breaks = seq(0, 315, by = 45),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
  theme_minimal() +
  labs(title = "Wind rose",
       y = "Number of days", x = NULL)

# Distribution of east winds per season
seuil_est <- quantile(data_meteo$wind_zonal_mean, 0.05, na.rm = TRUE) 

data_meteo %>%
  mutate(
    season = case_when(
      month(date) %in%  c(12,1,2)  ~ "Hiver",
      month(date) %in%  c(3,4,5)  ~ "Printemps",
      month(date) %in%  c(6,7,8)  ~ "Été",
      TRUE                         ~ "Automne"
    ),
    strong_east = wind_zonal_mean <= seuil_est
  ) %>%
  count(season, strong_east)


# Save
wind_scaled <- wind_biweekly_around
wind_scaled %>% as.data.frame()%>%select(-target_date)%>%scale()%>%as.data.frame()%>%write_tsv(file="code_CCM/wind_scaled.tsv")
