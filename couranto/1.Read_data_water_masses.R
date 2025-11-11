# Load the library to read .nc files  
library(ncdf4)
library(maps)
library(mapdata)

# Where is the file
nc_file <- "couranto/med-aroundB.nc"  

nc <- nc_open(nc_file)

# Step 0 : extract what we want #
u     <- ncvar_get(nc, "uo")      # dims: lat x depth x time
v     <- ncvar_get(nc, "vo")
lat   <- ncvar_get(nc, "latitude")
lon   <- ncvar_get(nc, "longitude")
depth <- ncvar_get(nc, "depth")
time  <- ncvar_get(nc, "time")

origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
dates  <- as.Date(origin + time)

dim(u)

# Step 1 : Take the good latitude #
valid_by_lat_u <- apply(u, 1, function(x) any(is.finite(x)))
valid_by_lat_u
# [1] TRUE FALSE  -> only first one has data

ilat    <- which(valid_by_lat_u)[1]   # 1
lat_sel <- lat[ilat]

# Extract this latitude-> matrix [depth x time]
u_prof <- u[ilat, , ]   # depth x time
v_prof <- v[ilat, , ]

dim(u_prof)

# Step 2 : Build the dataframe with all the depths #
ntime  <- ncol(u_prof)
ndepth <- nrow(u_prof)

# grid
df_curr <- expand.grid(
  depth = depth,             # longueur = ndepth
  itime = seq_len(ntime)     # 1..ntime
)

df_curr$u <- as.vector(u_prof)
df_curr$v <- as.vector(v_prof)

df_curr$date  <- dates[df_curr$itime]
df_curr$speed <- sqrt(df_curr$u^2 + df_curr$v^2)

head(df_curr)
summary(df_curr$speed)


# Step 3 : Look at a plot #
curr_allz <- df_curr %>%
  mutate(
    year  = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m")),
    month_factor = factor(month, levels = 1:12,
                          labels = c("Jan","Fév","Mar","Avr","Mai","Jun",
                                     "Jul","Aoû","Sep","Oct","Nov","Déc")),
    season = case_when(
      month %in% c(12, 1, 2) ~ "Hiver",
      month %in% c(3, 4, 5)  ~ "Printemps",
      month %in% c(6, 7, 8)  ~ "Été",
      TRUE                   ~ "Automne"
    )
  )

clim_month_depth <- curr_allz %>%
  group_by(depth, month_factor) %>%
  summarise(
    speed_mean = median(speed, na.rm = TRUE),
    speed_sd   = sd(speed, na.rm = TRUE),
    u_mean     = median(u, na.rm = TRUE),
    v_mean     = median(v, na.rm = TRUE),
    .groups = "drop"
  )
ggplot(clim_month_depth,
       aes(x = month_factor, y = depth, fill = speed_mean)) +
  geom_tile() +
  scale_y_reverse() +
  xlab("Month") +
  ylab("Depth (m)") +
  labs(fill = "Mean speed \n(m/s)") +
  theme_minimal()

# And seasonal variations
clim_season_depth <- curr_allz %>%
  group_by(season, depth) %>%
  summarise(speed_mean = mean(speed, na.rm = TRUE), .groups = "drop")

ggplot(clim_season_depth,
       aes(x = speed_mean, y = depth, colour = season)) +
  geom_path(size = 1.1) + geom_point() +
  scale_y_reverse() +
  xlab("Mean speed (m/s)") +
  ylab("Depth (m)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# And now the direction
clim_season_uv <- curr_allz %>%
  group_by(season, depth) %>%
  summarise(
    u_mean = mean(u, na.rm = TRUE),
    v_mean = mean(v, na.rm = TRUE),
    speed_mean = mean(speed, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(clim_season_uv,
       aes(x = u_mean, y = depth, colour = season)) +
  geom_path() + geom_point() +
  scale_y_reverse() +
  xlab("u moyen (m/s)  (Est > 0, Ouest < 0)") +
  ylab("Profondeur (m)") +
  theme_minimal()

##### Step 4 : get the mean per date #####
head(df_curr)
curr_mean_by_date <- df_curr %>%
  group_by(date) %>%
  summarise(
    u_mean     = mean(u,     na.rm = TRUE),  
    v_mean     = mean(v,     na.rm = TRUE),
    speed_mean = mean(speed, na.rm = TRUE),
    .groups = "drop"
  )


# Take the mean on 2 weeks around the target date
dates_biweekly <- readr::read_tsv("code_CCM/ref.tsv") %>%
  dplyr::select(target_date)

dates_vec <- dates_biweekly$target_date   # date vector
currents  <- curr_mean_by_date  

currents_biweekly_around <- tibble(
  target_date = dates_vec,
  
  curr_speed_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      currents$speed_mean[
        currents$date >= (.x - 7) &
          currents$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  ),
  
  curr_u_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      currents$u_mean[
        currents$date >= (.x - 7) &
          currents$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  ),
  
  curr_v_mean_around = map_dbl(
    dates_vec,
    ~ mean(
      currents$v_mean[
        currents$date >= (.x - 7) &
          currents$date <= (.x + 7)
      ],
      na.rm = TRUE
    )
  )
)

head(currents_biweekly_around)

# Check it worked out properly
d_test <- as.Date("1992-01-15")

# Calculate around the date
check_curr_around <- curr_mean_by_date %>%
  filter(
    date >= d_test - 7,
    date <= d_test + 7
  ) %>%
  summarise(
    n_days      = n(),
    mean_speed  = mean(speed_mean, na.rm = TRUE),
    mean_u      = mean(u_mean,     na.rm = TRUE),
    mean_v      = mean(v_mean,     na.rm = TRUE),
    first_date  = min(date),
    last_date   = max(date)
  )

check_curr_around

# Compare with the data table
currents_biweekly_around %>%
  filter(target_date == d_test)

# Save it
currents_scaled <- currents_biweekly_around

currents_scaled %>%
  as.data.frame() %>%
  select(-target_date) %>%     
  scale() %>%     # Scale for CCM later             
  as.data.frame() %>%
  write_tsv("code_CCM/currents_scaled.tsv")
