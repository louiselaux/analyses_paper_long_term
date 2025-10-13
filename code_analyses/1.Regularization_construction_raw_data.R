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
library(missForest)
library(doParallel)
library(lubridate)


# Load all the datasets
std <- read_tsv("data/std.tsv")

##### Step 1 :  Define a regular sequence of dates #####
head(std)


### Estimate the length of gaps ### 
ref <- tibble(
  target_date = seq(from = as.Date("1967-01-04"), to = max(std$date)+15, by = 14),
  year = year(target_date))

#pbs <- ref %>% count(year) %>% filter(n > 26) , not necessary here
#ref <- ref %>% filter(!(year %in% pbs$year & yday(target_date) > 360)) %>% select(-year)

avail_dates <- unique(std$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(target_date, avail_dates),
    date_diff = as.numeric(abs(closest_date - target_date)) # days 
  )

table_reg <- ref %>%
  left_join(std, by = c("closest_date" = "date"), relationship = "many-to-many")

##### Step 2 :  Identification of gap length #####

# Identification of gaps length

ts_all <- table_reg %>%
  arrange(name, target_date, depth) %>%
  group_by(name, depth) %>%
  mutate(
    value_obs  = if_else(date_diff <= 3, value, NA_real_),
    is_na      = is.na(value_obs),
    is_na_next = lead(is_na, default = FALSE),
    int        = if_else(is_na, 1L, 0L),
    start_gap  = is_na & !lag(is_na, default = FALSE),
    nb_gap     = cumsum(start_gap),
    gap_id     = if_else(is_na, nb_gap, NA_integer_)
  ) %>%
  group_by(name, gap_id, depth) %>%
  mutate(size_gap = if_else(is.na(gap_id), 0L, n())) %>%
  ungroup()

# Check it seemed to have worked as expected
ts_all%>% select(size_gap, closest_date, name, value, target_date, depth)%>% filter(size_gap %in% c(2:20)) %>% filter(name == "T")

##### Step 3 : Linear interpolation for small gaps ##### 

interp_linear <- ts_all %>%
  group_by(name, depth) %>%
  mutate(
    raw_value   = value,
    value_obs   = if_else(date_diff > 3, NA_real_, raw_value),
    value_interp = castr::interpolate(
      x    = closest_date[!is.na(raw_value)],
      y    = raw_value[!is.na(raw_value)],
      xout = target_date
    ),
    value_final = dplyr::coalesce(value_obs, value_interp),
    value_final = if_else(size_gap > 5, NA_real_, value_final)
  ) %>%
  ungroup()


# Plot to investigate

interp_linear %>% mutate(year=year(target_date))%>% filter(year==2009)%>% filter(depth =="50")%>%ggplot() + geom_point(aes(x=closest_date, y=raw_value), shape=15, colour="blue", size=2)+
  geom_point(aes(x=target_date,y=value_final), shape=12, color="red", size=2) + theme_bw() + facet_wrap ( ~name, scales = "free_y") 

# Taking into account the different depths 
interp_linear %>% mutate(year=year(target_date))%>% filter(year==2009)%>% filter(name=="T")%>%ggplot() + geom_point(aes(x=closest_date, y=raw_value), shape=15, colour="blue", size=2)+
  geom_point(aes(x=target_date,y=value_final), shape=12, color="red", size=2) + theme_bw() + facet_wrap ( ~depth, scales = "free_y") 


##### Step 4 : Prepare the data table for after  #####

interp_wide <- interp_linear %>%
  select(target_date, closest_date, name, value_final, depth) %>%
  pivot_wider(names_from = name, values_from = value_final,
              names_prefix = "lin_")

tablo_merged <- table_reg %>%
  pivot_wider(names_from = name, values_from = value) %>%
  left_join(interp_wide, by = c("target_date","closest_date","depth"))


tablo_merged <- tablo_merged %>%
  mutate(
    T    = coalesce(T, lin_T),
    CHLA = coalesce(CHLA, lin_CHLA),
    NO3  = coalesce(NO3, lin_NO3),
    S    = coalesce(S, lin_S),
    O    = coalesce(O, lin_O),
    SIOH4= coalesce (SIOH4, lin_SIOH4),
    MES= coalesce(MES, lin_MES)
  ) %>%
  select(-starts_with("lin_"))

tab <- tablo_merged %>%
  transmute(
    target_date,
    depth,
    month = lubridate::month(target_date),
    T, CHLA, NO3, S, O, SIOH4, MES
  )

##### Step 5 : Calculation of smoothed median per week = clim #####

# Median per month per variable

vars <- c("T","CHLA","NO3","S","O","SIOH4","MES")

tab_long <- tab %>%
  pivot_longer(all_of(vars), names_to = "name", values_to = "value") %>%
  mutate(week = isoweek(target_date)) 

clim_wk <- tab_long %>%
  group_by(week, depth, name) %>%
  summarise(
    n       = n(),
    n_miss  = sum(is.na(value)),
    med     = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percent_miss = n_miss / n,
    #med = if_else(percent_miss > 0.8 | n < 10, NA_real_, med) 
  )

#clim <- tab %>% group_by(month, depth) %>% summarise(across(c(T, CHLA, NO3, S, O, SIOH4, MES), ~median(., na.rm = TRUE), .names = "{.col}_clim"), .groups = "drop") # not necessary here


# Smoothing
clim_wk_smooth <- clim_wk %>%
  arrange(depth, name, week) %>%
  group_by(depth, name) %>%
  mutate(
    med_avg = castr::slide(
      med, k = 1, n = 3,
      fun = weighted.mean, w = c(1,2,1),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  select(week, depth, name, season = med_avg)

# Weekly anomaly
tab_anom_long <- tab_long %>%
  left_join(clim_wk_smooth, by = c("week","depth","name")) %>%
  mutate(
    clim = season,                
    anom = value - clim
  )

clim_wide <- tab_anom_long %>%
  transmute(target_date,
            var_depth = paste0(name, "_", depth),
            clim) %>%
  pivot_wider(names_from = var_depth, values_from = clim, names_prefix = "clim_")



#tabc <- tab %>% left_join(clim, by = c("month","depth")) 
#nrow(tabc) == nrow(tab)
#anti_join(tab, clim, by = c("month","depth")) %>% nrow()

# Anomalies= value - clim per month
#tab_anom <- tabc %>%
#  mutate( T_anom= T-T_clim,
#          CHLA_anom= CHLA- CHLA_clim,
#          NO3_anom= NO3-NO3_clim,
#          S_anom= S-S_clim,
#          O_anom= O-O_clim,
#          SIOH4_anom= SIOH4-SIOH4_clim,
#          MES_anom= MES-MES_clim)

tab_anom <- tab_anom_long %>%
  select(target_date, depth, name, anom) %>%
  pivot_wider(
    names_from = name,
    values_from = anom,
    names_glue = "{.name}_anom"
  )

tab_final <- tab %>%
  left_join(
    tab_anom_long %>%
      select(target_date, depth, name, anom) %>%
      pivot_wider(
        names_from = name, values_from = anom,
        names_glue = "{.name}_anom"
      ),
    by = c("target_date","depth")
  )

# We need to add the depth to the variables 
tab_wide <- tab_anom_long %>%
  transmute(
    target_date,
    var_depth = paste0(name, "_", depth),
    anom,
    clim
  ) %>%
  pivot_wider(
    names_from = var_depth,
    values_from = c(anom, clim),
    names_sep = "_"
  )

# One column by var depth
X_anom <- tab_wide %>%
  select(target_date, starts_with("anom_"))
X_anom_pca <- X_anom %>% select(-target_date)

##### Step 6 : Impute large holes from Random Forest #####
X_anom_pca <- as.data.frame(X_anom_pca)

set.seed(9)
cl <- makeCluster(6) 
registerDoParallel(cl) 

imp_rf <- missForest(X_anom_pca, ntree=300, nodesize=c(2,5), parallelize="variables", variablewise=TRUE, maxiter=20) 
stopCluster(cl)



error <- tibble(
  var= imp_rf$ximp |> names(),
  RMSE=sqrt(imp_rf$OOBerror),
  abs_mean=colMeans(abs(imp_rf$ximp), na.rm=TRUE),
  # NB: some values are negative => to estimate the average "size" of the variable we take the absolute value
  rel_error=RMSE/abs_mean
) |>
  arrange(rel_error) |>
  print()

## Get a data_table with all values
imp_rf_anom <- dplyr::bind_cols(
  target_date = X_anom$target_date,
  as.data.frame(imp_rf$ximp)
)

## Re-add the anomaly
final_wide_rf <- imp_rf_anom %>%
  left_join(clim_wide, by = "target_date")

for (base in gsub("^anom_", "", setdiff(names(imp_rf_anom), "target_date"))) {
  final_wide_rf[[base]] <- final_wide_rf[[paste0("anom_", base)]] +
    final_wide_rf[[paste0("clim_", base)]]
}

#write_tsv(final_wide_rf, file="data/final_wide_rf.tsv")

## Rechange the format #####
vars <- c("T","CHLA","NO3","S","O","SIOH4","MES")

obs_long <- tablo_merged %>%
  select(target_date, depth, all_of(vars)) %>%
  tidyr::pivot_longer(cols = all_of(vars), names_to = "var", values_to = "obs")

base_cols <- names(final_wide_rf)
base_cols <- base_cols[ base_cols != "target_date" &
                          !grepl("^(anom_|clim_)", base_cols) &
                          grepl("_(\\d+)$", base_cols) ]

final_wide_rf_base <- final_wide_rf %>%
  select(target_date, all_of(base_cols))

df_long_rf <- final_wide_rf_base %>%
  pivot_longer(
    cols = -target_date,
    names_to = c("var","depth"),
    names_pattern = "^(.*)_(\\d+)$", 
    values_to = "recon"
  ) %>%
  mutate(depth = as.numeric(depth)) %>%                
  left_join(obs_long %>% mutate(depth = as.numeric(depth)),
            by = c("target_date","depth","var")) %>%
  group_by(var, depth) %>%
  mutate(
    first       = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
    last        = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
    has_obs     = is.finite(first) & is.finite(last),
    in_range    = has_obs & target_date >= first & target_date <= last,
    was_missing = is.na(obs),
    final       = ifelse(was_missing & in_range, recon, obs),
    show_red    = was_missing & in_range
  ) %>%
  ungroup()

head(df_long_rf)

ggplot(df_long_rf %>% filter(var=="T"), aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(data = subset(df_long_rf, !was_missing),
             aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = subset(df_long_rf, show_red),
             aes(y = final), color = "red", size = 0.7) +
  facet_grid(var ~ depth, scales = "free_y") +
  theme_bw()

# Save the output for imputation with random forest somewhere 
head(df_long_rf)

readr::write_tsv(df_long_rf, "data/df_long_rf_all.tsv")
final_panel_rf <- df_long_rf %>%
  select(target_date, depth, var, final) %>%
  pivot_wider(names_from = var, values_from = final) %>%
  arrange(depth, target_date) 

head(final_panel_rf)
write_tsv(final_panel_rf, "data/final_panel_rf_all.tsv")


# A plot to confirm
# Focus year 1994
df_T_1994 <- df_long_rf %>%
  filter(var == "T", lubridate::year(target_date) %in% c(1993, 1994, 1995))

ggplot(df_T_1994, aes(x = target_date)) +
  geom_line(aes(y = final)) +
  geom_point(aes(y = obs), size = 1.0, alpha = 0.85) +
  geom_point(data = subset(df_T_1994, show_red),
             aes(y = final), color = "red", size = 1.2) +
  facet_wrap(~depth, scales = "free_y") +
  theme_bw()



