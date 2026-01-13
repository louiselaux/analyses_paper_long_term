# Check impact of regularization with log data rather than raw data

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

# Other functions
extract_wanted_level_from_tablo_raw <- function(
    tablo,
    level_col="family_like",
    value_col="concentration",
    depth=10
) {
  tablo %>%
    group_by(date = as.Date(object_date), name = .data[[level_col]], depth) %>%
    summarise(value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(name))
}

# log input
df_log <- extract_wanted_level_from_tablo(
  tablo, level_col="family_like", value_col="concentration"
) %>% 
  filter(!is.na(name))

df_log_filled <- fill_missing_dates_taxa(df_log)

# raw input
df_raw <- extract_wanted_level_from_tablo_raw(
  tablo, level_col="family_like", value_col="concentration"
)

df_raw_filled <- fill_missing_dates_taxa(df_raw)


# Regularization
res_log <- run_to_clim(df_long = df_log_filled,
                       start_date = as.Date("1967-01-04"),
                       by_days = 14,
                       tolerance_days = 3,
                       max_gap_interp = 3)

res_raw <- run_to_clim(df_long = df_raw_filled,
                       start_date = as.Date("1967-01-04"),
                       by_days = 14,
                       tolerance_days = 3,
                       max_gap_interp = 3)

panel_log <- res_log$final_panel_rf   # en log
panel_raw <- res_raw$final_panel_rf   # en raw

taxa <- setdiff(names(panel_log), c("target_date","depth"))

panel_log_bt <- panel_log %>%
  mutate(across(all_of(taxa), ~ pmax(10^.x - 1, 0)))

# Comparison
common_taxa <- intersect(
  setdiff(names(panel_raw), c("target_date","depth")),
  setdiff(names(panel_log_bt), c("target_date","depth"))
)


A <- panel_raw %>%
  arrange(depth, target_date) %>%
  select(target_date, depth, all_of(common_taxa))

B <- panel_log_bt %>%
  arrange(depth, target_date) %>%
  select(target_date, depth, all_of(common_taxa))

# check dimensions
stopifnot(nrow(A) == nrow(B))

# example: un taxon
tax <- "Chaetognatha"
cor(A[[tax]], B[[tax]], use="complete.obs")
summary(abs(log10(A[[tax]] + 1) - log10(B[[tax]] + 1)))

# Only on added dates
obs_dates <- sort(unique(as.Date(df_raw$date)))

cmp <- tibble(
  target_date = A$target_date,
  depth = A$depth,
  a = A[[tax]],
  b = B[[tax]],
  is_added = !(A$target_date %in% obs_dates)
)

# separated stats
cmp %>%
  mutate(err_log10 = abs(log10(a + 1) - log10(b + 1)),
         rel = abs(a - b) / pmax(a, 1e-6)) %>%
  group_by(is_added) %>%
  summarise(
    n = n(),
    cor = cor(a, b, use="complete.obs"),
    med_err_log10 = median(err_log10, na.rm=TRUE),
    p95_err_log10 = quantile(err_log10, 0.95, na.rm=TRUE),
    med_rel = median(rel, na.rm=TRUE),
    p95_rel = quantile(rel, 0.95, na.rm=TRUE)
  )

# Plot
ggplot(cmp, aes(target_date)) +
  geom_line(aes(y = a, colour = "raw→reg")) +
  geom_line(aes(y = b, colour = "log→reg→back"), linetype="dashed") +
  geom_point(data = cmp %>% filter(!is_added),
             aes(y = a), size=0.8, alpha=0.5)  + scale_y_log10()+
  theme_bw() + 
  labs(x="Date", y="Concentration", colour=NULL,
       title=paste(tax, "- comparaison raw vs log (back-transform)"),
       subtitle="Points = dates observées; lignes = séries régularisées")
