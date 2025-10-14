##### Step 1 : Regularization #####

# To regularize data 
regularize_panel <- function(df_long,
                             start_date     = as.Date("1967-01-04"),
                             by_days        = 14,
                             tolerance_days = 3) {
  stopifnot(all(c("date","depth","name","value") %in% names(df_long)))
  df <- df_long %>% mutate(date = as.Date(date))
  
  ref <- tibble(
    target_date = seq(from = start_date,
                      to   = max(df$date, na.rm = TRUE) + 15,
                      by   = by_days)
  )
  avail_dates <- unique(df$date)
  
  ref %>%
    mutate(
      closest_date = castr::closest(target_date, avail_dates),
      date_diff    = as.numeric(abs(closest_date - target_date))
    ) %>%
    left_join(df, by = c("closest_date" = "date"), relationship = "many-to-many")
}

##### Step 2 : Identify gap length #####
tag_gaps <- function(table_reg){
  table_reg %>%
    arrange(name, target_date, depth) %>%
    group_by(name, depth) %>%
    mutate(
      value_obs  = if_else(date_diff <= 3, value, NA_real_),
      is_na      = is.na(value_obs),
      start_gap  = is_na & !lag(is_na, default = FALSE),
      nb_gap     = cumsum(start_gap),
      gap_id     = if_else(is_na, nb_gap, NA_integer_)
    ) %>%
    group_by(name, depth, gap_id) %>%
    mutate(size_gap = if_else(is.na(gap_id), 0L, dplyr::n())) %>%
    ungroup()
}

##### Step 3 : Linear interpolation for small gaps  #####
interp_small_gaps <- function(ts_all, max_gap = 5){
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
      value_final  = if_else(size_gap > max_gap, NA_real_, value_final)
    ) %>%
    ungroup()
}

##### Step 4 : Calculate the smoothed median #####
weekly_climatology <- function(tab_long){
  clim_wk <- tab_long %>%
    mutate(week = isoweek(target_date)) %>%
    group_by(week, depth, name) %>%
    summarise(
      n      = dplyr::n(),
      n_miss = sum(is.na(value)),
      med    = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  clim_wk %>%
    arrange(depth, name, week) %>%
    group_by(depth, name) %>%
    mutate(
      season = castr::slide(
        med, k = 1, n = 3,
        fun = weighted.mean, w = c(1,2,1), na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    select(week, depth, name, season)
}

##### Step 5 : Compute anomalies = value - median #####
compute_anomalies <- function(tab_long, clim_wk_smooth){
  tab_long %>%
    mutate(week = isoweek(target_date)) %>%
    left_join(clim_wk_smooth, by = c("week","depth","name")) %>%
    mutate(clim = season,
           anom = value - clim)
}

##### Step 6 : Imputation RF on anomalies and reconstruction #####
impute_rf_and_reconstruct <- function(tab_anom_long, clim_wide,
                                      n_cores = max(1, parallel::detectCores() %/% 2),
                                      ntree = 300, maxiter = 20){
  tab_wide <- tab_anom_long %>%
    transmute(target_date,
              var_depth = paste0(name, "_", depth),
              anom, clim) %>%
    pivot_wider(names_from = var_depth,
                values_from = c(anom, clim),
                names_sep = "_")
  
  X_anom <- tab_wide %>% select(target_date, starts_with("anom_"))
  X_mat  <- as.data.frame(X_anom %>% select(-target_date))
  
  set.seed(9)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  imp_rf <- missForest(X_mat, ntree = ntree, nodesize = c(2,5),
                       parallelize = "variables", variablewise = TRUE,
                       maxiter = maxiter)
  stopCluster(cl)
  
  imp_rf_anom <- bind_cols(target_date = X_anom$target_date,
                           as.data.frame(imp_rf$ximp))
  
  final_wide_rf <- imp_rf_anom %>%
    left_join(clim_wide, by = "target_date")
  
  for (base in gsub("^anom_", "", setdiff(names(imp_rf_anom), "target_date"))) {
    final_wide_rf[[base]] <- final_wide_rf[[paste0("anom_", base)]] +
      final_wide_rf[[paste0("clim_", base)]]
  }
  
  err_tbl <- tibble(
    var      = names(imp_rf$ximp),
    RMSE     = sqrt(imp_rf$OOBerror),
    abs_mean = colMeans(abs(imp_rf$ximp), na.rm = TRUE),
    rel_error= RMSE / abs_mean
  ) %>% arrange(rel_error)
  
  list(final_wide_rf = final_wide_rf, rf_errors = err_tbl)
}

##### Step 7 : Reconstruction long #####
reconstruct_long_with_obs <- function(final_wide_rf, tablo_merged, vars_keep = NULL){
  if (is.null(vars_keep)){
    vars_keep <- names(tablo_merged) %>%
      setdiff(c("target_date","closest_date","depth"))
  }
  
  obs_long <- tablo_merged %>%
    select(target_date, depth, all_of(vars_keep)) %>%
    pivot_longer(cols = all_of(vars_keep), names_to = "var", values_to = "obs")
  
  base_cols <- names(final_wide_rf)
  base_cols <- base_cols[ base_cols != "target_date" &
                            !grepl("^(anom_|clim_)", base_cols) &
                            grepl("_(\\d+)$", base_cols) ]
  
  final_wide_rf_base <- final_wide_rf %>%
    select(target_date, all_of(base_cols))
  
  df_long_rf <- final_wide_rf_base %>%
    pivot_longer(cols = -target_date,
                 names_to = c("var","depth"),
                 names_pattern = "^(.*)_(\\d+)$",
                 values_to = "recon") %>%
    mutate(depth = safe_as_num(depth)) %>%
    left_join(obs_long %>% mutate(depth = safe_as_num(depth)),
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
  
  final_panel_rf <- df_long_rf %>%
    select(target_date, depth, var, final) %>%
    pivot_wider(names_from = var, values_from = final) %>%
    arrange(depth, target_date)
  
  list(df_long_rf = df_long_rf, final_panel_rf = final_panel_rf)
}

##### Step 8 : Whole pipeline #####
run_to_rf <- function(df_long,
                      start_date     = as.Date("1967-01-04"),
                      by_days        = 14,
                      tolerance_days = 3,
                      max_gap_interp = 5,
                      n_cores        = max(1, parallel::detectCores() %/% 2),
                      vars_fixed_order = NULL # ex: c("T","CHLA","NO3","S","O","SIOH4","MES")
){
  # Step 1–3 : Regularisation + small gaps
  table_reg <- regularize_panel(df_long, start_date, by_days, tolerance_days)
  ts_all    <- tag_gaps(table_reg)
  interp_lin <- interp_small_gaps(ts_all, max_gap = max_gap_interp)
  
  # Step 4 
  interp_wide <- interp_lin %>%
    select(target_date, closest_date, name, value_final, depth) %>%
    pivot_wider(names_from = name, values_from = value_final, names_prefix = "lin_")
  
  tablo_merged <- table_reg %>%
    pivot_wider(names_from = name, values_from = value) %>%
    left_join(interp_wide, by = c("target_date","closest_date","depth"))
  
  # Put NA 
  lin_vars <- grep("^lin_", names(tablo_merged), value = TRUE)
  core_vars <- gsub("^lin_", "", lin_vars)
  for (v in core_vars){
    if (v %in% names(tablo_merged))
      tablo_merged[[v]] <- dplyr::coalesce(tablo_merged[[v]], tablo_merged[[paste0("lin_", v)]])
  }
  tablo_merged <- tablo_merged %>% select(-starts_with("lin_"))
  
  
  vars_all <- setdiff(names(tablo_merged), c("target_date","closest_date","depth","date_diff"))
  if (!is.null(vars_fixed_order)) vars_all <- intersect(vars_fixed_order, vars_all)
  
  tab <- tablo_merged %>%
    transmute(target_date, depth, !!!rlang::syms(vars_all))
  
  # Step 5 : Clim + anomalies
  tab_long <- tab %>%
    pivot_longer(all_of(vars_all), names_to = "name", values_to = "value") %>%
    mutate(week = isoweek(target_date))
  
  clim_wk_smooth <- weekly_climatology(tab_long)
  write_tsv(clim_wk_smooth, file = "data/clim_wk_smooth.tsv")
  
  tab_anom_long <- compute_anomalies(tab_long, clim_wk_smooth)
  
  clim_wide <- tab_anom_long %>%
    transmute(target_date,
              var_depth = paste0(name, "_", depth),
              clim) %>%
    pivot_wider(names_from = var_depth, values_from = clim, names_prefix = "clim_")
  
  # Step 6 : imputation RF and reconstruction
  imp <- impute_rf_and_reconstruct(tab_anom_long, clim_wide, n_cores = n_cores)
  recon <- reconstruct_long_with_obs(imp$final_wide_rf, tablo_merged, vars_keep = vars_all)
  
  # What is out
  list(
    df_long_rf      = recon$df_long_rf,
    final_panel_rf  = recon$final_panel_rf,
    clim_wk_smooth  = clim_wk_smooth,
    rf_errors       = imp$rf_errors,
    debug_interp    = interp_lin %>% select(target_date, closest_date, depth, name, value, value_final, size_gap)
  )
}


safe_as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}
##### A function to plot it #####
library(dplyr)
library(lubridate)
library(ggplot2)

# Imputations in red
plot_rf_series <- function(df_long_rf,
                           vars,               
                           years    = NULL,    
                           depths   = NULL) {
  df <- df_long_rf %>% filter(var %in% vars, lubridate::year(target_date) %in% years)
  
  df_red <- df %>% dplyr::filter(show_red)
  ggplot(df, aes(x = target_date)) +
    geom_line(aes(y = final)) +
    #geom_point(aes(y = obs), size = 0.9, alpha = 0.85) + # TO DO : investigate why this line does not work
    geom_point(data = df_red,
               aes(x = target_date, y = final),
               color = "red", size = 1.1, inherit.aes = FALSE) +
    facet_wrap(~depth, scales = "free_y") +
    theme_bw() +
    labs(x = "Date", y = NULL)
}

# Function to remove median and smooth
weekly_clim_smooth <- function(df_long_idx) {
  # df_long_idx: cols = target_date, depth, name, value
  df_long_idx %>%
    dplyr::mutate(week = lubridate::isoweek(target_date)) %>%
    dplyr::group_by(week, depth, name) %>%
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(depth, name, week) %>%
    dplyr::group_by(depth, name) %>%
    dplyr::mutate(
      season = castr::slide(
        med, k = 1, n = 3,
        fun = weighted.mean, w = c(1, 2, 1), na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(week, depth, name, season)
}
