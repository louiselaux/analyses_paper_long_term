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
    RMSE     = (imp_rf$OOBerror),
    abs_mean = colMeans(abs(imp_rf$ximp), na.rm = TRUE),
    rel_error= RMSE / abs_mean
  ) %>% arrange(rel_error)
  
  list(final_wide_rf = final_wide_rf, rf_errors = err_tbl)
}

##### Step 7 : Reconstruction long #####
reconstruct_long_with_obs <- function(final_wide_rf, tablo_merged, obs_wide_raw,
                                      lin_cells, big_gap_cells, vars_keep = NULL){
  if (is.null(vars_keep)){
    vars_keep <- names(tablo_merged) %>% setdiff(c("target_date","closest_date","depth"))
  }
  
  # Observations only when within tolerance (value_obs)
  obs_long <- obs_wide_raw %>%
    dplyr::select(target_date, depth, dplyr::all_of(vars_keep)) %>%
    tidyr::pivot_longer(dplyr::all_of(vars_keep), names_to = "var", values_to = "obs")
  
  base_cols <- names(final_wide_rf)
  base_cols <- base_cols[ base_cols != "target_date" &
                            !grepl("^(anom_|clim_)", base_cols) &
                            grepl("_(\\d+)$", base_cols) ]
  
  final_wide_rf_base <- final_wide_rf %>% dplyr::select(target_date, dplyr::all_of(base_cols))
  
  df_long_rf <- final_wide_rf_base %>%
    tidyr::pivot_longer(cols = -target_date,
                        names_to   = c("var","depth"),
                        names_pattern = "^(.*)_(\\d+)$",
                        values_to  = "recon") %>%
    dplyr::mutate(depth = safe_as_num(depth)) %>%
    # Join obs, linear flags/values, big-gap flags
    dplyr::left_join(obs_long %>% dplyr::mutate(depth = safe_as_num(depth)),
                     by = c("target_date","depth","var")) %>%
    dplyr::left_join(lin_cells,     by = c("target_date","depth","var")) %>%
    dplyr::left_join(big_gap_cells, by = c("target_date","depth","var")) %>%
    dplyr::group_by(var, depth) %>%
    dplyr::mutate(
      first       = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
      last        = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
      has_obs     = is.finite(first) & is.finite(last),
      in_range    = has_obs & target_date >= first & target_date <= last,
      was_missing = is.na(obs),
      
      # Flags
      show_red_lin = !is.na(lin_value),                                # small-gap linear
      show_red_rf  = is.na(lin_value) & was_missing & in_range,        # big-gap RF
      
      # Final value preference: observed > linear > RF
      final = dplyr::case_when(
        !is.na(obs)                     ~ obs,
        !is.na(lin_value)               ~ lin_value,
        was_missing & in_range          ~ recon,
        TRUE                            ~ recon   # outside obs-range: keep recon as fallback
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(show_red = show_red_lin | show_red_rf)  # if you still want a combined flag
  
  final_panel_rf <- df_long_rf %>%
    dplyr::select(target_date, depth, var, final) %>%
    tidyr::pivot_wider(names_from = var, values_from = final) %>%
    dplyr::arrange(depth, target_date)
  
  list(df_long_rf = df_long_rf, final_panel_rf = final_panel_rf)
}



##### Step 8 : Whole pipeline  for hydro #####
run_to_rf <- function(df_long,
                      start_date     = as.Date("1967-01-04"),
                      by_days        = 14,
                      tolerance_days = 3,
                      max_gap_interp = 5,
                      n_cores        = max(1, parallel::detectCores() %/% 2),
                      vars_fixed_order = NULL # ex: c("T","CHLA","NO3","S","O","SIOH4","MES")
){
  # Step 1–3 : Regularization + small gaps
  table_reg <- regularize_panel(df_long, start_date, by_days, tolerance_days)
  ts_all    <- tag_gaps(table_reg)
  obs_wide_raw <- ts_all %>%
    select(target_date, depth, name, value_obs) %>% # value _obs  = value if in the tolerance interval otherwise NA
    tidyr::pivot_wider(names_from = name, values_from = value_obs)
  interp_lin <- interp_small_gaps(ts_all, max_gap = max_gap_interp)
  
  
  # Cells filled by linear interpolation (small gaps)
  lin_cells <- interp_lin %>%
    dplyr::filter(is.na(value_obs), !is.na(value_final), size_gap > 0, size_gap <= max_gap_interp) %>%
    dplyr::transmute(target_date, depth, var = name, lin_value = value_final)
  
  # Cells that were big gaps (left NA at this stage → will be RF later)
  big_gap_cells <- interp_lin %>%
    dplyr::filter(is.na(value_obs), size_gap > max_gap_interp) %>%
    dplyr::transmute(target_date, depth, var = name, is_big_gap = TRUE)
  
  
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
  recon <- reconstruct_long_with_obs(
    imp$final_wide_rf, tablo_merged,
    obs_wide_raw = obs_wide_raw,    
    lin_cells = lin_cells,          
    big_gap_cells = big_gap_cells,  
    vars_keep = vars_all
  )
  
  # What is out
  list(
    df_long_rf      = recon$df_long_rf,
    final_panel_rf  = recon$final_panel_rf,
    clim_wk_smooth  = clim_wk_smooth,
    rf_errors       = imp$rf_errors,
    debug_interp    = interp_lin %>% select(target_date, closest_date, depth, name, date_diff, value, value_obs,value_final, size_gap)
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

# Other function that works for plankton
plot_rf_series_2 <- function(df_long_rf, vars, years = NULL, depths = NULL) {
  df <- df_long_rf %>% dplyr::filter(var %in% vars)
  
  if (!is.null(years))  df <- df %>% dplyr::filter(lubridate::year(target_date) %in% years)
  if (!is.null(depths)) df <- df %>% dplyr::filter(depth %in% depths)
  
  if (!"show_red_lin" %in% names(df)) df$show_red_lin <- FALSE
  if (!"show_red_rf"  %in% names(df)) df$show_red_rf  <- df$show_red %||% FALSE
  
  ggplot(df, aes(x = target_date)) +
    geom_line(aes(y = final)) +
  
    # geom_point(aes(y = obs), size = 0.8, alpha = 0.6) +
    geom_point(data = df %>% dplyr::filter(show_red_lin),
               aes(y = final), shape = 16, size = 1.2, color = "orange") +
    geom_point(data = df %>% dplyr::filter(show_red_rf),
               aes(y = final), shape = 16, size = 1.2, color = "red") +
    facet_wrap(~ depth, scales = "free_y") +
    theme_bw() +
    labs(x = "Date", y = NULL,
         title = paste("Série régularisée –", paste(vars, collapse = ", ")),
         subtitle = "Orange = linear interpolation, Rouge = RF (big gaps)")
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


# Functions for plankton
tag_gaps_2 <- function(table_reg, tolerance_days_obs = 3, by_days = 14){
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


# Interpolation
interp_small_gaps_2 <- function(ts_all, max_gap = 3){
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


###### Whole pipeline plankton#####

run_to_clim <- function(df_long,
                        start_date     = as.Date("1967-01-04"),
                        by_days        = 14,
                        tolerance_days = 3,
                        max_gap_interp = 3,  
                        vars_fixed_order = NULL) {
  
  stopifnot(all(c("date","depth","name","value") %in% names(df_long)))
  
  # regularization
  table_reg <- regularize_panel(df_long, start_date, by_days, tolerance_days)
  ts_all    <- tag_gaps_2(table_reg, tolerance_days_obs = tolerance_days)
  
  # Observed values
  obs_wide_raw <- ts_all %>%
    dplyr::select(target_date, depth, name, value_obs) %>%
    tidyr::pivot_wider(names_from = name, values_from = value_obs)
  
  # linear interpolation on small gaps
  interp_lin <- interp_small_gaps_2(ts_all, max_gap = max_gap_interp)
  
  # Tag cells
  lin_cells <- interp_lin %>%
    dplyr::filter(is.na(value_obs), !is.na(value_final),
                  size_gap_tol > 0, size_gap_tol <= max_gap_interp) %>%
    dplyr::transmute(target_date, depth, var = name, lin_value = value_final)
  
  big_gap_cells <- interp_lin %>%
    dplyr::filter(is.na(value_obs), size_gap_tol > max_gap_interp) %>%
    dplyr::transmute(target_date, depth, var = name, is_big_gap = TRUE)
  
  # Wide data table
  interp_wide <- interp_lin %>%
    dplyr::select(target_date, closest_date, name, value_final, depth) %>%
    tidyr::pivot_wider(names_from = name, values_from = value_final, names_prefix = "lin_")
  
  tablo_merged <- table_reg %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    dplyr::left_join(interp_wide, by = c("target_date","closest_date","depth"))
  
  lin_vars  <- grep("^lin_", names(tablo_merged), value = TRUE)
  core_vars <- gsub("^lin_", "", lin_vars)
  for (v in core_vars) {
    if (v %in% names(tablo_merged))
      tablo_merged[[v]] <- dplyr::coalesce(tablo_merged[[v]], tablo_merged[[paste0("lin_", v)]])
  }
  tablo_merged <- tablo_merged %>% dplyr::select(-dplyr::starts_with("lin_"))
  
  # Keep variables
  vars_all <- setdiff(names(tablo_merged), c("target_date","closest_date","depth","date_diff"))
  if (!is.null(vars_fixed_order)) vars_all <- intersect(vars_fixed_order, vars_all)
  
  tab <- tablo_merged %>%
    dplyr::transmute(target_date, depth, !!!rlang::syms(vars_all))
  
  # Clim
  tab_long <- tab %>%
    tidyr::pivot_longer(dplyr::all_of(vars_all), names_to = "name", values_to = "value") %>%
    dplyr::mutate(week = lubridate::isoweek(target_date))
  
  clim_wk_smooth <- weekly_climatology(tab_long)
  readr::write_tsv(clim_wk_smooth, file = "data/clim_wk_smooth.tsv")
  
  # Reconstruction of the data table
  df_long_core <- tab %>%
    tidyr::pivot_longer(dplyr::all_of(vars_all), names_to = "var", values_to = "val_tmp") %>%
    dplyr::mutate(week = lubridate::isoweek(target_date)) %>%
    dplyr::left_join(clim_wk_smooth, by = c("week","depth","var" = "name")) %>%
    dplyr::select(-week)
  
  # Reput everything
  df_long <- df_long_core %>%
    dplyr::left_join(
      obs_wide_raw %>%
        tidyr::pivot_longer(-c(target_date, depth), names_to = "var", values_to = "obs"),
      by = c("target_date","depth","var")
    ) %>%
    dplyr::left_join(lin_cells,     by = c("target_date","depth","var")) %>%
    dplyr::left_join(big_gap_cells, by = c("target_date","depth","var")) %>%
    dplyr::group_by(var, depth) %>%
    dplyr::mutate(
      first       = suppressWarnings(min(target_date[!is.na(obs)], na.rm = TRUE)),
      last        = suppressWarnings(max(target_date[!is.na(obs)], na.rm = TRUE)),
      has_obs     = is.finite(first) & is.finite(last),
      in_range    = has_obs & target_date >= first & target_date <= last,
      was_missing = is.na(obs),
      
      # Show tags
      show_red_lin  = !is.na(lin_value),                          # small gaps being interpolated
      show_red_clim = is.na(lin_value) & was_missing & in_range,  # big gaps --> clim
      show_red_rf   = show_red_clim
    ) %>%
    dplyr::ungroup()
  
  # Finale value
  df_long$final <- dplyr::case_when(
    !is.na(df_long$obs)                  ~ df_long$obs,
    !is.na(df_long$lin_value)            ~ df_long$lin_value,
    df_long$was_missing & df_long$in_range ~ df_long$season,  
    TRUE                                 ~ df_long$season     
  )
  
  # To have as previously
  final_panel <- df_long %>%
    dplyr::select(target_date, depth, var, final) %>%
    tidyr::pivot_wider(names_from = var, values_from = final) %>%
    dplyr::arrange(depth, target_date)

  debug_interp <- interp_lin %>%
    dplyr::select(target_date, closest_date, depth, name, date_diff,
                  value, value_obs, value_final, size_gap_tol)
  
  list(
    df_long_rf      = df_long %>%  
      dplyr::mutate(show_red = show_red_lin | show_red_rf),
    final_panel_rf  = final_panel,
    clim_wk_smooth  = clim_wk_smooth,
    rf_errors       = NULL,         # no RF here 
    debug_interp    = debug_interp
  )
}

