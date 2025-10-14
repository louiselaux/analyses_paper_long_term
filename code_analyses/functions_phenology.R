
signif_stars <- function(x) {
  dplyr::case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}


# Function to detrend time series

compute_detrended <- function(dstl, stats) {
  dstl %>%
    dplyr::left_join(stats %>% dplyr::select(name, gls_intercept, gls_slope),
                     by = "name") %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(
      trend_gls     = gls_intercept + gls_slope * as.numeric(target_date),
      detrended_raw = value - trend_gls
    ) %>%
    dplyr::ungroup()
}


# Function to perform daily interpolation

interp_daily <- function(df_long,
                         group_col = name,
                         date_col  = target_date,
                         value_col = detrended_raw) {
  gsym <- rlang::ensym(group_col)
  dsym <- rlang::ensym(date_col)
  vsym <- rlang::ensym(value_col)
  
  df_long %>%
    dplyr::filter(!is.na(!!dsym)) %>%
    dplyr::group_by(!!gsym) %>%
    dplyr::group_modify(~{
      x <- .x %>%
        dplyr::filter(!is.na(!!vsym)) %>%
        dplyr::arrange(!!dsym)
      
      if (nrow(x) < 2) {
        return(tibble(date = as.Date(character()),
                      !!gsym := .x[[rlang::as_name(gsym)]][1],
                      value = numeric()))
      }
      
      dmin <- min(x[[rlang::as_name(dsym)]])
      dmax <- max(x[[rlang::as_name(dsym)]])
      grid <- tibble(date = seq(dmin, dmax, by = "1 day"))
      
      y_interp <- castr::interpolate(
        x    = x[[rlang::as_name(dsym)]],
        y    = x[[rlang::as_name(vsym)]],
        xout = grid$date
      )
      
      tibble(date = grid$date,
             !!gsym := .x[[rlang::as_name(gsym)]][1],
             value = y_interp)
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(year = lubridate::year(date),
                  doy  = lubridate::yday(date))
}


# Function to smooth 

prep_days_above <- function(interp_daily_df,
                            k_smooth = 20,
                            thr_fun = function(v) stats::quantile(v, 0.75, na.rm = TRUE)) {
  interp_daily_df %>%
    dplyr::group_by(name) %>%
    dplyr::arrange(date, .by_group = TRUE) %>%
    dplyr::mutate(
      value_smooth = if (k_smooth > 0)
        zoo::rollmean(value, k = k_smooth, fill = NA, align = "center") else value,
      thresh_rel   = thr_fun(value_smooth)
    ) %>%
    dplyr::ungroup()
}


# Function to count the number of days above a threshold
count_days_above_by_year <- function(prep_df) {
  prep_df %>%
    dplyr::group_by(name, year) %>%
    dplyr::summarise(
      days_above = sum(value_smooth > thresh_rel, na.rm = TRUE),
      .groups = "drop"
    )
}


# Function to have a statistical model (glm) to look for a significant effect of time of nb of days above a given threshold

fit_days_glm <- function(days_tbl) {
  days_tbl %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~{
      x <- .x %>% dplyr::filter(is.finite(days_above), is.finite(year))
      if (nrow(x) < 5 || length(unique(x$days_above)) < 2) return(tibble())
      m <- glm(days_above ~ year, family = quasipoisson(link = "log"), data = x)
      s <- summary(m)
      tibble(
        term     = "year",
        estimate = unname(coef(m)["year"]),
        p.value  = s$coefficients["year","Pr(>|t|)"],
        stars    = signif_stars(s$coefficients["year","Pr(>|t|)"]),
        n        = nrow(x)
      )
    }) %>%
    dplyr::ungroup()
}


# COG for plankton only

compute_cog <- function(interp_daily_df, min_total = 0) {
  interp_daily_df %>%
    dplyr::group_by(name, year) %>%
    dplyr::summarise(
      total = sum(value, na.rm = TRUE),
      cog   = ifelse(total > 0, sum(doy * value, na.rm = TRUE) / total, NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.na(total) | total >= min_total)
}

fit_cog_lm <- function(cog_tbl) {
  cog_tbl %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~{
      x <- .x %>% dplyr::filter(is.finite(cog), is.finite(year))
      if (nrow(x) < 5) return(tibble())
      m <- lm(cog ~ year, data = x)
      sm <- summary(m)
      tibble(
        term     = "year",
        estimate = unname(coef(m)["year"]),
        p.value  = sm$coefficients["year","Pr(>|t|)"],
        stars    = signif_stars(sm$coefficients["year","Pr(>|t|)"]),
        n        = nrow(x)
      )
    }) %>%
    dplyr::ungroup()
}


# Plots 

plot_days_models <- function(days_tbl) {
  ggplot(days_tbl, aes(x = year, y = days_above)) +
    geom_point(alpha = 0.7) +
    geom_smooth(
      method = "glm",
      method.args = list(family = quasipoisson(link = "log")),
      se = TRUE, color = "red"
    ) +
    facet_wrap(~ name, scales = "free_y") +
    theme_bw() +
    labs(x = "Année", y = "Jours > seuil (Q75)")
}

plot_cog_models <- function(cog_tbl) {
  ggplot(cog_tbl, aes(year, cog)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    facet_wrap(~ name, scales = "free_y") +
    theme_bw() +
    labs(x = "Année", y = "COG (jour de l’année)")
}
