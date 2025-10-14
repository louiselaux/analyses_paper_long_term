# Function 1

glance.gls <- function(m) {
  s <- summary(m)
  
  # r.squared
  f <- predict(m)
  mss <- sum((f - mean(f))^2)
  rss <- sum(residuals(m)^2)
  rsq <- mss / (mss + rss)
  
  # residuals
  shap <- shapiro.test(m$residuals)
  a <- pacf(residuals(m, type="normalized"), plot=FALSE)
  
  tibble(
    r.squared = rsq,
    
    statistic = s$tTable[2, "t-value"],
    p.value = s$tTable[2, "p-value"],
    
    intercept = m$coefficients[1],
    slope = m$coefficients[2],
    
    shapiro.p.value = shap$p.value,
    cor.struct=class(m$modelStruct$corStruct)[1],
    acf1 = a$acf[1],
    acf2 = a$acf[2]
  )
}

# Function 2

signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}


# Function 3, to remove the season

build_dstl <- function(final_panel,
                       clim_wk,
                       vars = c("T","CHLA","NO3","S","O","SIOH4","MES"),
                       depths_keep = NULL) {
  fp <- final_panel
  if (!is.null(depths_keep)) fp <- fp %>% dplyr::filter(depth %in% depths_keep)
  
  dstl <- fp %>%
    dplyr::select(target_date, depth, dplyr::all_of(vars)) %>%
    tidyr::pivot_longer(dplyr::all_of(vars), names_to = "name", values_to = "value") %>%
    dplyr::mutate(week = lubridate::isoweek(target_date)) %>%
    dplyr::left_join(clim_wk, by = c("week","depth","name")) %>%
    dplyr::mutate(
      clim    = season,
      deseason = value - clim
    ) %>%
    dplyr::select(target_date, depth, name, value, week, clim, deseason)
  
  dstl
}

# Function 4, to compute the gls trends

fit_trends <- function(dstl, depth_one = "10", acf_threshold = 0.2, min_n = 5) {
  
  stats <- dstl %>%
    dplyr::filter(depth == depth_one) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~{
      x <- .x %>% dplyr::filter(!is.na(deseason))
      if (nrow(x) < min_n) return(tibble())
      
      # MK
      mkt <- trend::mk.test(x$deseason)
      
      # GLS
      m <- nlme::gls(deseason ~ target_date, data = x)
      a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
      
      if (abs(a$acf[1]) > acf_threshold) {
        m <- nlme::gls(deseason ~ target_date, data = x,
                       correlation = nlme::corAR1(value = round(a$acf[1], 1)))
        a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
        if (length(a$acf) >= 2 && abs(a$acf[2]) > acf_threshold) {
          phi <- round(a$acf[1:2], 1)
          if (sum(phi) > 0.9) { phi <- phi - 0.1 }
          m <- nlme::gls(deseason ~ target_date, data = x,
                         correlation = nlme::corARMA(value = phi, p = 2, q = 0, form = ~ target_date))
          a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
        }
      }
      
      dplyr::bind_cols(
        broom::glance(mkt) |>
          dplyr::select(mk_p.value = p.value),
        glance.gls(m) |>
          dplyr::select(r.squared, p.value, intercept, slope, cor.struct, acf = acf1) |>
          dplyr::rename_with(~paste0("gls_", .))
      )
    }) %>%
    dplyr::ungroup()
  
  # Add observation period
  start_stop <- dstl %>%
    dplyr::filter(depth == depth_one) %>%
    dplyr::mutate(year = lubridate::year(target_date)) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(start = min(year), end = max(year), .groups = "drop")
  
  stats %>%
    dplyr::mutate(
      gls_acf    = abs(gls_acf),
      mk_signif  = signif_stars(mk_p.value),
      gls_signif = signif_stars(gls_p.value)
    ) %>%
    dplyr::left_join(start_stop, by = "name")
}


# Function 5 : To compute residuals and predictions to then have the detrended version 

compute_residuals <- function(dstl, depth_one = "10", acf_threshold = 0.2, min_n = 5) {
  
  residuals <- dstl %>%
    dplyr::filter(depth == depth_one) %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~{
      x <- .x %>% dplyr::filter(!is.na(deseason))
      if (nrow(x) < min_n) {
        return(dplyr::bind_cols(
          tibble(
            mk_p.value    = NA_real_,
            gls_r.squared = NA_real_,
            gls_p.value   = NA_real_,
            gls_intercept = NA_real_,
            gls_slope     = NA_real_,
            gls_cor.struct= NA_character_,
            gls_acf       = NA_real_
          ),
          tibble(diag = list(tibble(
            target_date       = as.Date(character()),
            depth             = numeric(),
            pred_gls_deseason = numeric(),
            resid_gls         = numeric(),
            resid_norm        = numeric()
          )))
        ))
      }
      
      mkt <- trend::mk.test(x$deseason)
      
      m <- nlme::gls(deseason ~ target_date, data = x)
      a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
      
      if (abs(a$acf[1]) > acf_threshold) {
        m <- nlme::gls(deseason ~ target_date, data = x,
                       correlation = nlme::corAR1(value = round(a$acf[1], 1)))
        a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
        if (length(a$acf) >= 2 && abs(a$acf[2]) > acf_threshold) {
          phi <- round(a$acf[1:2], 1)
          if (sum(phi) > 0.9) phi <- phi - 0.1
          m <- nlme::gls(deseason ~ target_date, data = x,
                         correlation = nlme::corARMA(value = phi, p = 2, q = 0, form = ~ target_date))
        }
      }
      
      pred    <- predict(m, newdata = x)
      r_norm  <- residuals(m, type = "normalized")
      r_raw   <- x$deseason - pred
      
      diag_tbl <- tibble(
        target_date       = x$target_date,
        depth             = x$depth,
        pred_gls_deseason = as.numeric(pred),
        resid_gls         = as.numeric(r_raw),
        resid_norm        = as.numeric(r_norm)
      )
      
      dplyr::bind_cols(
        broom::glance(mkt) |> dplyr::select(mk_p.value = p.value),
        glance.gls(m) |>
          dplyr::select(r.squared, p.value, intercept, slope, cor.struct, acf = acf1) |>
          dplyr::rename_with(~paste0("gls_", .)),
        tibble(diag = list(diag_tbl))
      )
    }) %>% dplyr::ungroup()
  
  # two tables
  stats_2 <- residuals %>% dplyr::select(-diag)
  residuals_df <- residuals %>% dplyr::select(name, diag) %>% tidyr::unnest(diag)
  
  # detrended time series
  dstl2 <- dstl %>%
    dplyr::filter(depth == as.numeric(depth_one) | depth == depth_one) %>%
    dplyr::left_join(residuals_df, by = c("name","target_date","depth")) %>%
    dplyr::mutate(
      pred_gls_deseason = as.numeric(pred_gls_deseason),
      resid_gls         = as.numeric(resid_gls),
      detrended_A       = as.numeric(value) - pred_gls_deseason,
      detrended_B       = clim + resid_gls
    )
  
  list(stats_2 = stats_2, residuals_df = residuals_df, dstl2 = dstl2)
}


# Function 6 : To compute the plots of the trends

plot_trends <- function(dstl, stats, vars_show = c("CHLA","T","S","O","NO3")) {
  ggplot() +
    facet_wrap(name ~ ., scales = "free_y", ncol = 1) +
    geom_line(
      aes(x = target_date, y = deseason),
      data = dstl %>% dplyr::filter(name %in% vars_show)
    ) +
    geom_abline(
      aes(slope = gls_slope, intercept = gls_intercept),
      data = stats %>% dplyr::filter(gls_signif %in% c("*","**","***"),
                                     name %in% vars_show),
      colour = "red", size = 1.5, alpha = 0.7
    ) +
    xlab("Date") + ylab("Deseasonalized component") +
    theme_bw()
}

