
# Estimate the number of BP
calculate_nb_BP_try <- function(df, var = "resid_norm", date_col = "target_date", max_breaks = 5) {
  fml <- as.formula(paste(var, "~ 1"))
  bp  <- breakpoints(fml, data = df, breaks = max_breaks)
  idx <- bp$breakpoints
  if (length(idx) == 0 || all(is.na(idx))) {
    dates <- as.Date(character())
  } else {
    dates <- df[[date_col]][idx]
  }
  list(BP_results = bp, break_dates = dates)
}


    

# Plot a graph with the BP lines

plot_BP <- function(df, var, break_dates, date_col = "target_date") {
  p <- ggplot(df, aes(x = .data[[date_col]], y = .data[[var]])) +
    geom_line() +
    theme_minimal() +
    labs(x = "Date", y = var, title = paste("", var))
  
  # add vertical lines only if we have any
  if (length(break_dates) > 0) {
    p <- p + geom_vline(xintercept = as.Date(break_dates),
                        linetype = "dashed", color = "hotpink")
  }
  p
}


# Estimate the number of BP with confidence intervals
calc_bp_with_ci <- function(df, var = "resid_norm", date_col = "target_date",
                            h = 0.15, max_breaks = 5, conf_level = 0.95) {
  
  df <- df %>% dplyr::arrange(.data[[date_col]])
  fml <- as.formula(paste(var, "~ 1"))
  
  # All brzaks possible until th max number
  bp0 <- strucchange::breakpoints(fml, data = df, h = h, breaks = max_breaks)
  
  # k that minimises the BIC
  k_optimal <- strucchange::breakpoints(bp0) 
  #if k_optimal$breakpoints is NA, then k = 0
  k <- if(is.na(k_optimal$breakpoints[1])) 0 else length(k_optimal$breakpoints)
  
  # extract the model
  bp <- strucchange::breakpoints(bp0, breaks = k)
  
  idx <- bp$breakpoints
  break_dates <- if (k == 0) {
    as.Date(character())
  } else {
    as.Date(df[[date_col]][idx])
  }
  
  ci_dates <- NULL
  if (k > 0) {
    # confidence interval
    ci_obj <- confint(bp0, breaks = k, level = conf_level)
    
    # get lower and upper
    ci_idx <- ci_obj$confint[, c(1, 3), drop = FALSE]
    
    ci_dates <- data.frame(
      lower = as.Date(df[[date_col]][ci_idx[, 1]]),
      upper = as.Date(df[[date_col]][ci_idx[, 2]])
    )
  }
  
  list(bp0 = bp0, bp = bp, k = k,
       break_dates = break_dates,
       ci_dates = ci_dates)
}


# Plot with IC
plot_bp_final <- function(df, var, res, date_col = "target_date") {
  
  # graph
  p <- ggplot(df, aes(x = .data[[date_col]], y = .data[[var]])) +
    geom_line(color = "grey70") +
    theme_minimal() +
    labs(title = paste("Analyse de rupture :", var),
         subtitle = paste("Nombre de ruptures optimal (BIC) :", res$k),
         x = "Années", y = "Valeur")
  
  # if break points are found  (k > 0)
  if (res$k > 0) {
    
    # add the mean per segment
    fac <- strucchange::breakfactor(res$bp)
    # adjust a simple linear model to get the means
    m_segmented <- lm(df[[var]] ~ fac)
    df$fitted <- fitted(m_segmented)
    
    # get the means 
    p <- p + geom_line(data = df, aes(y = fitted), color = "royalblue", size = 1.2)
    
    # Get the IC
    # Have the darked areas
    p <- p + geom_rect(data = res$ci_dates, 
                       aes(xmin = lower, xmax = upper, ymin = -Inf, ymax = Inf),
                       fill = "hotpink", alpha = 0.2, inherit.aes = FALSE)
    
    # Add the breaking points
    # get the BP date
    p <- p + geom_vline(xintercept = as.numeric(res$break_dates), 
                        linetype = "dashed", color = "hotpink", size = 0.8)
  } else {
    # if k= 0; then global mean
    p <- p + geom_hline(yintercept = mean(df[[var]], na.rm = TRUE), 
                        color = "darkgreen", linetype = "dotted")
  }
  
  return(p)
}

# Test with Nile data from the package
data("Nile")

# Convert in a used data frame for the function
df_nile <- data.frame(
  target_date = as.Date(paste0(floor(time(Nile)), "-01-01")),
  flow = as.numeric(Nile)
)

res_nile <- calc_bp_with_ci(df_nile, var = "flow", date_col = "target_date", h = 0.15)

# Check number of BP found
print(paste("Nombre de ruptures trouvées :", res_nile$k))
# If k=1, date should be 1898-01-01
print(res_nile$break_dates)

# Get the plot
p_nile <- plot_bp_final(df_nile, "flow", res_nile)

# Print it
print(p_nile)

# Compare with before
res_ancienne <- calculate_nb_BP_try(df_nile, var = "flow", max_breaks = 5)
length(res_ancienne$break_dates)
res_ancienne$break_dates
