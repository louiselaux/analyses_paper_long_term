lib_sizes <- c(100, 200, 300, 400,  500, 600, 700)

make_ccm_plot <- function(columns, target,
                          df_real, df_shuf,
                          embed_stats_tot,
                          lib_sizes = lib_sizes,
                          n = 100,
                          Tp = -1,
                          tau = -1) {
  # 1) Get the right E for each column = causee variable
  Ei <- embed_stats_tot$E_best[embed_stats_tot$variables == columns][1]
  if (is.na(Ei)) {
    stop(sprintf("Pas de E_best pour '%s' dans embed_stats_tot", columns))
  }
  
  # 2) CCM : columns xmap target
  ref_ccm <- rEDM::CCM(
    dataFrame   = df_real,
    columns     = columns,
    target      = target,
    E           = as.integer(Ei),
    Tp          = Tp,
    tau         = tau,
    libSizes    = lib_sizes,
    sample      = n,
    includeData = TRUE,
    showPlot    = FALSE,
    noTime      = TRUE
  )
  
  df_a <- ref_ccm$CCM1_PredictStat %>%
    transmute(LibSize = as.integer(LibSize),
              rho,
              type = "Real")
  
  # 3) CCM shuffled : columns_shuffled xmap target
  shuf_name <- paste0(columns, "_shuffled")
  if (!shuf_name %in% names(df_shuf)) {
    stop(sprintf("Colonne absente dans df_shuf : %s", shuf_name))
  }
  
  shuf_ccm <- rEDM::CCM(
    dataFrame   = df_shuf,
    columns     = shuf_name,
    target      = target,
    E           = as.integer(Ei),
    Tp          = Tp,
    tau         = tau,
    libSizes    = lib_sizes,
    sample      = n,
    includeData = TRUE,
    showPlot    = FALSE,
    noTime      = TRUE
  )
  
  df_b <- shuf_ccm$CCM1_PredictStat %>%
    transmute(LibSize = as.integer(LibSize),
              rho,
              type = "Shuffled")
  
  # 4) Bind the two
  df_test <- bind_rows(df_a, df_b) %>%
    mutate(type = factor(type, levels = c("Real", "Shuffled")))
  
  df_mean <- df_test %>%
    group_by(type, LibSize) %>%
    summarise(mean_rho = mean(rho), .groups = "drop")
  
  # 5) Plots : mean line
  ggplot(df_test, aes(x = LibSize, y = rho, color = type)) +
    geom_point(alpha = 0.7,
               position = position_jitter(width = 0.5, height = 0)) +
    geom_line(data = df_mean,
              aes(y = mean_rho, group = type),
              linewidth = 1) +
    scale_color_manual(values = c("Real" = "blue", "Shuffled" = "grey40")) +
    labs(
      title = paste0(columns, " xmap ", target),
      x = "LibSize",
      y = "rho (cross-map skill)",
      color = NULL
    ) +
    theme_bw(base_size = 12)
}
