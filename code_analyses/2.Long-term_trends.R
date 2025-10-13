##### Load the libraries #####

library(castr)
library(morphr)
library(stlplus)
library(readr)
library(dplyr)
library(tidyverse)
library(trend)
library(zoo)
library(missMDA)
library(broom)
library(nlme)

##### Load the data table#####

df_long <- read_tsv("data/df_long_rf_all.tsv")
final_panel <- read_tsv("data/final_panel_rf_all.tsv")
clim_wk <- read_tsv("data/clim_wk_smooth.tsv")

final_panel<- as.data.frame(final_panel)
head(final_panel)

#final_wide <- final_wide %>% select(target_date,T_10,CHLA_10,NO3_10,S_10,O_10,SIOH4_10,MES_10)

final_panel <- final_panel %>% filter(depth=="10")
table_reg_stl <- final_panel %>% pivot_longer(cols=c("T","CHLA","NO3","S","O","SIOH4","MES"), names_to="name", values_to="value")
head(table_reg_stl)

##### Step 1 : Remove the seasonality #####
tab_anom_long <- final_panel %>% mutate(week=isoweek(target_date))%>%pivot_longer(cols=c("T","CHLA","NO3","S","O","SIOH4","MES"),names_to="name", values_to="value")%>%
  left_join(clim_wk, by = c("week","depth","name")) %>%
  mutate(
    clim = season,                
    anom = value - clim) %>%
  rename(deseason=anom)

dstl <- tab_anom_long

##### Step 2 : Gls regression #####

# Function to have the calculation of statistics

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

# Calculation of all trends
stats <- dstl %>%
  filter(depth == "10") %>%
  mutate(deseason = deseason) %>% # not very necessary
  group_by(name) %>%
  group_modify(~{
    x <- .x %>% filter(!is.na(deseason)) # should not be the case because regularized previously
    
    if (nrow(x) < 5) return(tibble())
    
    # Mann-Kendall
    mkt <- trend::mk.test(x$deseason)
    
    # GLS 
    m <- gls(deseason ~ target_date, data = x)
    a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
    
    # If autocorrelation -> AR(1)
    if (abs(a$acf[1]) > 0.2) {
      m <- gls(deseason ~ target_date, data = x,
               correlation = corAR1(value = round(a$acf[1], 1)))
      # Recalculation of Partial autocorrelation on residuals of the model
      a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
      
      # If still autocorrelation --> AR2
      if (length(a$acf) >= 2 && abs(a$acf[2]) > 0.2) {
        phi <- round(a$acf[1:2], 1)
        if (sum(phi) > 0.9) { phi <- phi - 0.1 } 
        m <- gls(deseason ~ target_date, data = x,
                 correlation = corARMA(value = phi, p = 2, q = 0, form = ~ target_date))
        a <- pacf(residuals(m, type = "normalized"), plot = FALSE)
      }
    }
    
    bind_cols(
      glance(mkt) |> select(mk_p.value = p.value),
      glance.gls(m) |>
        select(r.squared, p.value, intercept, slope, cor.struct, acf = acf1) |>
        rename_with(~paste0("gls_", .))
    )
  }) %>%
  ungroup()

# display the result
dstl<- dstl%>%mutate(year=year(target_date))
# add date range
start_stop <- dstl %>%
  group_by(name) %>%
  summarise(start=min(year), end=max(year))

# add significance stars
# mk signif + gls non-signif may mean a non linear trend
signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}
stats <- stats %>%
  mutate(
    gls_acf = abs(gls_acf),
    mk_signif = signif_stars(mk_p.value),
    gls_signif = signif_stars(gls_p.value)
  ) %>%
  left_join(start_stop)


# Save it 

write_tsv(stats, file="data/stats_all.tsv")

# What model was used
stats %>%
  mutate(model_used = case_when(
    gls_cor.struct == "corAR1"  ~ "AR1",
    gls_cor.struct == "corARMA" ~ "AR2", 
    TRUE                        ~ "GLS"
  )) %>%
  select(name, gls_slope, gls_p.value, model_used, mk_p.value, starts_with("start"), starts_with("end"))


##### Step 3: Plot to have the deseasonalised component as well as the trend component in red #####

ggplot() +
  facet_wrap(name~., scales="free_y", ncol=1) +
  geom_line(aes((order=target_date), deseason), data=dstl%>%filter(name %in% c("CHLA","T","S","O","NO3"))) +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats%>%filter( gls_signif %in% c("*", "**", "***"))%>%filter(name %in% c("CHLA","T","S","O","NO3")), colour="red", size=1.5, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ylab("Deseasonalized component") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    strip.background = element_rect(fill = "lightblue")
  )
