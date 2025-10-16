####Trends in hplc data#####
# Load the libraries
library(readr)
library(dplyr)
library(tidyverse)

hplccomb_allpoints <- read_tsv("data/hplc_coeffs.tsv")

hplc_final<- hplccomb_allpoints

#Select surface depth for now

hplc_final$depth<- as.factor(hplc_final$depth)
hplc_final <- hplc_final %>% select(date, depth, pico_chla, nano_chla, micro_chla, tchla)

#Put "surface" for 0 and 1
#hplc_trend <- hplc_final %>%
#  mutate(depth=if_else(depth %in% c(0, 1), "srf", depth))
hplc_trend<- hplc_final
table(hplc_trend$depth)

# TO DO : investigate why one depth is 40
hplc_trend %>% filter(depth==40)
hplc_trend %>% filter(date =="2023-10-10")
# Remove the problematic depth for now
hplc_trend <- hplc_trend %>% filter(!depth=="40")

hplc_trend%>%group_by(date)%>%filter(depth=="srf")%>%summarise(n=n())%>%filter(n>1)

hplc_trend<- hplc_trend %>%filter(depth=="10")

#Change the date format
#hplc_trend <- hplc_trend %>%
#  mutate(date = as.Date(`date`, format = "%m/%d/%y"))

hplc_trend <- hplc_trend %>%arrange(date)

#Some plots
hplc_trend%>%ggplot()+geom_path(aes(x=date,y=tchla))+scale_y_log10()+theme_bw()



hplc_trend %>%
  ggplot() +
  geom_path(aes(x = date, y = tchla)) +
  scale_y_log10() +
  scale_x_date(date_breaks = "12 month", date_labels = "%b %Y") +
  theme_bw()

#Investigate a problem
hplc_trend%>%mutate(year=year(date))%>%group_by(year)%>%filter(year=="2017")%>%select(date)

#hplc_trend%>%filter("date"=='1/31/17')
hplc_trend%>%filter("date"== "2017-01-31")
# TO DO : check this compare with before

#####Regularize the time serie #####

## Regularize series

# Define a 2 weeks regular date grid
ref <- tibble(
  target_date=seq(from=as.Date("2012-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)


# identify years in which the number of obs is larger than usual
pbs <- ref %>%
  count(year) %>%
  filter(n>26)
# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# match data based on these reference dates
avail_dates <- unique(hplc_trend$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates),
    date_diff = abs(closest_date - target_date) |> as.numeric()
  )

# insert the data based on the matched dates
hplct <- left_join(ref, hplc_trend, by=c("closest_date"="date"), relationship="many-to-many")


#Pivot it
hplct <- hplct %>% #select(-depth)%>%
  pivot_longer(
    cols = pico_chla:tchla,
    names_to = "name",
    values_to = "value"
  )
# erase data for matches that are too far from the target

hplct <- hplct %>%
  mutate(value = if_else(date_diff > 6, NA, value))

#Graph
hplct%>%filter(name %in% c("micro_chla","nano_chla","pico_chla","diatoms","dinoflagellates","tchla"))%>%ggplot() + facet_wrap(~name, scales="free_y") +
  geom_point(aes(x=target_date, y=value), size=0.2) + theme_bw()+labs(x="date", y="valeur") +  theme(
    strip.background = element_rect(fill = "lightblue", color = "darkblue"),
    strip.text = element_text(color = "black", face = "bold")
  )+scale_y_log10()

hplct%>%mutate(month=factor(month(target_date)))%>%filter(name=="diatoms")%>%group_by(month,year)%>%summarise(count=n()) # 2 per month so ok


#Some plots
hplctt<- hplct
hplctt$week<- factor(week(hplctt$target_date))

#Seasonal variations of diatoms
hplctt%>%filter(name=="diatoms")%>%ggplot(aes(x = year, y = week, fill = value)) +
  geom_tile(color = "white") +
  
  labs(x = "Année", y = "Semaine", fill = "Valeur") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(
    colors = c("lightblue", "red")
  )+
  #scale_y_discrete(labels = c("Janvier", "Février", "Mars", "Avril", "Mai", "Juin",
  #                           "Juillet", "Août", "Septembre", "Octobre", "Novembre", "Décembre"))+
  scale_x_continuous(breaks = seq(min(hplctt$year), max(hplctt$year), by = 5))+
  facet_wrap(~name, scales = "free_y")

#Seasonal variations of dinoflagellates
hplctt%>%filter(name=="dinoflagellates")%>%ggplot(aes(x = year, y = week, fill = value)) +
  geom_tile(color = "white") +
  
  labs(x = "Année", y = "Semaine", fill = "Valeur") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(
    colors = c("lightblue", "red"),
    limits = c(0.0, 0.5)
  )+
  #scale_y_discrete(labels = c("Janvier", "Février", "Mars", "Avril", "Mai", "Juin",
  #                           "Juillet", "Août", "Septembre", "Octobre", "Novembre", "Décembre"))+
  scale_x_continuous(breaks = seq(min(hplctt$year), max(hplctt$year), by = 5))+
  facet_wrap(~name, scales = "free_y")

hplctt%>%mutate(month=month(target_date))%>%filter(name=="dinoflagellates")%>%ggplot()+geom_point(aes(x=as.factor(month),y=log10(as.numeric(value)), color=as.factor(year)))+ theme_bw()

#####Distribution of the different variables#####
variables <- unique(hplct$name)

ggplot(hplct, aes(x = value)) +
  geom_histogram(binwidth = 0.04, fill = "blue", color = "white", alpha = 0.7) +
  facet_wrap(~name, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Variables distribution",
    x = "Value",
    y = "Frequency"
  )

#Log transformation before
hplct_log <- hplct %>%
  mutate(value_log = log10(value))

# See hist
ggplot(hplct_log, aes(x = value_log)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "white", alpha = 0.7) +
  facet_wrap(~name, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Distribution logarithmique des variables dans `name`",
    x = "Valeur (log-transformée)",
    y = "Fréquence"
  )

## STL decomposition ----

library("stlplus")

# apply on all variables at all depth
dstl_hplc <- hplct_log |> filter(name %in% c("nano_chla","micro_chla","pico_chla","tchla"))%>%
  group_by(name) |>
  group_modify(.f=function(x,y) {
    # message(y)
    x <- arrange(x, target_date)
    x$value_log_interp <- castr::interpolate(x = x$target_date, y = x$value_log, xout = x$target_date)
    
    # if all is missing, do not do anything
    if ( all(is.na(x$value_log)) ) {
      out <- NULL
      # else perform stl
    } else {
      dec <- stlplus(x$value_log_interp, x$target_date, n.p=26, s.window="periodic", t.window=8*26)
      out <- dec$data |> select(raw, seasonal, trend, remainder)
    }
    out <- bind_cols(select(x, target_date:date_diff), out)
  }) |>
  ungroup() |>
  
  # cut the part before the variable becomes available for the first time
  group_by(name) |>
  group_modify(.f=function(x, y) {
    if (all(is.na(x$raw))) {
      out <- NULL
    } else {
      x <- arrange(x, target_date)
      start_idx <- min(which(!is.na(x$raw)))
      out <- x[start_idx:nrow(x),]
    }
    return(out)
  }) |>
  ungroup()

# plot the result
dstl_hplc |>
  pivot_longer(raw:remainder, names_to="component") |>
  mutate(component=factor(component, levels=c("raw", "trend", "seasonal", "remainder"))) |>
  ggplot() + facet_wrap(~interaction(name, component), scale="free", nrow=4) +
  geom_path(aes(x=target_date, y=value)) + theme_bw()
# -> seems to work quite well



## GLS regression ----

library("broom")
library("nlme")

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

# compute all trends
stats_hplc <- dstl_hplc |>filter(name %in% c("diatoms","dinoflagellates","micro_chla","nano_chla","pico_chla","prokaryotes","green_algae","tchla"))%>%
  mutate(deseason = trend+remainder) |>
  group_by(name) |>
  group_modify(.f=function(x, y) {
    # message(y)
    if (all(is.na(x$raw))) {
      return(data.frame())
    }
    # 0. fill missing values through linear interpolation
    x$deseason_filled <- castr::interpolate(x=x$target_date, y=x$deseason, xout=x$target_date)
    # TODO check length of "holes"
    
    #Where are the NA in deseason
    #is_na <- is.na(x$deseason)
    
    #Length
    #na_runs <- rle(is_na)
    #na_lengths <- na_runs$lengths[na_runs$values == TRUE]
    #x$na_gap_length <- ifelse(is_na, na_lengths[1], 0)
    
    # return(x)
    # 1. Mann-Kendall trend test
    mkt <- trend::mk.test(x$deseason_filled)
    
    # 2. GLS regression
    # simple model
    m <- gls(deseason_filled ~ target_date, data=x)
    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    # if autocorrelation is too strong
    if (abs(a$acf[1]) > 0.2) {
      # add AR1 model on residuals
      m <- gls(deseason_filled ~ target_date, data=x, cor=corAR1(round(a$acf[1], 1)))
      a <- pacf(residuals(m, type="normalized"), plot=FALSE)
      
      
    }
    
    # extract diagnostic information for both approaches
    bind_cols(
      glance(mkt) |> select(mk_p.value=p.value),
      glance.gls(m) |> select(r.squared, p.value=p.value, intercept, slope, cor.struct, acf=acf1) |> rename_with(function(n) {str_c("gls_", n)})
    )
  }) |>
  ungroup()

# add date range
start_stop <- dstl_hplc |>
  group_by(name) |>
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
stats_hplc <- stats_hplc |>
  mutate(
    gls_acf = abs(gls_acf),
    mk_signif = signif_stars(mk_p.value),
    gls_signif = signif_stars(gls_p.value)
  ) |>
  left_join(start_stop)

#Graphs

dstl_hplc<- dstl_hplc%>%mutate(deseason=trend+remainder)

#The plot
ggplot() +
  facet_wrap(name~., scales="free_y", ncol=1) +
  geom_path(aes(target_date, deseason), data=dstl_hplc %>% filter(name %in% c("micro_chla","nano_chla","pico_chla")), colour="grey20") +
  geom_abline(aes(slope=gls_slope, intercept=gls_intercept), data=stats_hplc%>%filter(name %in% c("micro_chla","nano_chla","pico_chla") & gls_signif %in% c("*", "**", "***")), colour="red", size=0.75, alpha=0.7) +
  #geom_abline(aes(slope=slope.gls, intercept=intercept.gls), data=subset(statsp2, signif.gls %in% c("*", "**", "***")), colour="pink", size=0.75, alpha=0.7) + theme(axis.title.y=element_blank()) + #when plotting 2nd line for recent years
  xlab("Date") +
  ylab("Deseasonalized component")+
  ggtitle("") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "lightblue")
  )




####################
#Other plots
ggplot(hplc_trend, aes(x = `Sampling date`, y = fct_rev(as.factor(`Depth (m)`)), fill = dinoflagellates)) +
  geom_raster() +
  
  scale_fill_viridis_c(option = "magma", direction = -1)+
  theme_minimal()


