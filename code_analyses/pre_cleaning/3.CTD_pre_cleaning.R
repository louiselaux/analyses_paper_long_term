

 # Libraries
library(castr)


# Data table
point_B <- read_tsv(file="data/point_B.tsv")

point_B%>% mutate(year=year(date))%>% select(year,date)%>%distinct()%>%group_by(year)%>%summarize(count=n())%>%print(n=50)
final_point_B_ctd<- final_point_B_ctd%>%distinct()



final_point_B_ctd<-point_B %>% group_by(date)%>%
       dplyr::summarise(
             date=date,
             DCM=maxd(fluorescence, pressure, n.smooth=2, k=3),
             MLD=mld(sigma_theta, pressure, default.depth = max(pressure)),
             thermocline = clined(temperature, pressure, n.smooth=2, k=2),
             pycnocline = clined(sigma_theta, pressure),
             stratif_index=stratif(sigma_theta, pressure, min.depths=1:10, max.depths=70:100))

final_point_B_ctd_s <- final_point_B_ctd %>%
     mutate(
         thermocline = ifelse(stratif_index < 0.3, thermocline, 80),
         pycnocline = ifelse(stratif_index < 0.3, thermocline, 80)
       )
 head(final_point_B_ctd)

write_tsv(final_point_B_ctd,file="data/final_point_B_ctd.tsv")
