##### Plankton analyses #####

# Load the libraries

library(readr)
library(purrr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(trend)      
library(strucchange) 
library(stringr)
library(stlplus)
library(nlme)
library(broom)
library(ggstats)
library(vegan)
library(worrms)
library(tidyverse)


# load the file
zoo <- read_tsv("data/zoo_abundances.tsv.gz")


# Remove useless stuffs
zz <-zoo |> filter(!str_detect(label,"bubble"),
                   !str_detect(label,"artefact"),
                   !str_detect(label,"badfocus<artefact"),
                   !str_detect(label,"detritus"),
                   !str_detect(label,"^part"))


# Remove problematic date
zz <- zz %>% filter(!object_date=="1971-05-19")

# compute community matrix: date x taxon
comm <- zz |>
  group_by(object_date, label) |> summarise(conc=sum(concentration, na.rm=TRUE)) |>
  pivot_wider(names_from = label, values_from = conc, values_fill = 0) |>
  as.data.frame()
head(comm)

# List of taxo columns
taxa <- setdiff(names(comm), "object_date")

taxa_tbl <- tibble(raw_label = taxa) %>%
  mutate(
    has_angle = str_detect(raw_label, "<"),
    stage     = if_else(has_angle, str_replace(raw_label, "<.*", ""), NA_character_),
    base_taxon= if_else(has_angle, str_replace(raw_label, ".*<", ""), raw_label),
    base_taxon= str_replace_all(base_taxon, ">", "") |> str_trim()
  )

# Some mappings
manual_map <- tribble(
  ~raw_label,                  ~base_taxon_target,
  "nauplii<Crustacea",         "Crustacea",
  "larvae<Crustacea",          "Crustacea",
  "egg<Actinopterygii",        "Actinopterygii",
  "calyptopsis<Euphausiacea",  "Euphausiacea",
  "zoea<Brachyura",            "Brachyura",
  "tail<Appendicularia",       "Appendicularia",
  "trunk<Appendicularia",      "Appendicularia",
  "head<Chaetognatha",         "Chaetognatha",
  "multiple<Copepoda",         "Copepoda",
  "multiple<other",            NA_character_,
  "other<living",              NA_character_
)

taxa_tbl <- taxa_tbl %>%
  left_join(manual_map, by = "raw_label") %>%
  mutate(base_taxon = coalesce(base_taxon_target, base_taxon)) %>%
  select(-base_taxon_target)


# Small test
wm_records_name("Calanoida", marine_only = FALSE)
wm_classification(1100)


# Use worms
worms_lookup_one <- function(name) {
  recs <- tryCatch(wm_records_name(name, marine_only = FALSE), error = function(e) NULL)
  if (is.null(recs) || nrow(recs) == 0) return(NULL)
  acc <- dplyr::filter(recs, status == "accepted")
  if (nrow(acc) > 0) acc[1, ] else recs[1, ]
}

worms_classif_ranks <- function(aphia_id, ranks =  c("genus","family","order","superorder",
                                                     "class","superclass","subclass",
                                                     "subphylum","phylum","kingdom")) {
  if (is.na(aphia_id)) return(tibble())
  cl <- tryCatch(wm_classification(aphia_id), error = function(e) NULL)
  if (is.null(cl) || nrow(cl) == 0) return(tibble())
  cl <- cl %>% mutate(rank = tolower(rank))
  out <- map_chr(ranks, ~{
    hit <- cl %>% filter(rank == .x) %>% slice_tail(n = 1)
    if (nrow(hit) == 0) NA_character_ else hit$scientificname[1]
  })
  tibble(!!!setNames(as.list(out), ranks))
}

# 
cache_path <- "data/worms_taxo_lookup.rds"

if (file.exists(cache_path)) {
  taxo_lookup <- readRDS(cache_path)
} else {
  looked <- taxa_tbl %>%
    mutate(.rec = map(base_taxon, ~{
      if (is.na(.x)) return(NULL)
      Sys.sleep(0.12) 
      worms_lookup_one(.x)
    })) %>%
    mutate(
      aphia_id      = map_int(.rec, ~ if (is.null(.x)) NA_integer_ else .x$AphiaID[1]),
      status        = map_chr(.rec, ~ if (is.null(.x)) NA_character_ else .x$status[1]),
      valid_AphiaID = map_int(.rec, ~ if (is.null(.x)) NA_integer_ else .x$valid_AphiaID[1]),
      accepted_name = map_chr(.rec, ~ if (is.null(.x)) NA_character_ else {
        vn <- .x$valid_name[1]
        if (!is.na(vn) && nzchar(vn)) vn else .x$scientificname[1]
      })
    ) %>%
    mutate(aphia_use = if_else(!is.na(valid_AphiaID) & valid_AphiaID > 0, valid_AphiaID, aphia_id)) %>%
    select(raw_label, stage, base_taxon, accepted_name, aphia_use, status)
  
  ranks_needed <- c("genus","family","order","superorder",
                    "class","superclass","subclass",
                    "subphylum","phylum","kingdom")
  ranks_tbl <- looked %>%
    mutate(.ranks = map(aphia_use, ~{
      Sys.sleep(0.08)
      worms_classif_ranks(.x, ranks = ranks_needed)
    })) %>%
    unnest_wider(.ranks)
  
  taxo_lookup <- ranks_tbl
  
  saveRDS(taxo_lookup, cache_path)
}

# 
taxo_lookup <- taxo_lookup %>%
  relocate(raw_label, stage, base_taxon, accepted_name, status, aphia_use) %>%
  arrange(raw_label)


print(taxo_lookup, n = 50)

##### Add a functional group column #####
gelatinouspred_list <- c("Aglaura","Cnidaria","Diphyidae","Abylidae","Hydrozoa","Rhopalonema velatum")
gelatinousfil_list  <- c("Oikopleuridae","Fritillariidae","Doliolida","Salpida")

carniv_cops <- c("Candaciidae","Corycaeidae","Oncaeidae")
omni_cops   <- c("Centropagidae","Metridinidae")
herbi_cops  <- c("Calanoida","Acartiidae","Oithonidae","Temoridae","Calanidae","Harpacticoida")

taxo_lookup_fg <- taxo_lookup %>%
  mutate(
    functional_group = case_when(
      raw_label %in% gelatinouspred_list ~ "gelatinous_pred",
      raw_label %in% gelatinousfil_list  ~ "gelatinous_filt",
      raw_label %in% carniv_cops         ~ "carnivorous_cops",
      raw_label %in% omni_cops           ~ "omnivorous_cops",
      raw_label %in% herbi_cops          ~ "herbivorous_cops",
      TRUE ~ NA_character_
    )
  )

head(taxo_lookup_fg)

taxo_lookup_fg%>%select(raw_label, base_taxon, genus, family, order, class, phylum, functional_group)

# Change it on google sheet and redownload it then
#write_tsv(taxo_lookup_fg, file="data/taxo_lookup_fg.tsv")
tax <- read_tsv(file="data/taxo_gathered.csv.tsv")
head(tax)

##### Merge my data table with concentrations with this one #####

tablo <- merge(zz, tax, by.x="label", by.y="raw_label", all.x=TRUE)
dim(zz)
dim(tablo) # ok, seems fine


write_tsv(tablo, file="data/tablo.tsv")


