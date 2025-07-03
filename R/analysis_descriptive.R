#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Plant species diversity analysis
# Descriptive analysis ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke
# 2025

# Aim get an genral overview of the data


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(rstatix)


### Start ###
rm(list = ls())



## 1 - load data -------------------------------------------------------------------


### abundance data #### 
# used for diversity calculation of vegetation plots (A1-A4)
# in skript _prepare_data_3b_site_environment.R (aggregated species level)
load(file = here("R", "objects", "abundances_tot_agg.Rdata"))
load(file = here("R", "objects", "abundances_target_agg.Rdata"))


# abundance data total transect
abundances_t <- read_csv(
  here("outputs", "temp", "vegetation", "abundances_t_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) 


### site environment data ####
sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) %>%   
  select(id.site, site.type, rest.meth, land.use.hist, rest.age, region,
         hydrology, mngm.type) %>% 
  distinct()


### diversity data ####

diversity <- read_csv(
  here("data", "processed",
       "data_processed_plants_site_diversity_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Descriptive Analysis ######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 - All species -------------------------------------------------------------


data_tot <- abundances_tot_agg %>% 
  left_join(sites, by = "id.site")

# number of taxonomic entities 
# (identified up to genus, aggregates, sections, species, subspecies)
# in total on vegetation plots (A1-A4)
(no_spec_tot_A1A4 <- data_tot %>% 
  summarise(n = n_distinct(name.plant.agg)) %>% 
  pull(n)
)
# 572 taxonomic entities ("species")


# number of taxonomic entities 
# (identified up to genus, aggregates, sections, species, subspecies)
# in total on entire transect (A1-A4, T)
abundances_t %>% 
    summarise(n = n_distinct(name.plant.agg))
# 693 taxonomic entities ("species")



# use for exclusion of genus level:

# # don't count taxa at genus or family level to decrease of pseudoturnover
# # (e.g. Acer is most likely one of the other three found species and not a new species)
# levels(as.factor(data$TaxonRank))
# data <- data %>% 
#   filter(!TaxonRank %in% c("GAT"))
# 
# # number of species (identified up to species level)
# # in total on vegetation plots (A1-A4)
# data %>% 
#   summarise(n = n_distinct(name.plant.agg))
# # 544 species

# number of species per site type on vegetation plots (A1-A4)
data_tot %>%
  group_by(site.type) %>%
  summarise(entries_count = n_distinct(name.plant.agg)) %>% 
  mutate(perc_of_total = entries_count/no_spec_tot_A1A4)
# site.type entries_count perc_of_total
# 1 negative            206         0.360
# 2 positive            334         0.584
# 3 restored            517         0.904
# --> numbers are misleading because we have 4x more restored sites 
# than positive and negative reference sites!

# number of species per region on vegetation plots (A1-A4)
data_tot %>%
  group_by(region) %>%
  summarise(entries_count = n_distinct(name.plant.agg)) %>% 
  mutate(perc_of_total = entries_count/no_spec_tot_A1A4)
# region    entries_count perc_of_total
# 1 centre           402         0.703
# 2 north            238         0.416
# 3 south            343         0.600



## 2 - Target species ----------------------------------------------------------


data_tgt <- abundances_target_agg %>% 
  left_join(sites, by = "id.site")

# number of taxonomic entities 
# (identified up to genus, aggregates, sections, species, supspecies)
# in total on vegetation plots (A1-A4)
(no_spec_tgt_A1A4 <- data_tgt %>% 
  summarise(n = n_distinct(name.plant.agg)) %>% 
  pull(n)
)
# 335 taxonomic entities ("species")


# number of taxonomic entities 
# (identified up to genus, aggregates, sections, species, subspecies)
# in total on entire transect (A1-A4, T)
abundances_t %>% 
  filter(target.species == 1) %>% 
  summarise(n = n_distinct(name.plant.agg))
# 360 taxonomic entities ("species")

# use for exclusion of genus level:

# # don't count taxa at genus or family level to decrease of pseudoturnover
# # (e.g. Acer is most likely one of the other three found species and not a new species)
# levels(as.factor(data$TaxonRank))
# data <- data %>%
#   filter(!TaxonRank %in% c("GAT"))
# 
# # number of species (identified up to species level)
# # in total on vegetation plots (A1-A4)
# data %>%
#   summarise(n = n_distinct(name.plant.agg))
# # 335 species


# number of species per site type on vegetation plots (A1-A4)
data_tgt %>%
  group_by(site.type) %>%
  summarise(entries_count = n_distinct(name.plant.agg)) %>% 
  mutate(perc_of_total = entries_count/no_spec_tgt_A1A4)
# site.type entries_count perc_of_total
# 1 negative            126         0.376
# 2 positive            251         0.749
# 3 restored            309         0.922
# --> numbers are misleading because we have 4x more restored sites 
# than positive and negative reference sites!

# number of species per region on vegetation plots (A1-A4)
data_tgt %>%
  group_by(region) %>%
  summarise(entries_count = n_distinct(name.plant.agg)) %>% 
  mutate(perc_of_total = entries_count/no_spec_tgt_A1A4)
# region entries_count perc_of_total
# 1 centre           249         0.743
# 2 north            156         0.466
# 3 south            220         0.657



## 3 - Diversity ---------------------------------------------------------------


diversity %>% 
  get_summary_stats(tot.hill.0, type = "full")
# 37.0 +/- 14.9 species found on average (mean +/- sd)

diversity %>% 
  get_summary_stats(target.hill.0, type = "full")
# 30.5 +/- 13.3 characteristic species found on average (mean +/- sd)



## 4 - Sites -------------------------------------------------------------------


sites %>% 
  count(site.type)
# 1 negative     33
# 2 positive     33
# 3 restored    121

sites %>% 
  count(site.type, region)
# site.type region     n
# 1 negative  centre    13
# 2 negative  north     10
# 3 negative  south     10
# 4 positive  centre    12
# 5 positive  north     11
# 6 positive  south     10
# 7 restored  centre    41
# 8 restored  north     40
# 9 restored  south     40

sites %>% 
  count(mngm.type)
# 1 both         25
# 2 grazing      30
# 3 mowing      124
# 4 none          7
# 5 NA            1
sites %>% 
  filter(is.na(mngm.type))

sites %>% 
  count(rest.meth)
# 1 cus          21
# 2 dih          40
# 3 mga          22
# 4 res          38
# 5 NA           66

sites %>% 
  count(land.use.hist)
# 1 arable land      80
# 2 grassland        41
# 3 NA               66

sites %>% 
  count(land.use.hist, region)
# land.use.hist region     n
# 1 arable land   centre    26
# 2 arable land   north     27
# 3 arable land   south     27
# 4 grassland     centre    15
# 5 grassland     north     13
# 6 grassland     south     13
# 7 NA            centre    25
# 8 NA            north     21
# 9 NA            south     20

range(sites$rest.age, na.rm = T)
# 2 - 36 years

sites %>% 
  mutate(rest.age.cat = case_when(rest.age <= 5 ~ "<5",
                                  rest.age > 5 & rest.age <= 10 ~ "6-10",
                                  rest.age > 10 ~ ">10")) %>% 
  count(rest.age.cat, site.type, region)
# rest.age.cat site.type region     n
# 1 6-10         restored  centre    11
# 2 6-10         restored  north     11
# 3 6-10         restored  south      6
# 4 <5           restored  centre    12
# 5 <5           restored  north     18
# 6 <5           restored  south      9
# 7 >10          restored  centre    15
# 8 >10          restored  north      9
# 9 >10          restored  south     17


### EUNIS, FFH, BTT Calthion ####


# # list comparison
# eunis <- traits %>% 
#   filter(r.all.diagnostic == 1)
# ffh <- traits %>% 
#   filter(ffh.lrt == 1)
# btt <- traits %>% 
#   filter(calth.btt == 1)
# 
# 
# x <- list(
#   EUNIS = eunis$name.plant,
#   FFH = ffh$name.plant,
#   BTT_CALTH = btt$name.plant
# )
# 
# traits %>% 
#   filter(target.species == 1) %>% 
#   count()
# # total 389 target species
# 
# 
# library("ggVennDiagram")
# ggVennDiagram(x) +
#   scale_fill_gradient(low="lightblue",high = "blue") +
#   # add n manually (!)
#   annotate("text", x= 10, y= 3, label = "n= 389")
# ggsave("outputs/figures/plants_species_diversity/target_species_lists_venn.jpg")
# 
# 
# # number of species in different target species lists
# traits %>%
#   summarise(across(-name.plant, ~ sum(. > 0, na.rm = TRUE))) %>%
#   pivot_longer(everything(), names_to = "system", values_to = "count") %>% 
#   print(n = 30)
