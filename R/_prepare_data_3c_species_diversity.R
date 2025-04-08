#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Plant species diversity analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(hillR)

### Start ###
rm(list = ls())


## 1 - load data -------------------------------------------------------------------

### site environment data ####
sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) %>%   
  select(id.site, site.type, region, hydrology) %>% 
  distinct()


### species abundance data site based (A1-A4) ####
abundances_site <- read_csv(
  here(
    "data", "processed", "data_processed_species_plants_site_20250306.csv"
  ),
  col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
)


### species abundance data plot based (entire transect) ####
abundances_plot <- read_csv(
  here(
    "data", "processed", "data_processed_species_plants_plot_20250306.csv"
  ),
  col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
)


### species traits data ####
traits <- read_csv(
  here("data", "processed", "data_processed_traits_plants_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) %>%
  rename_with(tolower, .cols = everything()) 



## 2 - Define target species -------------------------------------------------------


### preparation  ---------------------------------------------------------------


data <- traits %>% 
  mutate(
    # group all R-diagnostic species into one variable
    sum = rowSums(select(., r11.diagnostic:r37.diagnostic), na.rm = T ),
    r.all.diagnostic = ifelse(sum > 0, 1, 0),
    # group all FFH LRT species into one variable
    sum = rowSums(across(starts_with("ffh")), na.rm = T),
    ffh.lrt = if_else(sum > 0, 1, 0),
    # replace NA in Calthion lists
    nie.gn.gf.kennarten = replace_na(nie.gn.gf.kennarten, 0),
    sh.calthion = replace_na(sh.calthion, 0),
    `btt.bv.gn typische arten` = replace_na(`btt.bv.gn typische arten`, 0),
    `bv.§30.nasswiesen` = replace_na(`bv.§30.nasswiesen`, 0),
    bv.gruenland.feucht = replace_na(bv.gruenland.feucht, 0),
    `btt.hh.gfr kennzeichnende arten` = replace_na(`btt.hh.gfr kennzeichnende arten`, 0),
    `btt.hh.gfr wertgebende arten (rl-hh)` = replace_na(`btt.hh.gfr wertgebende arten (rl-hh)`, 0),
    lsa.nasswiesen = replace_na(lsa.nasswiesen, 0),
    # group all Calthion species into one variable
    calth.btt = if_else(
      nie.gn.gf.kennarten + sh.calthion + `btt.bv.gn typische arten`
      + `btt.hh.gfr kennzeichnende arten` + `btt.hh.gfr wertgebende arten (rl-hh)`
      + lsa.nasswiesen 
      > 0, 1, 0),
    # woody species
    wood = if_else(funct.group.x == "wood", 1, 0, missing = 0)
    )

# check for NA
data %>% 
  filter(!complete.cases(data %>% select(r.all.diagnostic, ffh.lrt, calth.btt, wood)))
# no NA


### definition  ----------------------------------------------------------------

# + EUNIS-ESy species
# + FFH-LRT species 
# + Calthion BTT 
# - woody species

data <- data %>% 
  mutate(
    target.species = if_else(
      r.all.diagnostic + ffh.lrt + calth.btt > 0, 1, 0),
    target.species = if_else(
      target.species - wood == 1 , 1, 0)
    )

# number of species in different lists
data %>%
  select(name.plant, r.all.diagnostic, ffh.lrt, calth.btt, wood,
         target.species) %>%
  summarise(across(-name.plant, ~ sum(. > 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "system", values_to = "count")

# system           count
# <chr>            <int>
# 1 r.all.diagnostic   114
# 2 ffh.lrt            349
# 3 calth.btt           76
# 4 wood                72
# 5 target.species     389


## 3 - combine tables -------------------------------------------------------------


traits <- data %>% 
  select(name.plant, target.species, r.all.diagnostic, ffh.lrt, calth.btt, rlg, funct.group.x) %>% 
  rename(funct.group = funct.group.x)

### abundances A1-A4 ----------------------------------------------------------

# add traits to abundances
abundances <- abundances_site %>% 
  left_join(traits, by = "name.plant")

# check NA
abundances %>% 
  filter(is.na(target.species))
# only unknown species or x Festulolium (doesn't have any trait): then ok

# assign target species status for x Festulolium (= 0)
abundances <- abundances %>% 
  mutate(target.species = case_when(name.plant == "x Festulolium" ~ 0,
                                    .default = target.species))
# check
abundances %>% 
  filter(name.plant == "x Festulolium")

### abundances transect ----------------------------------------------------------

# add traits to abundances
abundances_t <- abundances_plot %>% 
  left_join(traits, by = "name.plant")

# check NA
abundances_t %>% 
  filter(is.na(target.species))
# only unknown species or x Festulolium (doesn't have any trait): then ok

# assign target species status for x Festulolium (= 0)
abundances_t <- abundances_t %>% 
  mutate(target.species = case_when(name.plant == "x Festulolium" ~ 0,
                                    .default = target.species))
# check
abundances_t %>% 
  filter(name.plant == "x Festulolium")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Calculation Abundance Data (A1-A4) ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("abundances", "abundances_t", "sites", "traits")))

## 0 - Preparation -------------------------------------------------------------

#' 1. handle cf and not cf of same species on the the same site
#' --> are considered as same species
#' 2. handle family level taxa
#' --> delete to decrease pseudoturnover 
#' 3. handle genus level species
#' --> keep them if they don't occur with species of the same genus on site
#' --> delete others if only few cases or irrelevant in cover
#' 4. handle aggregates
#' 4.1 assign aggregate level for each species (new column "name.plant.agg")
#' --> species is in an aggregate: aggregate name
#' --> species is not in an aggregate: species name
#' 4.2 sum up cover of agg-level species on the same site (new column "cover.mean.agg")
#' 4.3 handle cases where target and non-target species are in the same aggregate
#' --> one species is target species --> all species of agg become target species
#' 
#' 
#' calculate diversity with aggregate level because of different 
#' considerations of agg level in regions (decrease of pseudoturnover by aggregating,
#' sensu Boch et al. 2022)
#' 


### handle cf ----
## are there cf and not cf of same species at same site?
abundances %>%
  select(id.site, name.plant, cf, cover.mean) %>% 
  group_by(id.site, name.plant) %>%
  filter(n_distinct(cf) > 1) %>%
  arrange(id.site)
# 34 entries

# combine and take sum of cover, because they are considered as the same species
abundances_1 <- abundances %>% 
  group_by(id.site, name.plant) %>%
  filter(n_distinct(cf) > 1) %>%
  mutate(cover.mean = sum(cover.mean, na.rm = TRUE)) %>%
  slice(1) %>% 
  ungroup() %>% 
  bind_rows(abundances %>%
              group_by(id.site, name.plant) %>%
              filter(n_distinct(cf) == 1)
  ) 


### handle family taxa ----

## are there taxa identified up to family level?
abundances_1 %>% 
  filter(TaxonRank == "FAM")
# two cases (M_KAL: Poaceae, M_TAN: Brassicaceae)

abundances_1 %>% 
  filter(id.site == "M_KAL") %>% 
  print(n = 50)
# other Poaceae at same site

abundances_1 %>% 
  filter(id.site == "M_TAN") %>% 
  print(n = 50)
# other Brassiceae at same site

# --> delete cases because unclear if distinct species and very little cases,
# therefore irrelevant

abundances_1a <- abundances_1 %>% 
  filter(TaxonRank != "FAM")



### handle genus taxa ----

## are there genus-level and minimum one species-level species at the same site?

# Create a dataframe with only the genus taxa
genus_data <- abundances_1a %>%
  filter(str_count(name.plant, "\\s") == 0)

# Create a dataframe with only the species taxa
species_data <- abundances_1a %>%
  filter(str_count(name.plant, "\\s") >= 1)

# Join the two dataframes on matching id.site and name.plant
matched_data <- genus_data %>%
  inner_join(species_data, by = "id.site",
             suffix = c(".genus", ".species"),
             relationship = "many-to-many") %>%
  filter(name.plant.genus == word(name.plant.species, 1)) %>%
  select(id.site, name.plant.genus, cover.mean.genus,
         name.plant.species, cover.mean.species) %>% 
  distinct(id.site, .keep_all = TRUE)
# 16 cases where minimum one species and only-genus at the same site (26 cases in total)

# delete cases 
# because unclear if same species 
# and because irrelevant (very little cases & cover < 3%)
abundances_2 <- abundances_1a %>% 
  anti_join(matched_data, by = c("id.site", "name.plant" = "name.plant.genus"))


### aggregation ----

# clean mistakes of data entry (NA in IsChildTaxonOf)
abundances_2 %>% 
  filter(is.na(IsChildTaxonOf))
abundances_3 <- abundances_2 %>% 
  mutate(IsChildTaxonOf = case_when(
    id.site == "N_UET" & name.plant == "Bromus hordeaceus" ~ "Bromus hordeaceus agg.",
    id.site == "M_ALT" & name.plant == "Stellaria apetala" ~ "Stellaria media agg.",
    id.site %in% c("M_SSB", "M_TTH") & name.plant == "Galatella linosyris" ~ "Galatella",
    .default = IsChildTaxonOf
  ))


#' assign aggregate level for each species (new column "name.plant.agg")
#' --> species is in an aggregate: aggregate name
#' --> species is not in an aggregate: species name
abundances_3 <- abundances_3 %>% 
  mutate(name.plant.agg = case_when(
    str_detect(IsChildTaxonOf, "agg.") ~ IsChildTaxonOf,
    .default = name.plant
  ))


# change names manually (aggregates not detected automatically)
# Epilobium lamyi -> tetragonum s. l.
# Festuca pulchra -> valesiaca s. l. -> ovina agg.
# Luzula multiflora s. str. -> s. l. -> campestris agg.
# Ornithogalum umbellatum is part of Ornithogalum umbellatum agg. (sensu Rothmaler)
# Vicia sativa is part of Vicia sativa agg. (sensu Rothmaler)
abundances_3 <- abundances_3 %>% 
  mutate(name.plant.agg = case_when(
    name.plant.agg == "Epilobium lamyi" ~ "Epilobium tetragonum s. l.",
    name.plant.agg == "Festuca pulchra" ~ "Festuca ovina agg.",
    name.plant.agg == "Luzula multiflora s. str." ~ "Luzula campestris agg.",
    name.plant.agg == "Vicia sativa" ~ "Vicia sativa agg.",
    name.plant.agg == "Ornithogalum umbellatum" ~ "Ornithogalum umbellatum agg.",
    .default = name.plant.agg
  ))


# check doubles on sites due to aggregation
abundances_3 %>% 
  select(id.site, name.plant, name.plant.agg, cover.mean) %>% 
  # arrange(id.site, name.plant.agg) %>% 
  arrange(name.plant) %>% 
  add_count(id.site, name.plant.agg) %>% 
  filter(n > 1) %>% 
  print(n = 74)
# several cases (74 entries)


# take sum of cover into new variable "cover.mean.agg"
abundances_3 <- abundances_3 %>% 
  group_by(id.site, name.plant.agg) %>%
  mutate(cover.mean.agg = sum(cover.mean)) %>% 
  ungroup()

# check
check <- abundances_3 %>% 
  filter(cover.mean != cover.mean.agg) %>% 
  select(id.site, name.plant, name.plant.agg, cover.mean, cover.mean.agg) %>% 
  arrange(id.site, name.plant.agg) %>% 
  print(n = 80)
# sum is correct and not more cases than doubles (n = 74)


# handle cases where target and non-target species are in the same aggregate
# filter aggregates with target and non-target species
agg_target <- abundances_3 %>%
  select(id.site, name.plant, name.plant.agg, target.species, cover.mean, cover.mean.agg) %>% 
  group_by(name.plant.agg) %>%
  filter(n_distinct(target.species) > 1) %>%
  ungroup() %>% 
  arrange(name.plant.agg) %>% 
  distinct(name.plant, name.plant.agg, .keep_all = T)
# several cases (n = 40)


# define all species of an aggregate the same target species value
# one species is target species --> all species of agg are target species
abundances_3 <- abundances_3 %>% 
  mutate(target.species.agg = case_when(
    name.plant.agg %in% agg_target$name.plant.agg ~ 1,
    .default = target.species
  ))

# check
abundances_3 %>%
  select(id.site, name.plant, name.plant.agg, target.species, target.species.agg,
         cover.mean, cover.mean.agg) %>% 
  group_by(name.plant.agg) %>%
  filter(n_distinct(target.species) > 1) %>%
  ungroup() %>% 
  distinct(name.plant, name.plant.agg, .keep_all = T) %>%
  arrange(name.plant.agg, target.species) %>% 
  print(n= 40)
   #%>%
  # write_csv(
  #   here(
  #     "outputs", "tables", "list_conflict_target_20241111.csv"
  #   ))
# ok
  # write_csv(
  #   here(
  #     "outputs", "tables", "target_plants_agg_20241111.csv"
  #   )
  # )







abundances_3 <- abundances_3 %>% 
  select(id.site, name.plant, name.plant.agg, cover.mean, cover.mean.agg,
         target.species, target.species.agg, everything())
# for diversity analysis: species at agg. level is used
# --> name.plant.agg


# list with all taxa on species level
name.plant_list_total <- abundances_3 %>% 
  distinct(name.plant, .keep_all = T) %>% 
  select(name.plant, IsChildTaxonOf, name.plant.agg, target.species, target.species.agg) %>%
  arrange(name.plant) #%>% 
  # write_csv(
  #   here(
  #     "outputs", "temp", "vegetation", "name.plant_taxa_20241211.csv"
  #   )
  # )

# list with all taxa on agg level
name.plant.agg_list_total <- abundances_3 %>% 
  distinct(name.plant.agg, .keep_all = T) %>% 
  select(name.plant.agg, target.species.agg) %>% 
  arrange(name.plant.agg) #%>% 
  # write_csv(
  #   here(
  #     "outputs", "temp", "vegetation", "name.plant.agg_taxa_20241211.csv"
  #   )
  # )




rm(list = setdiff(ls(), c("abundances_3", "abundances_t", "sites", "traits",
                          "name.plant_list_total", "name.plant.agg_list_total")))



## 1 - Hill numbers: Total diversity -------------------------------------------

# use aggregated species and aggregated species cover
# species at genus and family level are counted (they don't occur on same site with 
# taxa at species-level)

# combine where more than one species of an aggregate on a site
abundances_tot_agg <- abundances_3 %>% 
  distinct(id.site, name.plant.agg, .keep_all = T)

# save data for descriptive analysis
# save(abundances_tot_agg, file = here(
#   "R", "objects", "abundances_tot_agg.Rdata"))


## transform into species_site matrix for hillR package:
# converting long into wide data (columns = species, rows = sites)
abundances_wide <- pivot_wider(
    abundances_tot_agg,
    id_cols = id.site,
    names_from = name.plant.agg,
    values_from = cover.mean.agg,
    values_fill = 0
  )

# calculate Hill numbers with HillR package
div_data <- abundances_wide["id.site"]
# check: number of columns correct? set to number of variables of abundance table
div_data$tot.hill.0 <- hill_taxa(comm = abundances_wide[,2:573],  q = 0)
div_data$tot.hill.1 <- hill_taxa(comm = abundances_wide[,2:573],  q = 1)
div_data$tot.hill.2 <- hill_taxa(comm = abundances_wide[,2:573],  q = 2)

# add Hill numbers to site id
diversity <- left_join(sites %>% select(id.site), div_data, by = "id.site")



## 2 - Hill numbers: Target species ----------------------------------------------------------------

# keep only target species
abundances_target <- abundances_3 %>% 
  filter(target.species.agg == 1)

# combine where more than one species of an aggregate on a site
abundances_target_agg <- abundances_target %>% 
  distinct(id.site, name.plant.agg, .keep_all = T)

# save data for descriptive analysis
# save(abundances_target_agg, file = here(
#   "R", "objects", "abundances_target_agg.Rdata"))


## transform into species_site matrix for hillR package:
# converting long into wide data (columns = species, rows = sites)
abundances_target_wide <-
  pivot_wider(
    abundances_target_agg,
    id_cols = id.site,
    names_from = name.plant.agg,
    values_from = cover.mean.agg,
    values_fill = 0
  )



# calculate Hill numbers with HillR package
div_data_target <- abundances_target_wide["id.site"]
# check: number of columns correct? set to number of variables of abundance table
div_data_target$target.hill.0 <- hill_taxa(comm = abundances_target_wide[,2:336],  q = 0)
div_data_target$target.hill.1 <- hill_taxa(comm = abundances_target_wide[,2:336],  q = 1)
div_data_target$target.hill.2 <- hill_taxa(comm = abundances_target_wide[,2:336],  q = 2)


# add Hill numbers
diversity <- left_join(diversity, div_data_target, by = "id.site")



## 3 - FCSi - Index of Favourable Conservation Status --------------------------

# test if derived species = 0 in community
# if there are no derived species: set characteristic and derived species +1
diversity %>% 
  mutate(test = tot.hill.0 - target.hill.0) %>% 
  filter(test == 0)
# zero rows in tibble = ok!
# here: N_KRI and N_BEN have no derived species

data <- diversity %>% 
  # log function: default is natural logarithm. makes sense?
  mutate(fcsi.hill.0 = log(target.hill.0 + 1 / (tot.hill.0 - target.hill.0 + 1)),
         fcsi.hill.1 = log(target.hill.1 + 1 / (tot.hill.1 - target.hill.1 + 1)),
         fcsi.hill.2 = log(target.hill.2 + 1 / (tot.hill.2 - target.hill.2 + 1))
  )
diversity <- data

## 4 - Hill numbers: Red List Germany Species ------------------------------------------

# # keep only red list Germany species
# abundances_rlg <- abundances_3 %>% 
#   filter(rlg %in% c(1, 2, 3, "V")) %>% 
#   select(id.site, name.plant.agg, cover.mean) %>% 
#   # add sites again that don't contain any red list species
#   full_join(abundances_3 %>% distinct(id.site), by = "id.site")
# 
# 
# ## transform into species_site matrix for hillR package:
# # converting long into wide data (columns = species, rows = sites)
# abundances_rlg_wide <-
#   pivot_wider(
#     abundances_rlg,
#     id_cols = id.site,
#     names_from = name.plant.agg,
#     values_from = cover.mean,
#     values_fill = 0
#   )
# 
# 
# 
# # calculate Hill numbers with HillR package
# div_data_rlg <- abundances_rlg_wide["id.site"]
# # check: number of columns correct? set to number of variables of abundance table
# div_data_rlg$rlg.hill.0 <- hill_taxa(comm = abundances_rlg_wide[,2:142],  q = 0)
# div_data_rlg$rlg.hill.1 <- hill_taxa(comm = abundances_rlg_wide[,2:142],  q = 1)
# div_data_rlg$rlg.hill.2 <- hill_taxa(comm = abundances_rlg_wide[,2:142],  q = 2)
# 
# 
# # add Hill numbers
# diversity <- left_join(diversity, div_data_rlg, by = "id.site")
# 
# 



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Calculation Transect based Data ###########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("abundances_3", "abundances_t", "sites", "traits", "diversity")))


## 0 - Preparation -------------------------------------------------------------

# Create `id.site` by removing the last three characters
abundances_t_1 <- abundances_t %>%
  mutate(id.site = substr(id.plot, 1, 5)) %>% 
  # condense to site based list
  # distinct(id.site, name.plant, cf,  .keep_all = T) %>% 
  select(id.site, id.plot, everything(), -cover.adapt.brbl, -remarks)


# prepare name plants plot based

#' 1. handle cf and not cf of same species on the the same plot
#' --> are considered as same species
#' 2. handle family level taxa
#' --> delete to decrease pseudoturnover 
#' 3. handle genus level species
#' --> keep them if they don't occur with species of the same genus on site
#' --> delete others if only few cases
#' 4. handle aggregates
#' 4.1 assign aggregate level for each species (new column "name.plant.agg")
#' --> species is in an aggregate: aggregate name
#' --> species is not in an aggregate: species name
#' 4.3 handle cases where target and non-target species are in the same aggregate
#' --> one species is target species --> all species of agg become target species
#' 
#' 
#' calculate diversity with aggregate level because of different 
#' considerations of agg level in regions (decrease of pseudoturnover by aggregating,
#' sensu Boch et al. 2022)
#' 


### handle cf ----
## are there cf and not cf of same species at same plot?
abundances_t_1 %>%
  select(id.plot, name.plant, cf) %>% 
  group_by(id.plot, name.plant) %>%
  filter(n_distinct(cf) > 1) %>%
  arrange(id.plot)
# no entries

abundances_t_1 %>%
  select(id.site, name.plant, cf) %>% 
  group_by(id.site, name.plant) %>%
  filter(n_distinct(cf) > 1) %>%
  arrange(id.site)
# 58 entries
# cf are considered as same species, plot based no need to combine


# check other duplicates
abundances_t_1 %>% 
  group_by(id.plot, name.plant) %>% 
  filter(n_distinct(name.plant) > 1)
# no duplicates


### handle family taxa ----

## are there taxa identified up to family level?
abundances_t_1 %>% 
  filter(TaxonRank == "FAM")
# four cases (M_KAL: Poaceae, M_TAN: Brassicaceae, N_BLU: Brassicaceae, N_NLG: Poaceae)

abundances_t_1 %>% 
  filter(id.site == "M_KAL") %>% 
  print(n = 50)
# other Poaceae at same site

abundances_t_1 %>% 
  filter(id.site == "M_TAN") %>% 
  print(n = 50)
# other Brassiceae at same site

abundances_t_1 %>% 
  filter(id.site == "N_BLU") %>% 
  print(n = 50)
# no other Brassiceae at same site

abundances_t_1 %>% 
  filter(id.site == "N_NLG") %>% 
  print(n = 50)
# other Poaceae at same site


# --> delete cases because unclear if distinct species and very little cases,
# therefore irrelevant

abundances_t_2 <- abundances_t_1 %>% 
  filter(TaxonRank != "FAM")



### handle genus taxa ----

## are there genus-level and minimum one species-level species at the same site?

# Create a dataframe with only the genus taxa
genus_data <- abundances_t_2 %>%
  filter(str_count(name.plant, "\\s") == 0)

# Create a dataframe with only the species taxa
species_data <- abundances_t_2 %>%
  filter(str_count(name.plant, "\\s") >= 1)

# Join the two dataframes on matching id.site and name.plant
matched_data <- genus_data %>%
  inner_join(species_data, by = "id.site",
             suffix = c(".genus", ".species"),
             relationship = "many-to-many") %>%
  filter(name.plant.genus == word(name.plant.species, 1)) %>%
  select(id.site, name.plant.genus, name.plant.species) %>% 
  distinct(id.site, .keep_all = TRUE)
# 34 cases where minimum one species and only-genus at the same site

# delete cases 
# because unclear if same species 
# and because irrelevant (very little cases)
abundances_t_3 <- abundances_t_2 %>% 
  anti_join(matched_data, by = c("id.site", "name.plant" = "name.plant.genus"))


### aggregation ----

# clean mistakes of data entry (NA in IsChildTaxonOf)
abundances_t_3 %>% 
  filter(is.na(IsChildTaxonOf))
abundances_t_4 <- abundances_t_3 %>% 
  mutate(IsChildTaxonOf = case_when(
    id.site == "N_UET" & name.plant == "Bromus hordeaceus" ~ "Bromus hordeaceus agg.",
    id.site == "M_ALT" & name.plant == "Stellaria apetala" ~ "Stellaria media agg.",
    id.site %in% c("M_SSB", "M_TTH", "M_TTP") & name.plant == "Galatella linosyris" ~ "Galatella",
    id.site %in% c("M_HAH", "M_TTH") & name.plant == "Inula conyzae" ~ "Inula",
    .default = IsChildTaxonOf
  ))

#' assign aggregate level for each species (new column "name.plant.agg")
#' --> species is in an aggregate: aggregate name
#' --> species is not in an aggregate: species name
abundances_t_4 <- abundances_t_4 %>% 
  mutate(name.plant.agg = case_when(
    str_detect(IsChildTaxonOf, "agg.") ~ IsChildTaxonOf,
    .default = name.plant
  ))


# change names manually (aggregates not detected automatically)
# Epilobium lamyi -> tetragonum s. l.
# Festuca pulchra -> valesiaca s. l. -> ovina agg.
# Luzula multiflora s. str. -> s. l. -> campestris agg.
# Gymnadenia densiflora -> conopsea s. l.
# Gymnadenia conopsea s. st. -> conopsea s. l.
# Ornithogalum umbellatum is part of Ornithogalum umbellatum agg. (sensu Rothmaler)
# Vicia sativa is part of Vicia sativa agg. (sensu Rothmaler)
abundances_t_4 <- abundances_t_4 %>% 
  mutate(name.plant.agg = case_when(
    name.plant.agg == "Epilobium lamyi" ~ "Epilobium tetragonum s. l.",
    name.plant.agg == "Festuca pulchra" ~ "Festuca ovina agg.",
    name.plant.agg == "Luzula multiflora s. str." ~ "Luzula campestris agg.",
    name.plant.agg == "Gymnadenia densiflora" ~ "Gymnadenia conopsea s. l.",
    name.plant.agg == "Gymnadenia conopsea s. str." ~ "Gymnadenia conopsea s. l.",
    name.plant.agg == "Vicia sativa" ~ "Vicia sativa agg.",
    name.plant.agg == "Ornithogalum umbellatum" ~ "Ornithogalum umbellatum agg.",
    .default = name.plant.agg
  ))


# check "doubles" on sites due to aggregation 
# (> 5 --> not all doubles included!)
abundances_t_4 %>% 
  select(id.site, name.plant, name.plant.agg) %>% 
  arrange(id.site, name.plant.agg) %>% 
  add_count(id.site, name.plant.agg) %>% 
  filter(n > 5)
# several cases (74 entries)

# check "doubles" on plots due to aggregation
abundances_t_4 %>% 
  select(id.plot, name.plant, name.plant.agg) %>% 
  arrange(id.plot, name.plant.agg) %>% 
  add_count(id.plot, name.plant.agg) %>% 
  filter(n > 1)
# several cases (105 entries)


# handle cases where target and non-target species are in the same aggregate
# filter aggregates with target and non-target species
agg_target <- abundances_t_4 %>%
  select(id.site, name.plant, name.plant.agg, target.species) %>% 
  group_by(name.plant.agg) %>%
  filter(n_distinct(target.species) > 1) %>%
  ungroup() %>% 
  arrange(name.plant.agg) %>% 
  distinct(name.plant, name.plant.agg, .keep_all = T)
# several cases (n = 45)


# define all species of an aggregate the same target species value
# one species is target species --> all species of agg are target species
abundances_t_4 <- abundances_t_4 %>% 
  mutate(target.species.agg = case_when(
    name.plant.agg %in% agg_target$name.plant.agg ~ 1,
    .default = target.species
  ))

# check
abundances_t_4 %>%
  select(id.site, name.plant, name.plant.agg, target.species, target.species.agg) %>% 
  group_by(name.plant.agg) %>%
  filter(n_distinct(target.species) > 1) %>%
  ungroup() %>% 
  distinct(name.plant, name.plant.agg, .keep_all = T) %>%
  arrange(name.plant.agg, target.species) %>% 
  print(n= 45)
#%>%
# write_csv(
#   here(
#     "outputs", "tables", "list_conflict_target_20241111.csv"
#   ))
# ok
# write_csv(
#   here(
#     "outputs", "tables", "target_plants_agg_20241111.csv"
#   )
# )



abundances_t_4 <- abundances_t_4 %>% 
  select(id.site, name.plant, name.plant.agg,
         target.species, target.species.agg, everything())


# list with all taxa on species level
name.plant_list_total <- abundances_t_4 %>% 
  distinct(name.plant, .keep_all = T) %>% 
  select(name.plant, IsChildTaxonOf, name.plant.agg, target.species, target.species.agg) %>%
  arrange(name.plant)# %>% 
# write_csv(
#   here(
#     "outputs", "tables", "name.plant_transect_taxa_20241126.csv"
#   )
# )

# list with all taxa on agg level
name.plant.agg_list_total <- abundances_t_4 %>% 
  distinct(name.plant.agg, .keep_all = T) %>% 
  select(name.plant.agg, target.species.agg) %>% 
  arrange(name.plant.agg)# %>% 
# write_csv(
#   here(
#     "outputs", "tables", "name.plant.agg_transect_taxa_20241126.csv"
#   )
# )



## 1 - Total species richness -------------------------------------------

# use aggregated species
# species at genus and family level are counted (they don't occur on same site with 
# taxa at species-level)


richness_tot <- abundances_t_4 %>% 
  group_by(id.site) %>% 
  summarise(transect_tot = n_distinct(name.plant.agg))


## 2 - Target species richness -------------------------------------------

# use aggregated species
# species at genus and family level are counted (they don't occur on same site with 
# taxa at species-level)


richness_target <- abundances_t_4 %>% 
  filter(target.species.agg == 1) %>% 
  group_by(id.site) %>% 
  summarise(transect_target = n_distinct(name.plant.agg))


## 3 - Forb number -------------------------------------------

# check NA
abundances_t_4 %>% 
  filter(is.na(funct.group))
# no herbs with NA

forb_no.site <- abundances_t_4 %>% 
  filter(funct.group == "herb") %>% 
  group_by(id.site) %>% 
  summarise(forb.no.site = n_distinct(name.plant.agg))


## 4 - Diversity table --------------------------------------------------------

diversity_transect <- richness_tot %>% 
  left_join(richness_target, by = "id.site") %>% 
  left_join(forb_no.site, by = "id.site")


## 5 - Forb Index -------------------------------------------

# = forb species richness * proportion of forbs on cumulative total plant cover
# calculate per plot, then mean per site


# take sum of cover into new variable "cover.mean.agg"
forb_index <- abundances_t_4 %>% 
  group_by(id.plot, name.plant.agg) %>%
  mutate(cover.agg = sum(cover)) %>% 
  ungroup()

# check
forb_index %>% 
  select(id.plot, name.plant, name.plant.agg, cover, cover.agg) %>% 
  filter(cover != cover.agg) %>% 
  arrange(id.plot, name.plant.agg) %>% 
  print(n =100)
# ok

# delete "doubles" at plots
forb_index <- forb_index %>% 
  distinct(id.plot, name.plant.agg, .keep_all = T)


#plant cumulative cover + rel.cover per vegetation plot/subtransect
forb_index <- forb_index %>%
  filter(!str_ends(id.plot, "_T")) %>%
  group_by(id.plot) %>%
  mutate(
    # plant cumulative cover per plot
    cum.cover.plot = sum(cover.agg),
    # rel.cover of plant per plot
    rel.cover.plot = round(cover.agg/cum.cover.plot, 2)) 


# forb no + forb proportion of cum.cover per vegetation plot
forb_index <- forb_index %>% 
  filter((funct.group %in% c("herb"))) %>%
  group_by(id.plot) %>% 
  mutate(forb.no.plot = n_distinct(name.plant.agg),
         forb.rel.cover.plot = sum(rel.cover.plot),
         forb.index.plot = forb.no.plot * forb.rel.cover.plot) 

# calculate forb.index per site (based on only subtransect)
forb_index <- forb_index %>% 
  ungroup() %>%
  select(id.site, id.plot, forb.index.plot, forb.rel.cover.plot) %>%
  distinct() %>%
  group_by(id.site) %>%
  # 4 plots don't have herbs --> 744 plots --> calculate mean by sum divided by 4, not with mean()
  summarise(forb.index.site = round(sum(forb.index.plot)/4,2),
            forb.rel.cover.site = round(sum(forb.rel.cover.plot)/4,2)) 






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Export ####################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("sites", "traits", "diversity", "diversity_transect",
                          "forb_index", "abundances_3", "abundances_t_4")))



# Species Diversity A1-A4

diversity %>% 
  write_csv(
    here(
      "data", "processed", "data_processed_plants_site_diversity_20250306.csv"
    ))


# Species Diversity Transect (A1-A4 and T)

diversity_transect %>% 
  write_csv(
    here(
      "data", "processed", "data_processed_plants_transect_diversity_20250306.csv"
    ))


# Forb Index (A1-A4)

forb_index %>% 
  write_csv(
    here(
      "data", "processed", "data_processed_plants_forbindex_20250306.csv"
    ))



# abundance table A1-A4
abundances_3 %>%
  write_csv(
    here(
      "outputs", "temp", "vegetation", "abundances_A1A4_20250306.csv"
    )
  )


# abundance table transect
abundances_t_4 %>%
  select(id.plot, id.site, everything()) %>%
  write_csv(
    here(
      "outputs", "temp", "vegetation", "abundances_t_20250306.csv"
    )
  )



# Characteristic Species list

# data <- abundances_3 %>% 
#   filter(target.species.agg == 1) %>% 
#   distinct(name.plant, .keep_all = T) %>% 
#   select(name.plant, name.plant.agg, target.species, target.species.agg,
#          r.all.diagnostic, ffh.lrt, calth.btt) %>% 
#   arrange(name.plant.agg, name.plant) %>% 
#   write_csv(
#     here(
#       "outputs", "tables", "target_plants_list_20241119.csv"
#     )
#   )


# for Alina (RII paper)
# list of species and number of sites they're occurring
# in regions and hydrology

list_plants <- abundances_tot_agg %>%
  left_join(sites, by = "id.site") %>% 
  select(id.site, name.plant, name.plant.agg, target.species, target.species.agg,
         region, hydrology, site.type)  %>%
  filter(site.type == "restored") %>%
  # subset of 77 sites
  filter(id.site %in% c("M_ASE", "M_BAD", "M_BUH", "M_CAL", "M_DOB", "M_DOR",
                        "M_FRB", "M_HAI", "M_HIR", "M_HNB", "M_ILB", "M_JER",
                        "M_KOT", "M_KUH", "M_NEL", "M_PRA", "M_SKF", "M_STA",
                        "M_STF", "M_STN", "M_TAN", "M_TTH", "M_TTP", "M_VOC",
                        "M_WEF", "M_WIM", "N_ALA", "N_BEB", "N_BIT", "N_BRD",
                        "N_BRE", "N_BRO", "N_CNI", "N_DAM", "N_GAM", "N_GAR",
                        "N_GUT", "N_HAI", "N_HOH", "N_HOY", "N_JAS", "N_KAP",
                        "N_LUB", "N_MOO", "N_OCH", "N_PAN", "N_RET", "N_TUT",
                        "N_VOG", "N_WAL", "N_WEN", "N_ZDF", "S_ACR", "S_AUH",
                        "S_BUC", "S_ECH", "S_FCH", "S_FHZ", "S_GIG", "S_GSH",
                        "S_GSR", "S_GTZ", "S_HLZ", "S_HTH", "S_HUB", "S_KRT",
                        "S_MBM", "S_OST", "S_OSW", "S_PEI", "S_SBN", "S_SHH",
                        "S_SNF", "S_STB", "S_WGD", "S_WNK", "S_WTZ")) %>% 
  arrange(region, hydrology, name.plant.agg, name.plant) %>%
  group_by(region, hydrology, name.plant.agg, target.species.agg) %>%
  summarise(unique_sites = n_distinct(id.site), .groups = "drop") %>% 
  pivot_wider(
    names_from = c(hydrology, region),
    values_from = unique_sites,
    values_fill = 0
  ) %>% 
  arrange(name.plant.agg) %>% 
  write_csv(
    here(
      "outputs", "temp", "vegetation", "list_plants_sites_by_region_hydrology_20241205.csv"
    )
  )





