#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Taxonomic plant diversity ####
# Prepare data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Markus Bauer
# 2025-04-08



### Packages ###
library(renv)
library(here)
library(tidyverse)

### Start ###
rm(list = ls())
# installr::updateR(
#   browse_news = FALSE,
#   install_R = TRUE,
#   copy_packages = TRUE,
#   copy_site_files = TRUE,
#   keep_old_packages = FALSE,
#   update_packages = FALSE,
#   start_new_R = FALSE,
#   quit_R = TRUE
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



sites <- read_csv(
  here("data", "raw", "data_processed_environment_nms_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, hydrology, region, rest.meth, land.use.hist, rest.age
  ) %>%
  mutate(
    region = fct_relevel(region, "north", "centre", "south"),
    hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
    rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih")
  )

diversity <- read_csv(
  here("data", "raw", "data_processed_plants_site_diversity_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Prepare data ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data <- sites %>%
  distinct() %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(
    across(where(is.numeric), ~ as.numeric(scale(.)), .names = "{col}.std")
    ) %>%
  left_join(diversity, by = "id.site")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data #######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



write_csv(
  data, here("data", "processed", "data_processed.csv")
)
