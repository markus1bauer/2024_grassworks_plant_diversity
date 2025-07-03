#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Taxonomic plant diversity
# Prepare data ####
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
    .default = "?",
    rest.age = "d"
  )) %>%
  rename(fg.ratio = site.forb.legu.grass.ratio) %>%
  mutate(
    region = fct_relevel(region, "north", "centre", "south"),
    hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
    rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih")
  ) %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(
    across(c(rest.age, fg.ratio), ~ as.numeric(scale(.)), .names = "{col}.std")
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
  select(
    id.site, site.type, hydrology, region, fg.ratio, fg.ratio.std,
    rest.meth, land.use.hist, rest.age, rest.age.std
    ) %>%
  distinct() %>%
  left_join(diversity, by = "id.site")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data #######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



write_csv(
  data, here("data", "processed", "data_processed.csv")
)
