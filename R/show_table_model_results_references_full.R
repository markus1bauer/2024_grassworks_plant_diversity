#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration effects
# Output table: Model results restfact full model
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A PREPARATION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(modelsummary)
library(pandoc)


### Start ###
rm(list = ls())

## load models -------------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio.Rdata"))


# remove dispersion component from tidy dataframe
mod_list <- modelsummary(restfact_tothill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_tothill0 <- mod_list

mod_list <- modelsummary(restfact_targethill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_targethill0 <- mod_list

mod_list <- modelsummary(restfact_tothill1, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_tothill1 <- mod_list

mod_list <- modelsummary(restfact_targethill1, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_targethill1 <- mod_list

mod_list <- modelsummary(restfact_fcsihill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_fscihill0 <- mod_list

mod_list <- modelsummary(restfact_fgratio, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_fgratio <- mod_list

# Identify objects to keep (those starting with "ml")
objects_to_keep <- grep("^ml", ls(), value = TRUE)

# Remove all other objects
rm(list = setdiff(ls(), objects_to_keep))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B CREATE TABLE ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

m <- list(
  "Total species richness" = ml_tothill0,
  "Char. species richness" = ml_targethill0,
  "Total Hill-Shannon" = ml_tothill1,
  "Char. Hill-Shannon" = ml_targethill1,
  "FCSi" = ml_fscihill0,
  "Forb-grass ratio" = ml_fgratio)

cm <- c("(Intercept)"    = "(Intercept)",
        "rest.methmga"    = "Management adaptation",
        "rest.methres"    = "Regional seed mixture",
        "rest.methdih" = "Direct harvesting",
        "land.use.histgrassland" = "Grassland",
        "rest.age.std" = "Restoration age",
        "SD (Intercept region)" = "Region (random)",
        "SD (Intercept hydrology)" = "Hydrology (random)",
        "R2 Marg." = "R2 marginal",
        "R2 Cond." = "R2 conditional")

# get_gof(ml_tothill0)
gm <- list(
  list("raw" = "nobs", "clean" = "N", "fmt" = 0),
  list("raw" = "r2.marginal", "clean" = "R2 marginal", "fmt" = 2),
  list("raw" = "r2.conditional", "clean" = "R2 conditional", "fmt" = 2),
  list("raw" = "aic", "clean" = "AIC", "fmt" = 0))


# Version 1
modelsummary(m,
                    fmt = fmt_statistic(estimate = 2, std.error = 2, statistic = 2, p.value = 3),
                    stars = TRUE,
                    estimate = c("Est.+/-S.E." = "{estimate}+/-{std.error}{stars}"),
                    # estimate = c("Est. (+/-S.E.)" = "{estimate}{stars} (+/-{std.error})"),
                    statistic = c("z" = "statistic"),
                    shape = term ~ model + statistic,
                    # shape = term + statistic ~ model + statistic,
                    coef_map = cm,
                    gof_map = gm,
                    output = here(
                      "outputs", "statistics", "vegetation", "model_summary_restfact_full.docx")
)



# Version 2 - not used
modelsummary(m,
             fmt = fmt_statistic(estimate = 2, std.error = 2, statistic = 2, p.value = 3),
             stars = TRUE,
             statistic = c("z" = "statistic"),
             coef_map = cm,
             gof_map = gm,
             output = here(
               "outputs", "statistics", "vegetation", "model_summary_restfact_full_2.docx")
)



# end skript






