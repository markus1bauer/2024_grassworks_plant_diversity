#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Output table: Model results
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

load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_fcsihill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_fgratio.Rdata"))


# remove dispersion component from tidy dataframe
mod_list <- modelsummary(restref_tothill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_tothill0 <- mod_list

mod_list <- modelsummary(restref_targethill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_targethill0 <- mod_list

mod_list <- modelsummary(restref_tothill1, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_tothill1 <- mod_list

mod_list <- modelsummary(restref_targethill1, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_targethill1 <- mod_list

mod_list <- modelsummary(restref_fcsihill0, output = "modelsummary_list")
mod_list$tidy <- mod_list$tidy %>% 
  filter(component != "dispersion")
ml_fscihill0 <- mod_list

mod_list <- modelsummary(restref_fgratio, output = "modelsummary_list")
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
        "site.typenegative"    = "Negative reference",
        "site.typerestored"    = "Restored sites",
        "site.typepositive" = "Positive reference",
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
               "outputs", "statistics", "vegetation", "model_summary_restref.docx")
             )

            

# Version 2 - not used
modelsummary(m,
             fmt = fmt_statistic(estimate = 2, std.error = 2, statistic = 2, p.value = 3),
             stars = TRUE,
             statistic = c("z" = "statistic"),
             coef_map = cm,
             gof_map = gm,
             output = here(
               "outputs", "statistics", "vegetation", "model_summary_restref_2.docx")
             )



# end skript




# code -------------------------------------------------------------------------



tab <- modelsummary(m,
                    output = "gt",
                    shape = effect + term ~ model + statistic,
                    statistic = c("std.error", "statistic"
                                  # , "{p.value}{stars}"
                    ),
                    gof_omit = 'RMSE|BIC',
                    stars = TRUE,
                    coef_map = cm,
                    # include_reference = TRUE
)

# TODO google: change objects in lists
tab$`_boxhead`$column_label

gtsave(tab, here("outputs", "tables", "vegetation", "model_summary_restref.html"))
gtsave(tab, here("outputs", "tables", "vegetation", "model_summary_restref.png"))

library(gt)
tab %>% 
  # tab_style(
  #   locations = cells_column_labels(),
  #   style = cell_borders(
  #     sides = "top", color = "black", style = "solid", weight = px(1)
  #   )
  # ) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(
              columns = starts_with("Est."),
              # rows = "Est." == "0.538***"
              rows = contains("*")
              ))
tab %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = names(tab)[grep("Est.", names(tab))] %>% 
              lapply(FUN = \(col_name){
                cells_body(columns = col_name,
                           rows = grepl('/*', tab[[col_name]])
                )
              })
  )

tab <- tidy(m_tothill0) %>% 
  gt()

my_data |> 
  gt() |>
  tab_style(
    style = list(
      cell_fill(color = 'red'),
      cell_text(weight = 'bold')),
    locations = names(my_data)[grep('outcome', names(my_data))] |>
      lapply(FUN = \(col_name){
        cells_body(columns = col_name,
                   rows = grepl('bad', my_data[[col_name]])
        )
      })
  )
)


model <- lm(mpg ~ wt + hp, data = mtcars)

# Extract regression results
results <- tidy(model)

tidy(m_tothill0) %>%
  mutate(estimate = ifelse(p.value < 0.05, 
                           paste0("**", round(estimate, 3), "**"), 
                           round(estimate, 3))) %>%
  gt() %>% 
  fmt_markdown(columns = "estimate") # Format Markdown for bold text

modelsummary(m_tothill0,
             shape = term ~ model + statistic,
)



library(lme4)
restref_tothill0_1 <- glmer.nb(tot.hill.0 ~ site.type
                               + (1|region) + (1|hydrology), data = data_model_tot)
restref_targethill0_1 <- glmer.nb(target.hill.0 ~ site.type
                                  + (1|region) + (1|hydrology), data = data_model_target)


get_estimates(restref_tothill0)



f <- function(x) format(round(x, 3), big.mark=",")
gm <- list(
  list("raw" = "estimate", "clean" = "Estimate", "fmt" = f),
  list("raw" = "std.error", "clean" = "SE", "fmt" = f),
  list("raw" = "statistic", "clean" = "Z-value", "fmt" = f),
  list("raw" = "p.value", "clean" = "p", "fmt" = f)
  )


modelsummary(restref_tothill0, 
             # component = "conditional",
             # shape = term ~ model + statistic,
             # gof_map = gm
             # coef_omit = where(component == "dispersion"),
             coef_omit = c("Intercept"),
             coef_rename = coef_rename
             # estimate = "{estimate} {p.value}"
)

modelsummary(restref_targethill0)


library(broom.mixed)
tidy_model <- broom.mixed::tidy(restref_tothill0)
tidy_model <- as.data.frame(tidy_model)
modelsummary(tidy_model)

library(broom)
tidy_model <- broom::tidy(restref_tothill0)
library(parameters)
param <- parameters(restref_tothill0, component = "conditional") %>% 
  filter(Component == "conditional")
modelsummary(restref_tothill0,
             shape = term + component ~ model + statistic,
             # component = "conditional"

                          )
options(modelsummary_get = "easystats")


options(modelsummary_get = "broom.mixed")
summ <- modelsummary(restref_tothill0,
             # shape = term ~ model + statistic,
             # Component = "conditional"
)

modelsummary(m_tothill0)             
sigma(restref_tothill0)





# Custom tidy function to filter out the dispersion term
custom_tidy <- function(model) {
  broom.mixed::tidy(model) %>%
    filter(Component == "conditional")
}

# Use the custom tidy function in modelsummary
modelsummary(restref_tothill0,
             tidy_fun = custom_tidy,
             
             # shape = term ~ model + statistic
             
             )










m <- list(
  "Total Species Richness" = restref_tothill0,
  "Target Species Richness" = restref_targethill0)

cm <- c("(Intercept)"    = "(Intercept)",
        "site.typerestored"    = "Site type = restored",
        "site.typepositive" = "Site type = positive",
        "SD (Intercept region)" = "Region",
        "SD (Intercept hydrology)" = "Hydrology")

tab <- modelsummary(m,
                    output = "gt",
                    shape = term + effect ~ model + statistic,
                    statistic = c("std.error", "statistic"
                                  # , "{p.value}{stars}"
                    ),
                    gof_omit = 'RMSE|BIC',
                    stars = TRUE,
                    coef_map = cm,
                    # include_reference = TRUE
)






library(sjPlot)
library(sjmisc)
library(sjlabelled)

tab_model(restref_tothill0, restref_targethill0,
          show.ci = F, show.se = T, show.aic = T, show.stat = T,
          p.style = "stars",
          CSS = css_theme("cells"))















modelsummary(restref_tothill0,
             shape = term + component ~ model + statistic,
             # statistic = NULL
             # coef_omit = where(component == "dispersion"),
             # estimate = "{estimate} {p.value}"
             )

modelsummary(m,
             shape = term ~ model + statistic,
             statistic = "conf.int",
             gof_map = NA)
