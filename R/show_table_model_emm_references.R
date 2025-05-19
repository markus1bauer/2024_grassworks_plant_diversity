#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference
# Plot table: EMMs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(gt)
library(glmmTMB)




### Start ###
rm(list = ls())




## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_fgratio.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restref_fcsihill0.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculation of EMMs & SE ####################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
# type ="response" for back-transformation from log-scale

emm.rest.ref <- emmeans(restref_tothill0, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.tothill0.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")

emm.rest.ref <- emmeans(restref_tothill1, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.tothill1.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")

emm.rest.ref <- emmeans(restref_targethill0, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.targethill0.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")

emm.rest.ref <- emmeans(restref_targethill1, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.targethill1.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")

emm.rest.ref <- emmeans(restref_fgratio, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.fgratio.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")

emm.rest.ref <- emmeans(restref_fcsihill0, ~ site.type)
cld_results <- multcomp::cld(emm.rest.ref, adjust = "tukey", Letters = letters)
emm.restref.fcsihill0.df <- summary(emm.rest.ref, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(site.type, .group), by = "site.type")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# H - Joint table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

emm.table <- emm.restref.tothill0.df %>% 
  left_join(emm.restref.targethill0.df, by = "site.type") %>% 
  left_join(emm.restref.tothill1.df, by = "site.type") %>% 
  left_join(emm.restref.targethill1.df, by = "site.type") %>% 
  left_join(emm.restref.fcsihill0.df, by = "site.type") %>% 
  left_join(emm.restref.fgratio.df, by = "site.type")

emm.table.gt <- emm.table %>% 
  mutate(
    site.type = case_when(
      site.type == "negative" ~ "Negative reference",
      site.type == "restored" ~ "Restored sites",
      site.type == "positive" ~ "Positive reference",
      TRUE ~ site.type)) %>% 
  
  gt()

# Apply tab_spanner for each group
# change x and y manually according to sequence in creating table
colnames(emm.table)
column_groups <- list(
  "Total species richness" = c("response.x", "SE.x", ".group.x"),
  "Char. species richness" = c("response.y", "SE.y", ".group.y"),
  "Total Hill-Shannon" = c("response.x.x", "SE.x.x", ".group.x.x"),
  "Char. Hill-Shannon" = c("response.y.y", "SE.y.y", ".group.y.y"),
  "FCSi" = c("emmean", "SE.x.x.x", ".group.x.x.x"),
  "Forb-grass ratio" = c("response", "SE.y.y.y", ".group.y.y.y")
)


for (label in names(column_groups)) {
  emm.table.gt <- emm.table.gt %>%
    tab_spanner(
      label = label,
      columns = column_groups[[label]]
    )
}


# change column names
col_names <- colnames(emm.table)
new_col_names <- col_names %>%
  gsub("^resp(.*)", "EMM", .) %>%
  gsub("^SE(.*)", "S.E.", .) %>% 
  gsub("^emmean(.*)", "EMM", .) %>% 
  gsub("^.group(.*)", "sign.", .)
  
names(new_col_names) <- col_names  # Map old to new names
emm.table.gt <- emm.table.gt %>%
  cols_label(!!!new_col_names[2:40])

# more visual changes
emm.table.gt <- emm.table.gt %>% 
  # round decimal places
  fmt_number(
    columns = where(is.numeric),
    decimals = 1,
  ) %>% 
  cols_label(site.type = "")

emm.table.gt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# gt table
gtsave(emm.table.gt,
       here(
         "outputs", "statistics", "vegetation", "model_emm_restref.docx"
       ))
