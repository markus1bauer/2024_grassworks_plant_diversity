#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
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

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill1.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio.Rdata"))
load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculation of EMMs & SE ####################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
# type ="response" for back-transformation from log-scale

emm.rest.meth <- emmeans(restfact_tothill0, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.tothill0.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")

emm.rest.meth <- emmeans(restfact_tothill1, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.tothill1.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")

emm.rest.meth <- emmeans(restfact_targethill0, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.targethill0.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")

emm.rest.meth <- emmeans(restfact_targethill1, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.targethill1.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")

emm.rest.meth <- emmeans(restfact_fgratio, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.fgratio.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")

emm.rest.meth <- emmeans(restfact_fcsihill0, ~ rest.meth)
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)
emm.restmeth.fcsihill0.df <- summary(emm.rest.meth, infer = F, type = "response") %>% select(-df) %>% 
  left_join(cld_results %>% select(rest.meth, .group), by = "rest.meth")



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# H - Joint table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

emm.table <- emm.restmeth.tothill0.df %>% 
  left_join(emm.restmeth.targethill0.df, by = "rest.meth") %>% 
  left_join(emm.restmeth.tothill1.df, by = "rest.meth") %>% 
  left_join(emm.restmeth.targethill1.df, by = "rest.meth") %>% 
  left_join(emm.restmeth.fcsihill0.df, by = "rest.meth") %>% 
  left_join(emm.restmeth.fgratio.df, by = "rest.meth")
  
emm.table.gt <- emm.table %>% 
  mutate(
    rest.meth = case_when(
      rest.meth == "cus" ~ "Cultivar seed mixture",
      rest.meth == "mga" ~ "Management adaptation",
      rest.meth == "res" ~ "Regional seed mixture",
      rest.meth == "dih" ~ "Direct harvesting",
      TRUE ~ rest.meth)) %>% 
  
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
  cols_label(rest.meth = "")

emm.table.gt


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# gt table
gtsave(emm.table.gt,
       here(
         "outputs", "statistics", "vegetation", "model_emm_restfact_full.docx"
       ))
