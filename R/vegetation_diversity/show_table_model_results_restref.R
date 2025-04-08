#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Output table: Model results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PREPARATION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(broom.mixed) # Tidy up the model summary
library(performance)
library(gt)
library(glmmTMB)




### Start ###
rm(list = ls())



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 1 - Total Species Richness #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_tothill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_tothill0) %>% 
  add_column(nobs = nobs(restref_tothill0))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_1 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_tothill0.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 2 - Characteristic Species Richness ###################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_targethill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_targethill0) %>% 
  add_column(nobs = nobs(restref_targethill0))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_2 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_targethill0.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 3 - Total Hill-Shannon #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(file = here("outputs", "models", "vegetation", "model_plants_restref_tothill1.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_tothill1) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_tothill1) %>% 
  add_column(nobs = nobs(restref_tothill1))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_3 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_tothill1.csv"
#     ),
#     delim = ";"
#   )






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 4 - Characteristic Hill-Shannon ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill1.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_targethill1) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_targethill1) %>% 
  add_column(nobs = nobs(restref_targethill1))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_4 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_targethill1.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 5 - FCSi ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(file = here("outputs", "models", "vegetation", "model_plants_restref_fcsihill0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_fcsihill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_fcsihill0) %>% 
  add_column(nobs = nobs(restref_fcsihill0))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_5 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_fcsihill0.csv"
#     ),
#     delim = ";"
#   )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Table 6 - Forb-Grass Ratio ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

load(file = here("outputs", "models", "vegetation", "model_plants_restref_fgratio.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restref_fgratio) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 

# extract goodness-of-fit statistics
gof <- model_performance(restref_fgratio) %>% 
  add_column(nobs = nobs(restref_fgratio))

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_smry_6 <- m_est %>% 
  bind_rows(gof_t)

# save
# model_smry %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_smry_restref_fgratio.csv"
#     ),
#     delim = ";"
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Identify objects to keep (those starting with "model")
objects_to_keep <- grep("^model", ls(), value = TRUE)

# Remove all other objects
rm(list = setdiff(ls(), objects_to_keep))


model_smry_joint <- model_smry_1 %>% 
  left_join(model_smry_2 %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_smry_3 %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_smry_4 %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_smry_5 %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_smry_6 %>% select(term, estimate, std.error, statistic, p.value),
            by = "term")

# save table with all summary results
# model_smry_joint %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restref_all.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# gt table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# model_smry_joint <- read_delim(
#   here("outputs", "statistics", "vegetation", "model_summary_restref_all.csv"),
#   col_names = TRUE, 
#   delim = ";"
# )





# Create the gt table
m_smry_gt <- model_smry_joint %>%
  select(everything(), -effect, -component, -group) %>% 
  mutate(
    term = case_when(
      term == "site.typerestored" ~ "Restored sites",
      term == "site.typepositive" ~ "Positive reference",
      term == "sd__(Intercept) region" ~ "Region (random)",
      term == "sd__(Intercept) hydrology" ~ "Hydrology (random)",
      term == "R2_conditional" ~ "R2 conditional",
      term == "R2_marginal" ~ "R2 marginal",
      term == "nobs" ~ "No. of observations",
      TRUE ~ term)) %>% 
  filter(!term %in% c("Score_log", "Score_spherical", "AICc", "BIC")) %>% 
  gt()



# Create a list of columns for estimate and p.value
colnames(model_smry_joint)
estimate_columns <- c("estimate.x", "estimate.y"
                      , "estimate.x.x", "estimate.y.y"
                      , "estimate.x.x.x", "estimate.y.y.y"
                      # , "estimate.x.x.x.x", "estimate.y.y.y.y"
                      # , "estimate"
                      )
pvalue_columns <- c("p.value.x", "p.value.y"
                    , "p.value.x.x", "p.value.y.y"
                    , "p.value.x.x.x", "p.value.y.y.y"
                    # , "p.value.x.x.x.x", "p.value.y.y.y.y"
                    # , "p.value"
                    )

# Loop through each estimate and p.value pair to apply bold styling if p.value < 0.05
for (i in seq_along(estimate_columns)) {
  m_smry_gt <- m_smry_gt %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(!!sym(estimate_columns[i])),  # Dynamically use the estimate column
        rows = !!sym(pvalue_columns[i]) < 0.05 & !is.na(!!sym(pvalue_columns[i]))  # Apply condition to the corresponding p.value column
      )
    )
}

for (i in seq_along(std_error_columns)) {
  m_smry_gt <- m_smry_gt %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(!!sym(std_error_columns[i])),  # Dynamically use the estimate column
        rows = !!sym(pvalue_columns[i]) < 0.05 & !is.na(!!sym(pvalue_columns[i]))  # Apply condition to the corresponding p.value column
      )
    )
}



# Apply tab_spanner for each group
# change x and y manually according to sequence in creating model_smry_joint
colnames(model_smry_joint)
column_groups <- list(
  "Total species richness" = c("estimate.x", "std.error.x", "statistic.x", "p.value.x"),
  "Char. species richness" = c("estimate.y", "std.error.y", "statistic.y", "p.value.y"),
  "Total Hill-Shannon" = c("estimate.x.x", "std.error.x.x", "statistic.x.x", "p.value.x.x"),
  "Char. Hill-Shannon" = c("estimate.y.y", "std.error.y.y", "statistic.y.y", "p.value.y.y"),
  "FCSi" = c("estimate.x.x.x", "std.error.x.x.x", "statistic.x.x.x", "p.value.x.x.x"),
  "Forb-grass ratio" = c("estimate.y.y.y", "std.error.y.y.y", "statistic.y.y.y", "p.value.y.y.y")
)
for (label in names(column_groups)) {
  m_smry_gt <- m_smry_gt %>%
    tab_spanner(
      label = label,
      columns = column_groups[[label]]
    )
}

# change column names
col_names <- colnames(model_smry_joint)
new_col_names <- col_names %>%
  gsub("^esti(.*)", "Est.", .) %>%
  gsub("^std(.*)", "S.E.", .) %>% 
  gsub("^stat(.*)", "z", .) %>% 
  gsub("^p.value(.*)", "p.value", .) 

names(new_col_names) <- col_names  # Map old to new names
m_smry_gt <- m_smry_gt %>%
  cols_label(!!!new_col_names[4:40])

# more visual changes
m_smry_gt <- m_smry_gt %>% 
  # hide p.value column
  cols_hide(columns = matches("^p.value")) %>% 
  # round decimal places
  fmt_number(
    columns = where(is.numeric),
    decimals = 2,
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    rows = term == "No. of observations", # Apply only to the row where term == "No. of observations"
    decimals = 0
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    rows = term %in% c("AIC", "AICc", "BIC"), # Apply only to the specific rows
    decimals = 0,
    use_seps = F
  ) %>% 
  sub_missing(
    columns = everything(), # Apply to all columns
    missing_text = "" # Replace NA with blank
  ) %>% 
  cols_label(term = "")


  



gtsave(m_smry_gt,
       here(
         "outputs", "statistics", "vegetation", "model_summary_restref_all_gt.docx"
       ))





### try out with p value stars




library(gt)

# Get all estimate and p-value column names dynamically
estimate_cols <- names(model_smry_joint)[str_detect(names(model_smry_joint), "^estimate")]
pvalue_cols <- names(model_smry_joint)[str_detect(names(model_smry_joint), "^p\\.value")]

# Ensure corresponding p-value column exists for each estimate column
mapping <- setNames(pvalue_cols, estimate_cols)


# Create the gt table
m_smry_gt <- model_smry_joint %>%
  select(everything(), -effect, -component, -group) %>% 
  mutate(
    term = case_when(
      term == "site.typerestored" ~ "Site type = Restored",
      term == "site.typepositive" ~ "Site type = Positive Reference",
      term == "sd__(Intercept) region" ~ "Region (random)",
      term == "sd__(Intercept) hydrology" ~ "Hydrology (random)",
      term == "R2_conditional" ~ "R2 conditional",
      term == "R2_marginal" ~ "R2 marginal",
      TRUE ~ term)) %>% 
  filter(!term %in% c("Score_log", "Score_spherical")) %>% 
  mutate(across(matches("^estimate"), ~ round(.x, digits = 3)),
         across(matches("^std"), ~ round(.x, digits = 3)),
         across(matches("^statistic"), ~ round(.x, digits = 3)),
         ) %>% 
  gt() %>% 
  { 
    tbl <- . # Store gt table
    for (est_col in names(mapping)) {
      pval_col <- mapping[[est_col]]
      tbl <- tbl %>%
        fmt(
          columns = all_of(est_col),
          fns = function(x) {
            stars <- case_when(
              model_smry_joint[[pval_col]] < 0.001 ~ "***",
              model_smry_joint[[pval_col]] < 0.01  ~ "**",
              model_smry_joint[[pval_col]] < 0.05  ~ "*",
              TRUE ~ ""
            )
            paste0(x, stars)
          }
        )
    }
    tbl # Return modified gt table
  } 





# Create a list of columns for estimate and p.value
colnames(model_smry_joint)
estimate_columns <- c("estimate.x", "estimate.y"
                      , "estimate.x.x", "estimate.y.y"
                      , "estimate.x.x.x", "estimate.y.y.y"
                      # , "estimate.x.x.x.x", "estimate.y.y.y.y"
                      # , "estimate"
)
pvalue_columns <- c("p.value.x", "p.value.y"
                    , "p.value.x.x", "p.value.y.y"
                    , "p.value.x.x.x", "p.value.y.y.y"
                    # , "p.value.x.x.x.x", "p.value.y.y.y.y"
                    # , "p.value"
)

# Loop through each estimate and p.value pair to apply bold styling if p.value < 0.05
for (i in seq_along(estimate_columns)) {
  m_smry_gt <- m_smry_gt %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = vars(!!sym(estimate_columns[i])),  # Dynamically use the estimate column
        rows = !!sym(pvalue_columns[i]) < 0.05 & !is.na(!!sym(pvalue_columns[i]))  # Apply condition to the corresponding p.value column
      )
    )
}





# Apply tab_spanner for each group
# change x and y manually according to sequence in creating model_smry_joint
colnames(m_smry)
column_groups <- list(
  "Total Species Richness" = c("estimate.x", "std.error.x", "statistic.x", "p.value.x"),
  "Total Hill-Shannon" = c("estimate.y", "std.error.y", "statistic.y", "p.value.y"),
  "Char. Species Richness" = c("estimate.x.x", "std.error.x.x", "statistic.x.x", "p.value.x.x"),
  "Char. Hill-Shannon" = c("estimate.y.y", "std.error.y.y", "statistic.y.y", "p.value.y.y"),
  "Forb-Grass Ratio" = c("estimate.x.x.x", "std.error.x.x.x", "statistic.x.x.x", "p.value.x.x.x"),
  "FCSi" = c("estimate.y.y.y", "std.error.y.y.y", "statistic.y.y.y", "p.value.y.y.y")
  # ,
  # "" = c("estimate.x.x.x.x", "std.error.x.x.x.x", "statistic.x.x.x.x", "p.value.x.x.x.x"),
  # "" = c("estimate.y.y.y.y", "std.error.y.y.y.y", "statistic.y.y.y.y", "p.value.y.y.y.y"),
  # "Ellenberg L value" = c("estimate", "std.error", "statistic", "p.value")
  
)
for (label in names(column_groups)) {
  m_smry_gt <- m_smry_gt %>%
    tab_spanner(
      label = label,
      columns = column_groups[[label]]
    )
}

# change column names
col_names <- colnames(model_smry_joint)
new_col_names <- col_names %>%
  gsub("^esti(.*)", "Est.", .) %>%
  gsub("^std(.*)", "S.E.", .) %>% 
  gsub("^stat(.*)", "z", .) %>% 
  gsub("^p.value(.*)", "p.value", .) 

names(new_col_names) <- col_names  # Map old to new names
m_smry_gt <- m_smry_gt %>%
  cols_label(!!!new_col_names[4:40])

# more visual changes
m_smry_gt <- m_smry_gt %>% 
  # hide p.value column
  cols_hide(columns = matches("^p.value")) %>% 
  # round decimal places
  fmt_number(
    columns = where(is.numeric),
    rows = term == "nobs", # Apply only to the row where term == "nobs"
    decimals = 0
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    rows = term %in% c("AIC", "AICc", "BIC"), # Apply only to the specific rows
    decimals = 1,
    use_seps = F
  ) %>% 
  sub_missing(
    columns = everything(), # Apply to all columns
    missing_text = "" # Replace NA with blank
  )



