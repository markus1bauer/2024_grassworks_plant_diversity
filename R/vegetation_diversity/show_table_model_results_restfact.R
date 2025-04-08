#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Output: Show Table Model Summary
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATION ###############################################################
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
# B - Table 1 - Total Species Richness #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_tothill0_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_tothill0_20y) %>% 
  add_column(nobs = nobs(restfact_tothill0_20y))

# library(modelsummary)
# gof <- get_gof(restfact_tothill0_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_1a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_tothill0_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_tothill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_tothill0) %>% 
  add_column(nobs = nobs(restfact_tothill0))

# library(modelsummary)
# gof <- get_gof(restfact_tothill0) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_1b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_tothill0.csv"
#     ),
#     delim = ";"
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Table 2 - Total Hill-Shannon #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill1_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_tothill1_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_tothill1_20y) %>% 
  add_column(nobs = nobs(restfact_tothill1_20y))

# library(modelsummary)
# gof <- get_gof(restfact_tothill1_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_2a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_tothill1_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill1.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_tothill1) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_tothill1) %>% 
  add_column(nobs = nobs(restfact_tothill1))

# library(modelsummary)
# gof <- get_gof(restfact_tothill1) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_2b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_tothill1.csv"
#     ),
#     delim = ";"
#   )

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D - Table 3 - Total Diff q1-q0 ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_diffq1q0_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_diffq1q0_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_diffq1q0_20y) %>% 
  add_column(nobs = nobs(restfact_diffq1q0_20y))

# library(modelsummary)
# gof <- get_gof(restfact_diffq1q0_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_3a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_diffq1q0_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_diffq1q0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_diffq1q0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_diffq1q0) %>% 
  add_column(nobs = nobs(restfact_diffq1q0))

# library(modelsummary)
# gof <- get_gof(restfact_diffq1q0) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_3b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_diffq1q0.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E - Table 4 - Characteristic Species Richness ################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targethill0_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targethill0_20y) %>% 
  add_column(nobs = nobs(restfact_targethill0_20y))

# library(modelsummary)
# gof <- get_gof(restfact_targethill0_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_4a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targethill0_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0.Rdata"))



# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targethill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targethill0) %>% 
  add_column(nobs = nobs(restfact_targethill0))

# library(modelsummary)
# gof <- get_gof(restfact_targethill0) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_4b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targethill0.csv"
#     ),
#     delim = ";"
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# F - Table 5 - Characteristic Hill-Shannon ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill1_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targethill1_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targethill1_20y) %>% 
  add_column(nobs = nobs(restfact_targethill1_20y))

# library(modelsummary)
# gof <- get_gof(restfact_targethill1_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_5a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targethill1_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill1.Rdata"))



# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targethill1) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targethill1) %>% 
  add_column(nobs = nobs(restfact_targethill1))

# library(modelsummary)
# gof <- get_gof(restfact_targethill1) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_5b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targethill1.csv"
#     ),
#     delim = ";"
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# G - Table 6 - Characteristic Diff q1-q0 ######################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targetdiffq1q0_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targetdiffq1q0_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targetdiffq1q0_20y) %>% 
  add_column(nobs = nobs(restfact_targetdiffq1q0_20y))

# library(modelsummary)
# gof <- get_gof(restfact_targetdiffq1q0_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_6a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targetdiffq1q0_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targetdiffq1q0.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_targetdiffq1q0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_targetdiffq1q0) %>% 
  add_column(nobs = nobs(restfact_targetdiffq1q0))

# library(modelsummary)
# gof <- get_gof(restfact_targetdiffq1q0) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_6b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_targetdiffq1q0.csv"
#     ),
#     delim = ";"
#   )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# H - Table 7 - Forb-Grass Ratio ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_fgratio_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_fgratio_20y) %>% 
  add_column(nobs = nobs(restfact_fgratio_20y))

# library(modelsummary)
# gof <- get_gof(restfact_fgratio_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_7a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_fgratio_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio.Rdata"))



# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_fgratio) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_fgratio) %>% 
  add_column(nobs = nobs(restfact_fgratio))

# library(modelsummary)
# gof <- get_gof(restfact_fgratio) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_7b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_fgratio.csv"
#     ),
#     delim = ";"
#   )




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# I - Table 8 - FCSi ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_fcsihill0_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_fcsihill0_20y) %>% 
  add_column(nobs = nobs(restfact_fcsihill0_20y))

# library(modelsummary)
# gof <- get_gof(restfact_fcsihill0_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_8a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_fcsihill0_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0.Rdata"))



# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_fcsihill0) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_fcsihill0) %>% 
  add_column(nobs = nobs(restfact_fcsihill0))

# library(modelsummary)
# gof <- get_gof(restfact_fcsihill0) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_8b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_fcsihill0.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# J - Table 9 - Ellenberg L value #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list = ls())

## a - <20 years model ---------------------------------------------------------


load(file = here("outputs", "models", "vegetation", "model_plants_restfact_ellenbergL_20y.Rdata"))


# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_ellenbergL_20y) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_ellenbergL_20y) %>% 
  add_column(nobs = nobs(restfact_ellenbergL_20y))

# library(modelsummary)
# gof <- get_gof(restfact_ellenbergL_20y) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_9a <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_ellenbergL_20y.csv"
#     ),
#     delim = ";"
#   )


## b - Final model -------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_ellenbergL.Rdata"))



# extract estimate statisitcs
m_est <- broom.mixed::tidy(restfact_ellenbergL) %>% 
  mutate(term = case_when(effect == "ran_pars" ~ paste(term, group),
                          .default = term)) #%>% 
# mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
# mutate(across(c(estimate, std.error, statistic), \(x) round(x, 3))) #%>%   # Rounding values
# mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
# mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                            .default = p.value))


# extract goodness-of-fit statistics
gof <- model_performance(restfact_ellenbergL) %>% 
  add_column(nobs = nobs(restfact_ellenbergL))

# library(modelsummary)
# gof <- get_gof(restfact_ellenbergL) %>% 
#   mutate(across(c(aic, bic, icc), \(x) round(x, 1))) %>%   # Rounding values
#   mutate(across(c(r2.conditional, r2.marginal), \(x) round(x, 3)))   # Rounding values

gof_t <- as.data.frame(t(gof)) %>% 
  rownames_to_column(var = "term") %>% 
  rename(estimate = "V1")

# combine into one table
model_sum_9b <- m_est %>% 
  bind_rows(gof_t)

# save
# model_sum %>% 
#   write_delim(
#     here(
#       "outputs", "statistics", "vegetation", "model_summary_restfact_ellenbergL.csv"
#     ),
#     delim = ";"
#   )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# H - Joint table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Identify objects to keep (those starting with "model")
objects_to_keep <- grep("^model", ls(), value = TRUE)

# Remove all other objects
rm(list = setdiff(ls(), objects_to_keep))


# Table models < 20 years
model_smry_joint_20y <- model_sum_1a %>% 
  left_join(model_sum_4a %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_sum_2a %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  # left_join(model_sum_3a %>% select(term, estimate, std.error, statistic, p.value),
  #           by = "term") %>% 
  left_join(model_sum_5a %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  # left_join(model_sum_6a %>% select(term, estimate, std.error, statistic, p.value),
  #           by = "term") %>% 
  left_join(model_sum_8a %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_sum_7a %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") #%>% 
# left_join(model_sum_9a %>% select(term, estimate, std.error, statistic, p.value),
#           by = "term")

model_smry_joint_20y %>% 
  write_delim(
    here(
      "outputs", "statistics", "vegetation", "model_summary_restfact_20y_all.csv"
    ),
    delim = ";"
  )

  
model_smry_joint <- model_sum_1b %>% 
  left_join(model_sum_4b %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  left_join(model_sum_2b %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  # left_join(model_sum_3b %>% select(term, estimate, std.error, statistic, p.value),
  #           by = "term") %>% 
  left_join(model_sum_5b %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>% 
  # left_join(model_sum_6b %>% select(term, estimate, std.error, statistic, p.value),
  #           by = "term") %>% 
  left_join(model_sum_8b %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") %>%
  left_join(model_sum_7b %>% select(term, estimate, std.error, statistic, p.value),
            by = "term") #%>% 
   
  # left_join(model_sum_9b %>% select(term, estimate, std.error, statistic, p.value),
  #           by = "term") 
  
model_smry_joint %>% 
  write_delim(
    here(
      "outputs", "statistics", "vegetation", "model_summary_restfact_final_all.csv"
    ),
    delim = ";"
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# I - gt table #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## a - final model -------------------------------------------------------------


# model_smry_joint <- read_delim(
#   here("outputs", "statistics", "vegetation", "model_summary_restfact_final_all.csv"),
#   col_names = TRUE, 
#   delim = ";"
# )



# Create the gt table
m_smry_gt <- model_smry_joint %>%
  select(everything(), -effect, -component, -group) %>% 
  mutate(
    term = case_when(
      term == "rest.methmga" ~ "Rest. method = Management adaptation",
      term == "rest.methres" ~ "Rest. method = Regional seed mixture",
      term == "rest.methdih" ~ "Rest. method = Direct harvesting",
      term == "land.use.histgrassland" ~ "Previous land use = Grassland",
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
std_error_columns <- c("std.error.x", "std.error.y"
                       , "std.error.x.x", "std.error.y.y"
                       , "std.error.x.x.x", "std.error.y.y.y"
                       # , "std.error.x.x.x.x", "std.error.y.y.y.y"
                       # , "std.error"
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
  "(A) Total species richness" = c("estimate.x", "std.error.x", "statistic.x", "p.value.x"),
  "(B) Char. species richness" = c("estimate.y", "std.error.y", "statistic.y", "p.value.y"),
  "(C) Total Hill-Shannon" = c("estimate.x.x", "std.error.x.x", "statistic.x.x", "p.value.x.x"),
  "(D) Char. Hill-Shannon" = c("estimate.y.y", "std.error.y.y", "statistic.y.y", "p.value.y.y"),
  "(E) FCSi" = c("estimate.x.x.x", "std.error.x.x.x", "statistic.x.x.x", "p.value.x.x.x"),
  "(F) Forb-grass ratio" = c("estimate.y.y.y", "std.error.y.y.y", "statistic.y.y.y", "p.value.y.y.y")
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
    decimals = 3,
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    rows = term == "No. of observations", # Apply only to the row where term == "No. of observations"
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
  ) %>% 
  cols_label(term = "")





gtsave(m_smry_gt,
       here(
         "outputs", "statistics", "vegetation", "model_summary_restfact_final_all_gt.docx"
       ))




### try out with p value stars

## --> doesn't work well yet, stars are not in right positions



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
      term == "rest.methmga" ~ "Rest. method = Management adaptation",
      term == "rest.methres" ~ "Rest. method = Regional seed mixture",
      term == "rest.methdih" ~ "Rest. method = Direct harvesting",
      term == "land.use.histgrassland" ~ "Previous landuse = Grassland",
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




## b - model < 20 years --------------------------------------------------------


# model_smry_joint_20y <- read_delim(
#   here("outputs", "statistics", "vegetation", "model_summary_restfact_20y_all.csv"),
#   col_names = TRUE, 
#   delim = ";"
# )



# Create the gt table
m_smry_gt <- model_smry_joint_20y %>%
  select(everything(), -effect, -component, -group) %>% 
  mutate(
    term = case_when(
      term == "rest.methmga" ~ "Rest. method = Management adaptation",
      term == "rest.methres" ~ "Rest. method = Regional seed mixture",
      term == "rest.methdih" ~ "Rest. method = Direct harvesting",
      term == "land.use.histgrassland" ~ "Previous land use = Grassland",
      term == "rest.age.std" ~ "Restoration Age",
      term == "sd__(Intercept) region" ~ "Region (random)",
      term == "sd__(Intercept) hydrology" ~ "Hydrology (random)",
      term == "R2_conditional" ~ "R2 conditional",
      term == "R2_marginal" ~ "R2 marginal",
      term == "nobs" ~ "No. of observations",
      TRUE ~ term)) %>% 
  filter(!term %in% c("Score_log", "Score_spherical", "AICc", "BIC")) %>% 
  gt()



# Create a list of columns for estimate and p.value
colnames(model_smry_joint_20y)
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
# change x and y manually according to sequence in creating model_smry_joint_20y
colnames(model_smry_joint_20y)
column_groups <- list(
  "(A) Total Species Richness" = c("estimate.x", "std.error.x", "statistic.x", "p.value.x"),
  "(B) Char. Species Richness" = c("estimate.y", "std.error.y", "statistic.y", "p.value.y"),
  "(C) Total Hill-Shannon" = c("estimate.x.x", "std.error.x.x", "statistic.x.x", "p.value.x.x"),
  "(D) Char. Hill-Shannon" = c("estimate.y.y", "std.error.y.y", "statistic.y.y", "p.value.y.y"),
  "(E) FCSi" = c("estimate.x.x.x", "std.error.x.x.x", "statistic.x.x.x", "p.value.x.x.x"),
  "(F) Forb-Grass Ratio" = c("estimate.y.y.y", "std.error.y.y.y", "statistic.y.y.y", "p.value.y.y.y")
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
col_names <- colnames(model_smry_joint_20y)
new_col_names <- col_names %>%
  gsub("^esti(.*)", "Est.", .) %>%
  gsub("^std(.*)", "S.E.", .) %>% 
  gsub("^stat(.*)", "z", .) %>% 
  gsub("^p.value(.*)", "p.value", .) 

names(new_col_names) <- col_names  # Map old to new names
m_smry_gt <- m_smry_gt %>%
  cols_label(!!!new_col_names[4:40])

m_smry_gt <- m_smry_gt %>% 
  # hide p.value column
  cols_hide(columns = matches("^p.value")) %>% 
  # round decimal places
  fmt_number(
    columns = where(is.numeric),
    decimals = 3,
  ) %>% 
  fmt_number(
    columns = where(is.numeric),
    rows = term == "No. of observations", # Apply only to the row where term == "No. of observations"
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
  ) %>% 
  cols_label(term = "")



gtsave(m_smry_gt,
       here(
         "outputs", "statistics", "vegetation", "model_summary_restfact_20y_all_gt.docx"
       ))














## end script



sum_gt2 <- sum_gt %>% 
  tab_style(
    style = cell_text(weight = "bold"),  # Set the text to bold
    locations = cells_body(
      columns = vars(starts_with("estimate")),  # Target all 'estimate' columns
      rows = p.value < 0.05                  # Apply condition for rows where p.value < 0.05
    )
  )


sum_gt <- m_sum2 %>% 
  select(everything(), -effect, -component, -group) %>% 
  gt() %>% 
  tab_spanner(
    label = "Total Species Richness",
    columns = c(estimate.x, std.error.x, statistic.x, p.value.x)
  ) %>% 
  tab_spanner(
    label = "Total Hill-Shannon",
    columns = c(estimate.y, std.error.y, statistic.y, p.value.y)
  ) %>%
  tab_spanner(
    label = "Total Diff q1-q0",
    columns = c(estimate.x.x, std.error.x.x, statistic.x.x, p.value.x.x)
  ) %>%
  tab_spanner(
    label = "Char. Species Richness",
    columns = c(estimate.y.y, std.error.y.y, statistic.y.y, p.value.y.y)
  ) %>%
  tab_spanner(
    label = "Char. Hill-Shannon",
    columns = c(estimate.x.x.x, std.error.x.x.x, statistic.x.x.x, p.value.x.x.x)
  ) %>%
  tab_spanner(
    label = "Char. Diff q1-q0",
    columns = c(estimate.y.y.y, std.error.y.y.y, statistic.y.y.y, p.value.y.y.y)
  ) %>%
  tab_spanner(
    label = "Forb-Grass Ratio",
    columns = c(estimate.x, std.error.x, statistic.x, p.value.x)
  ) %>%
  tab_spanner(
    label = "FCSi",
    columns = c(estimate.x, std.error.x, statistic.x, p.value.x)
  ) %>%
  tab_spanner(
    label = "Ellenberg L value",
    columns = c(estimate.x, std.error.x, statistic.x, p.value.x)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = "estimate.x",
      rows = p.value.x < 0.05
    )
  )


