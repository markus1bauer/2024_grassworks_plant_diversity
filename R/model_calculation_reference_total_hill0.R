#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: Total Species Richness (Hill 0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A PREPARATION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(glmmTMB)
library(performance) # visual check of model assumptions
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(rstatix)


### Start ###
rm(list = ls())


## load data -------------------------------------------------------------------

### site environment data ####

sites <- read_csv(
  here("data", "raw", "data_processed_environment_nms_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, hydrology, region) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         site.type = fct_relevel(site.type, "negative", "restored", "positive"),
  )


### diversity data ####

diversity <- read_csv(
  here("data", "raw", "data_processed_plants_site_diversity_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))


## set model data --------------------------------------------------------------

data_all <- sites %>% 
  left_join(diversity, by = "id.site") %>% 
  dplyr::select(id.site, tot.hill.0, hydrology, region, site.type)


rm(list = setdiff(ls(), c("data_all")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## a Missing values ------------------------------------------------------------
colSums(is.na(data_all)) 


## b Outliers ------------------------------------------------------------------


# Outliers: check with Cleveland dotplot
dotchart(data_all$tot.hill.0, ylab = "Order of the data")

# rstatix test for outliers
data_all %>% 
  dplyr::select(id.site, tot.hill.0) %>% 
  identify_outliers(tot.hill.0)


## c inspect categorical covariates -----------------------------------------

table(data_all$site.type)
table(data_all$hydrology)
table(data_all$region)
table(data_all$site.type, data_all$hydrology)
table(data_all$site.type, data_all$region)



## d Check collinearity --------------------------------------------------------

# between continuous covariates

# no numerical variable in model data



## e Relationships -------------------------------------------------------------

#' Plot response variable versus each covariate.
ggplot(data_all, aes(x = site.type, y = tot.hill.0)) +
   geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = region, y = tot.hill.0)) +
    geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = hydrology, y = tot.hill.0)) +
    geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")


## f distribution --------------------------------------------------------------

# --> count data: poisson or non-binomial


## g Interactions --------------------------------------------------------------

# --> no numerical covariates


# check categorical coviariates
library(MASS)

## site.type vs. X

# region
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$tot.hill.0)
int_model <- glm(tot.hill.0 ~ region * site.type, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(tot.hill.0 ~ region * site.type, data = data_all)
anova(int_model)

# hydrology
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$tot.hill.0)
int_model <- glm(tot.hill.0 ~ hydrology * site.type, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(tot.hill.0 ~ hydrology * site.type, data = data_all)
anova(int_model)

## region vs. X

# hydrology
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$tot.hill.0)
int_model <- glm(tot.hill.0 ~ hydrology * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(tot.hill.0 ~ hydrology * region, data = data_all)
anova(int_model)


detach(package:MASS)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - FULL MODEL ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------

data_model <- data_all

B1 <- glmmTMB(tot.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = poisson
)
check_overdispersion(B1)

# change to negative binomial family due to overdispersion
B1 <- glmmTMB(tot.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = nbinom2
)
check_overdispersion(B1)
summary(B1)
performance(B1)



## 2 Model validation full model -----------------------------------------------


#' Get residuals and fitted values of model
data_model$E_B1 <- resid(B1, type = "pearson")   
data_model$F_B1 <- fitted(B1)


### a Plot residuals vs fitted values ---------------------------------------------


ggplot(data = data_model, aes(x = F_B1, y = E_B1)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus site.type
ggplot(data = data_model, aes(x = site.type, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus region
ggplot(data = data_model, aes(x = region, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus hydrology
ggplot(data = data_model, aes(x = hydrology, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### c Plot residuals vs covariates not in the model ------------------------------

# no other covariates




## 3 Model fitting -------------------------------------------------------------

# --> no need to fit model, no interaction terms



## 4 Model validation drop model -----------------------------------------------

# --> not needed, because same model


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D - FINAL MODEL ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))


## 1 Summary -------------------------------------------------------------------

data_model_tothill0 <- data_all

restref_tothill0 <- glmmTMB(tot.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model_tothill0,
              family = nbinom2
)

summary(restref_tothill0)
performance(restref_tothill0)


## 2 Post-hoc test -------------------------------------------------------------

# estimated marginal means (EMMs)
# Tukey-adjusted pariwise comparisons
emm.site.type <- emmeans(restref_tothill0, "site.type")
summary(emm.site.type, infer = F, type = "response")
pairs(emm.site.type, adjust = "tukey")


## 3 Save final model ----------------------------------------------------------

save(restref_tothill0, data_model_tothill0,
     file = here("outputs", "models",
                 "model_reference_total_hill0.Rdata"))






## end script








