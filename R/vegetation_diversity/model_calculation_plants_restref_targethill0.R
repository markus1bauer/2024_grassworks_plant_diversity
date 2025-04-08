#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: Characteristic Species Richness (Hill 0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A PREPARATION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
# library(DHARMa)
library(performance) # visual check of model assumptions
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(rstatix)
# library(broom.mixed) # Tidy up the model summary
# library(scales) # label_number()
library(glmmTMB)
library(ggbeeswarm)


### Start ###
rm(list = ls())


## load data -------------------------------------------------------------------

### site environment data ####

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, hydrology, region,
    longitude, latitude
  ) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         site.type = fct_relevel(site.type, "negative", "restored", "positive"),
  )


### diversity data ####

diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))




## transform input data --------------------------------------------------------

# standardise explanatory variable (only numerical variables)

# --> no numerical explanatory variabls



## set model data --------------------------------------------------------------


# join diversity data
data_all <- sites %>%
  left_join(diversity, by = "id.site")


data_all <- data_all %>%
  dplyr::select(
    id.site,
    target.hill.0,
    hydrology,
    region,
    site.type
  )







rm(list = setdiff(ls(), c("data_all")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Protocol of data exploration (Steps 1-8)
# used from Zuur et al. (2010) Methods Ecol Evol 
#[DOI: 10.1111/2041-210X.12577](https://doi.org/10.1111/2041-210X.12577)


## a Missing values ------------------------------------------------------------
colSums(is.na(data_all)) 
# vis_dat(data_model)
# gg_miss_var(data_model)  

# --> no missing values



## b Outliers, zero-inflation, transformations? (Step 1, 3, 4) -----------------


# Outliers: check with Cleveland dotplot
dotchart(data_all$target.hill.0,
         ylab = "Order of the data")
# outlier are points far right or left in plot
# doesn't look like there are outliers


## wie umgehen mit ouliers?? Messfehler: unrealistische WErte
# sort(data_all$target.hill.0)


# another test for outliers
data_all %>% 
  select(id.site, target.hill.0) %>% 
  identify_outliers(target.hill.0)
# no extreme outliers


## c inspect categorical covariates -----------------------------------------

table(data_all$site.type)
#' Unbalanced...but enough observations per level.

table(data_all$hydrology)
#' Unbalanced...but enough observations per level.

table(data_all$region)
#' Balanced.

library(lattice)
#' Was each sitetype measured in every hydrology?
table(data_all$site.type, data_all$hydrology)
histogram( ~ site.type | hydrology, data_all)
#' Unbalanced, do we have enough observations per combination?
#' only 4 observation in dry-negative

#' Was each sitetype measured in every region?
table(data_all$site.type, data_all$region)
histogram( ~ site.type | region, data_all)
#' Unbalanced, but enough observations




## d Check collinearity part 1 (Step 5) ----------------------------------------

# between continuous covariates

# no numerical variable in model data --> no need to check



## e Relationships --------------------------------------------------------------

#' Plot response variable versus each covariate.


ggplot(data_all, aes(x = site.type, y = target.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Site type")
ggplot(data_all, aes(x = region, y = target.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = target.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# difference between site types
# difference between regions --> use as random factor
# difference between hydrology --> use as random factor

## f distribution --------------------------------------------------------------

# --> count data: poisson or non-binomial


# library(lattice)
# histogram(data_all$target.hill.0)
# # gamma distribution?
# 
# x <- data_all$target.hill.0
# 
# library(fitdistrplus)
# library(logspline)
# 
# descdist(x, discrete = FALSE)
# 
# fit.norm <- fitdist(x, "norm")
# fit.beta <- fitdist(x, "beta")
# fit.gamma <- fitdist(x, "gamma")
# 
# 
# plot(fit.norm)
# plot(fit.gamma)
# 
# # gamma looks okeyish
# 
# fit.norm$aic
# fit.gamma$aic
# # gamma has lowest AIC
# 
# detach(package:fitdistrplus)


## g Interactions --------------------------------------------------------------
# --> no numerical covariates


# check categorical coviariates

## site.type vs. X

# region
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$target.hill.0)
ggplot(data_all, aes(x = interaction(region, site.type), y = target.hill.0))+ 
  geom_boxplot()
# --> interaction between region and site.type

# hydrology
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.0)
ggplot(data_all, aes(x = interaction(hydrology, site.type), y = target.hill.0))+ 
  geom_boxplot()
# --> interaction between hydrology and site.type


## region vs. X

# hydrology
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.0)
ggplot(data_all, aes(x = interaction(hydrology, region), y = target.hill.0))+ 
  geom_boxplot()
# --> interaction between hydrology and region


# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

## site.type vs. X

# region
int_model <- glm(target.hill.0 ~ region * site.type, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ region * site.type, data = data_all)
anova(int_model)
# interaction

# hydrology
int_model <- glm(target.hill.0 ~ hydrology * site.type, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ hydrology * site.type, data = data_all)
anova(int_model)
# no interaction

## region vs. X

# hydrology
int_model <- glm(target.hill.0 ~ hydrology * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ hydrology * region, data = data_all)
anova(int_model)
# interaction

detach(package:MASS)



## h Spatial dependency --------------------------------------------------------------

library(rnaturalearth)


#' Get Germany map (medium resolution)
Germany <- ne_countries(country = "germany",
                        scale = "medium", 
                        returnclass = "sf")


ggplot(data = Germany) +
  geom_sf(fill = "transparent") +
  geom_point(data = data_all, 
             aes(x = longitude, 
                 y = latitude,),
             alpha = 0.3)  +
  theme_minimal() +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title = NULL)) +
  labs(title = "Sampling locations") 



#' Use this if you don't have online access:
ggplot(data = data_all) +
  theme(text = element_text(size=13)) +
  geom_point(aes(x = longitude, y = latitude, col = region), size = 0.75) +
  xlab("Longitude") + ylab("Latitude")


#' And zoom in per region:
ggplot(data = data_all) +
  theme(text = element_text(size=13)) +
  geom_point(aes(longitude, latitude, col = "black"), size = 0.75) +
  xlab("Longitude") + ylab("Latitude") +
  facet_wrap(~region, scale = "free", ncol = 2)
#' spatial variation per region is different, south less wide in latitude


## i conclusions  --------------------------------------------------------------

#' no missing values
#' no outliers
#' very little observations in negative-dry
#' interaction between:
#'  - site.type and region
#'  - region and hydrology
#' interaction unclear between:
#'  - site.type and hydrology


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 



### a Random structure ---------------------------------------------------------

R1 <- glmmTMB(target.hill.0 ~ 1 + (1|region), data = data_model, family = nbinom2)
R2 <- glmmTMB(target.hill.0 ~ 1 + (1|hydrology), data = data_model, family = nbinom2)
R3 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = nbinom2)
R4 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = nbinom2)
R5 <- glmmTMB(target.hill.0 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = nbinom2)
R6 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = nbinom2)
Rnull <- glm(target.hill.0 ~ 1, data = data_model, family = "poisson") # right family??


AIC(R1, R2, R3, R4, R5, R6, Rnull) %>% 
  arrange(AIC)
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure


# test if use of binomial-distribution is ok (due to overdispersion)
R3a <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = poisson)
check_overdispersion(R3a)
# significant result shows overdispersion
# --> using of negative binomial-distribution is neccessary




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D - FULL MODEL ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------


#  mu_ij = Intercept + site.type_ij + region_i + hydrology_i

# beyound optimal model
## consider interactions between explanatory variables
# -> only one explanatory variable, no interactions

data_model <- data_all

# check NA
colSums(is.na(data_model)) 
# no NA


# B1 crossed random intercept model with region and hydrology as random factors
B1 <- glmmTMB(target.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = poisson
)
check_overdispersion(B1)
# significant -> overdispersion

# change to non-binomial
B1 <- glmmTMB(target.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = nbinom2                      # Negative binomial family
)

summary(B1)
performance(B1)


#' Check overispersion for the NB GLMM.
check_overdispersion(B1)
#' No overdispersion detected.



## 2 Model validation full model -----------------------------------------------


#' As part of the model validation, we need to:
#'  -Plot residuals versus fitted values.
#'  -Plot residuals versus each covariate in the model.
#'  -Plot residuals versus each covariate NOT in the model.
#'  -Plot residuals versus time (if relevant).
#'  -Plot residuals versus spatial coordinates (if relevant).


#' Get residuals and fitted values of model
data_model$E_B1 <- resid(B1, type = "pearson")   #' Observed eveness minus fitted values.
data_model$F_B1 <- fitted(B1)  #' Fitted values contain the random effects.


### a Plot residuals vs fitted values ---------------------------------------------


ggplot(data = data_model, aes(x = F_B1, y = E_B1)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' Is there a pattern in here?
#' negative fitted values?
#' maybe smaller variance in residuals with larger values?



### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus site.type
ggplot(data = data_model, aes(x = site.type, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus region
ggplot(data = data_model, aes(x = region, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus hydrology
ggplot(data = data_model, aes(x = hydrology, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine


### c Plot residuals vs covariates not in the model ------------------------------

# no other covariates



### d check for adjustments of model -------------------------------------------

# no adjustments needed. model is fine




## 3 Model fitting -------------------------------------------------------------

# fit the model for reducing interactions,
# keep all main terms

# --> no need to fit model, no interaction terms



## 4 Model validation drop model -----------------------------------------------

# --> not needed, because same model



## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all")))


### a Summary ----

# the model we use is:
# target.hill.0 ~ site.type + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill0.Rdata"))

data_model_targethill0 <- data_all

restref_targethill0 <- glmmTMB(target.hill.0 ~ site.type + (1|region) + (1|hydrology), 
                    data = data_model_targethill0,
                    family = nbinom2                      # Negative binomial family
)

summary(restref_targethill0)
#                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        2.87094    0.14796  19.403   <2e-16 ***
# site.typerestored  0.60836    0.06408   9.493   <2e-16 ***
# site.typepositive  0.73994    0.07640   9.685   <2e-16 ***

# --> site.type is significant

performance(restref_targethill0)
#' ICC = 0.415
#' The covariates explain 32 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 60 % of the variation. (R2 cond.)

### summary statistics ###
data_model_targethill0 %>%
  group_by(site.type) %>%
  get_summary_stats(target.hill.0, type = "full")


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.site.type <- emmeans(restref_targethill0, "site.type")
summary(emm.site.type, infer = F, type = "response")
# site.type response   SE  df
# negative      17.7 2.61 Inf
# restored      32.4 4.49 Inf
# positive      37.0 5.33 Inf
pairs(emm.site.type, adjust = "tukey", 
      # type = "response"
)
# contrast            estimate     SE  df z.ratio p.value
# negative - restored   -0.608 0.0641 Inf  -9.493  <.0001
# negative - positive   -0.740 0.0764 Inf  -9.685  <.0001
# restored - positive   -0.132 0.0555 Inf  -2.370  0.0467
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 3 estimates

pairs(regrid(emm.site.type), adjust = "tukey") # regrid() for back-transformation from log-scale
## --> different p-values, because of calculating them after back-transformation

# --> restored are significantly different from positive reference! unlike in target.hill.0


### c report final model ----

# see skript show_table_model_results_restref.R



### d save final model ----


save(restref_targethill0, data_model_targethill0,
     file = here("outputs", "models", "vegetation",
                 "model_plants_restref_targethill0.Rdata"))






## end script
  
  
  
  
  
  
  
  
