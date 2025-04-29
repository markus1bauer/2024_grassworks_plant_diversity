#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: Characteristic Hill-Shannon (q1)
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
  here("data", "raw", "data_processed_environment_nms_20250306.csv"),
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
  here("data", "raw", "data_processed_plants_site_diversity_20250306.csv"),
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
    target.hill.1,
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
dotchart(data_all$target.hill.1,
         ylab = "Order of the data")
# outlier are points far right or left in plot
# doesn't look like there are outliers


## wie umgehen mit ouliers?? Messfehler: unrealistische WErte
# sort(data_all$target.hill.1)


# another test for outliers
data_all %>% 
  select(id.site, target.hill.1) %>% 
  identify_outliers(target.hill.1)
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


ggplot(data_all, aes(x = site.type, y = target.hill.1)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Site type")
ggplot(data_all, aes(x = region, y = target.hill.1)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = target.hill.1)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# difference between site types
# difference between regions --> use as random factor
# difference between hydrology --> use as random factor


## f distribution --------------------------------------------------------------

library(lattice)
histogram(data_all$target.hill.1)
# gamma distribution?

x <- data_all$target.hill.1

library(fitdistrplus)
library(logspline)

descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.beta <- fitdist(x, "beta")
fit.gamma <- fitdist(x, "gamma")


plot(fit.norm)
plot(fit.gamma)

# gamma looks okeyish

fit.norm$aic
fit.gamma$aic
# gamma has lowest AIC

detach(package:fitdistrplus)


## g Interactions --------------------------------------------------------------

# --> no numerical covariates


# check categorical coviariates

## site.type vs. X

# region
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$target.hill.1)
ggplot(data_all, aes(x = interaction(region, site.type), y = target.hill.1))+ 
  geom_boxplot()
# --> interaction between region and site.type

# hydrology
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.1)
ggplot(data_all, aes(x = interaction(hydrology, site.type), y = target.hill.1))+ 
  geom_boxplot()
# --> interaction between hydrology and site.type


## region vs. X

# hydrology
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.1)
ggplot(data_all, aes(x = interaction(hydrology, region), y = target.hill.1))+ 
  geom_boxplot()
# --> interaction between hydrology and region


# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

## site.type vs. X

# region
int_model <- glm(target.hill.1 ~ region * site.type, data = data_all, family = Gamma(link="log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

# hydrology
int_model <- glm(target.hill.1 ~ hydrology * site.type, data = data_all, family = Gamma(link="log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

## region vs. X

# hydrology
int_model <- glm(target.hill.1 ~ hydrology * region, data = data_all, family = Gamma(link="log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

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
#' interaction unclear between:
#'  - site.type and hydrology
#'  - site.type and region
#'  - region and hydrology



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 



### a Random structure ---------------------------------------------------------

R1 <- glmmTMB(target.hill.1 ~ 1 + (1|region), data = data_model, family = Gamma(link="log"))
R2 <- glmmTMB(target.hill.1 ~ 1 + (1|hydrology), data = data_model, family = Gamma(link="log"))
R3 <- glmmTMB(target.hill.1 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = Gamma(link="log"))
R4 <- glmmTMB(target.hill.1 ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R5 <- glmmTMB(target.hill.1 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R6 <- glmmTMB(target.hill.1 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
Rnull <- glm(target.hill.1 ~ 1, data = data_model, family = Gamma(link="log")) # right family??


AIC(R1, R2, R3, R4, R5, R6, Rnull) %>% 
  arrange(AIC)
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure


check_overdispersion(R3)
# no overdispersion



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
B1 <- glmmTMB(target.hill.1 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = Gamma(link="log")
)
check_overdispersion(B1)
# no overdispersion

summary(B1)
performance(B1)



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
# target.hill.1 ~ site.type + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill1.Rdata"))

data_model_targethill1 <- data_all

restref_targethill1 <- glmmTMB(target.hill.1 ~ site.type + (1|region) + (1|hydrology), 
                            data = data_model_targethill1,
                            family = Gamma(link="log")
)

summary(restref_targethill1)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        1.75351    0.18932   9.262  < 2e-16 ***
#   site.typerestored  0.68769    0.08397   8.189 2.62e-16 ***
#   site.typepositive  0.94694    0.10469   9.045  < 2e-16 ***

# --> site.type is significant

performance(restref_targethill1)
#' ICC = 0.357
#' The covariates explain 26 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 53 % of the variation. (R2 cond.)

### summary statistics ###
data_model_targethill1 %>%
  group_by(site.type) %>%
  get_summary_stats(target.hill.1, type = "full")


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.site.type <- emmeans(restref_targethill1, "site.type")
summary(emm.site.type, infer = F, type = "response")
# site.type response   SE  df
# negative      5.77 1.09 Inf
# restored     11.49 2.04 Inf
# positive     14.89 2.80 Inf

pairs(emm.site.type, adjust = "tukey", 
      # type = "response"
)
# contrast            estimate     SE  df z.ratio p.value
# negative - restored   -0.688 0.0840 Inf  -8.189  <.0001
# negative - positive   -0.947 0.1047 Inf  -9.045  <.0001
# restored - positive   -0.259 0.0833 Inf  -3.112  0.0053
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 3 estimates  

pairs(regrid(emm.site.type), adjust = "tukey") # regrid() for back-transformation from log-scale
## --> different p-values, because of calculating them after back-transformation

# --> restored are sign different from pos and neg ref


### c report final model ----

# see skript show_table_model_results_restref.R



### d save final model ----


save(restref_targethill1, data_model_targethill1,
     file = here("outputs", "models",
                 "model_reference_target_hill1.Rdata"))






## end script








