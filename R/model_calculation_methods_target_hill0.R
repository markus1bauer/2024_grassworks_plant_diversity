#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Response variable: Characteristic Species Richness (Hill 0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATIION ###############################################################
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
library(mgcv)



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
    id.site, site.type, hydrology, region, rest.meth, land.use.hist, rest.age
  ) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih"),
  )



### diversity data ####

diversity <- read_csv(
  here("data", "raw", "data_processed_plants_site_diversity_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))


## transform input data --------------------------------------------------------

# standardize numeric explanatory variables
data <- sites %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)), .names = "{col}.std"))



## set model data --------------------------------------------------------------


# join diversity data
data_all <- data %>%
  left_join(diversity, by = "id.site") %>% 
  dplyr::select(id.site, target.hill.0, hydrology, region,
                rest.meth, rest.age.std, land.use.hist)



rm(list = setdiff(ls(), c("data_all")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## a Missing values ------------------------------------------------------------
colSums(is.na(data_all)) 
# 66 NAs in rest.meth, rest.age and land.use.hist due to reference sites
# 13 other NAs in rest.age



## b Outliers ------------------------------------------------------------------

library(lattice)
Z <- data_all %>% 
  dplyr::select(where(is.numeric))


dotplot(as.matrix(Z), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data from text file")


# rstatix test for outliers
data_all %>% 
  select(id.site, target.hill.0) %>% 
  identify_outliers(target.hill.0)
data_all %>% 
  select(id.site, rest.age.std) %>% 
  identify_outliers(rest.age.std)



## c inspect categorical covariates -----------------------------------------

table(data_all$site.type)
table(data_all$hydrology)
table(data_all$region)
table(data_all$site.type, data_all$hydrology)
table(data_all$site.type, data_all$region)



## d Check collinearity part 1 ----------------------------------------

# between continuous covariates
# only one numerical variable in model data --> no need to check


# between continuous variables and factors

ggplot(data_all, aes(x = rest.meth, y = rest.age.std)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
data_all %>% 
  anova_test(rest.age.std ~ rest.meth)
# collinearity
ggplot(data_all, aes(x = land.use.hist, y = rest.age.std)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = hydrology, y = rest.age.std)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = region, y = rest.age.std)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")


## e Relationships --------------------------------------------------------------

#' Plot response variable versus each covariate.

ggplot(data_all, aes(x = rest.meth, y = target.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = land.use.hist, y = target.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent") 
ggplot(data_all, aes(x = rest.age.std, y = target.hill.0)) +
  geom_point() + geom_smooth(method = "glm") 
ggplot(data_all, aes(x = region, y = target.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent") 
ggplot(data_all, aes(x = hydrology, y = target.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")


## f distribution --------------------------------------------------------------

# --> count data: poisson or non-binomial


## g Interactions --------------------------------------------------------------

library(MASS)

## rest.meth vs. X

interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$land.use.hist,
                 response = data_all$target.hill.0)
# interaction
int_model <- glm(target.hill.0 ~ rest.meth * land.use.hist, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * land.use.hist, data = data_all)
anova(int_model)
# no interaction


ggplot(data_all, aes(rest.age.std, target.hill.0, color = rest.meth)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# no interaction
int_model <- glm(target.hill.0 ~ rest.meth * rest.age.std, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * rest.age.std, data = data_all)
anova(int_model)
# no interaction


interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$region,
                 response = data_all$target.hill.0)
# interaction
int_model <- glm(target.hill.0 ~ rest.meth * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * region, data = data_all)
anova(int_model)
# interaction

interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.0)
# interaction
int_model <- glm(target.hill.0 ~ rest.meth * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * hydrology, data = data_all)
anova(int_model)
# interaction


## land.use.hist vs. X

ggplot(data_all, aes(rest.age.std, target.hill.0, color = land.use.hist)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction
int_model <- glm(target.hill.0 ~ land.use.hist * rest.age.std, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * rest.age.std, data = data_all)
anova(int_model)
# no interaction 

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$region,
                 response = data_all$target.hill.0)
# no interaction
int_model <- glm(target.hill.0 ~ land.use.hist * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * region, data = data_all)
anova(int_model)
# no interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$hydrology,
                 response = data_all$target.hill.0)
# interaction
int_model <- glm(target.hill.0 ~ land.use.hist * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * hydrology, data = data_all)
anova(int_model)
# no interaction



## rest.age.std vs. X

ggplot(data_all, aes(rest.age.std, target.hill.0, color = region)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction
int_model <- glm(target.hill.0 ~ rest.age.std * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.age.std * region, data = data_all)
anova(int_model)
# interaction

ggplot(data_all, aes(rest.age.std, target.hill.0, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction
int_model <- glm(target.hill.0 ~ rest.age.std * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.age.std * hydrology, data = data_all)
anova(int_model)
# interaction


## region vs. hydrology

interaction.plot(x.factor = data_all$hydrology, trace.factor = data_all$region,
                 response = data_all$target.hill.0)
# interaction
int_model <- glm(target.hill.0 ~ region * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ region * hydrology, data = data_all)
anova(int_model)
# interaction


detach(package:MASS)


## i conclusions--------------------------------------------------------------

#' missing values in rest.age and land.use.hist
#' detected collinearity between restoration age and restoration method --> consider removing rest.age
#' interaction between:
#' - rest.meth and region, hydrology
#' - rest.age and region, hydrology
#' - region and hydrology
#' interaction unclear between: 
#' - rest.meth and land.use.hist
#' - land.use.hist and rest.age
#' no interaction between:
#' - rest.meth and rest.age
#' - land.use.hist and region, hydrology



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))


# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 
# 108 sites left
# NA in rest.age


R1 <- glmmTMB(target.hill.0 ~ 1 + (1|region), data = data_model, family = nbinom2)
R2 <- glmmTMB(target.hill.0 ~ 1 + (1|hydrology), data = data_model, family = nbinom2)
R3 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = nbinom2)
R4 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = nbinom2)
R5 <- glmmTMB(target.hill.0 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = nbinom2)
R6 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = nbinom2)
R7 <- glmmTMB(target.hill.0 ~ 1 + (rest.age.std|region) + (1|hydrology), data = data_model, family = nbinom2)
R8 <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (rest.age.std|hydrology), data = data_model, family = nbinom2)
# -> singularity
R9 <- glmmTMB(target.hill.0 ~ 1 + (rest.age.std|region) + (rest.age.std|hydrology), data = data_model, family = nbinom2)
# -> singularity
Rnull <- glmmTMB(target.hill.0 ~ 1, data = data_model, family = nbinom2) # right family??


AIC(R1, R2, R3, R4, R5, R6, R7, R8, R9, Rnull) %>% 
  arrange(AIC)
#       df       AIC
# R7     6  797.1703
# R4     4  798.7042
# R6     5  798.9827
# R3     4  799.5521
# R1     3  813.7030
# R2     3  845.8456
# Rnull  1 1104.0169

# R7, R4, R6 and R3 are very similar (< delta 2 AIC)
# use model R3 because it is easiest and test if random slope is needed
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(target.hill.0 ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later

# test if use of binomial-distribution is ok (due to overdispersion)
R3a <- glmmTMB(target.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = poisson)
check_overdispersion(R3a)
# significant result shows overdispersion
# --> using of negative binomial-distribution is neccessary

# # how strong is hydrology?
# glmm_hydr <- glmer.nb(target.hill.0 ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(target.hill.0 ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(target.hill.0 ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(target.hill.0 ~ rest.meth * hydrology + (1|region), data = data_model)
# summary(glmm_hydr)
# summary(glmm_test)
# summary(glmm_test_1)
# summary(glmm_test_2)
# # rest.meth is not significant anymore --> hydrology covers all effects


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D - FULL MODEL - <20 YEARS ##################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------


#' To account for collinearity between restoration age and restoration method we
#' would have to exclude one of the factors from analysis. First we remove colinerarity
#' and check for influence of restoration age in a reduced dataset with sites 
#' younger than 20 years.


## subset data to account for collinearity in age and restoration method
# exclude all sites > 20 years
data_model_20y <- data_all %>% 
  filter(rest.age <= 20)

# check NA
colSums(is.na(data_model_20y)) 
# no NA


## check collinearity
ggplot(data_model_20y, aes(x = rest.meth, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
# no collinearity

data_model_20y %>% 
  anova_test(rest.age ~ rest.meth)
# not significant, no collinearity



# beyound optimal model
#  mu_ij = Intercept + rest.meth_ij * land.use.hist_ij + rest.age_ij * land.use.hist_ij + region_i + hydrology_i

# consider interactions between explanatory variables

# B1_y crossed random intercept model with region and hydrology as random factors
yB1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = poisson
)
#' Check overispersion for poisson
check_overdispersion(yB1)
# overdispersion detected --> change to neg-binomial

yB1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = nbinom2                      # Negative binomial family
)

# yB1 <- glmer.nb(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
#                + (1|region) + (1|hydrology), data = data_model_20y)

summary(yB1)
# restoration age is not significant

#' Check overispersion for the NB GLMM.
check_overdispersion(yB1)
#' No overdispersion detected.


## 2 Model validation full model -----------------------------------------------


#' As part of the model validation, we need to:
#'  -Plot residuals versus fitted values.
#'  -Plot residuals versus each covariate in the model.
#'  -Plot residuals versus each covariate NOT in the model.
#'  -Plot residuals versus time (if relevant).
#'  -Plot residuals versus spatial coordinates (if relevant).


#' Get residuals and fitted values of model
data_model_20y$E_yB1 <- resid(yB1, type = "pearson")   #' Observed eveness minus fitted values.
data_model_20y$F_yB1 <- fitted(yB1)  #' Fitted values contain the random effects.


### a Plot residuals vs fitted values ---------------------------------------------


ggplot(data = data_model_20y, aes(x = F_yB1, y = E_yB1)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' Is there a pattern in here?
#' negative fitted values?
#' looks good



### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okeyish

#' Is there a pattern in here? To answer this question, fit a smoother on 
#' the residuals:
E_yB1 <- resid(yB1, type = "pearson")  #' Observed eveness minus fitted values.
F_yB1 <- fitted(yB1) #' Fitted values contain the random effects.

T_yB1 <- gam(E_yB1 ~ s(rest.age.std), data = data_model_20y)
summary(T_yB1)
# smoother not significant

#' Plot the smoother that was applied on the residuals.
plot(T_yB1, ylim = range(E_yB1))
points(x = data_model_20y$rest.age.std,
       y = E_yB1)
abline(h = 0, col = 2, lty = 2)
# no pattern, smoother falls together with line, no GAMM needed

#' Plot residuals versus rest.age for each region:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1, col = region)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' slopes are different for each region
#' It seems that we may need a random intercept and slope model.

#' Plot residuals versus rest.age for each hydrology:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1, col = hydrology)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' slopes are different for each hydrology
#' It seems that we may need a random intercept and slope model.


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus land.use.hist
ggplot(data = data_model_20y, aes(x = land.use.hist, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus region
ggplot(data = data_model_20y, aes(x = region, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus hydrology
ggplot(data = data_model_20y, aes(x = hydrology, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine


### c Plot residuals vs covariates not in the model ------------------------------
# --> no other covariates


### d Model validation with DHARMa ---------------------------------------------


# ### a Plot residuals vs fitted values
# simulation_output <- simulateResiduals(yB1, plot = TRUE)
# # quantile deviations detected
# 
# 
# ### b Plot residuals vs covariates in the model 
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$rest.meth) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$land.use.hist) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$rest.age.std) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$hydrology)
# # within-group deviations from uniformity significant
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$region)
# # within-group deviations from uniformity significant



### e check for adjustments of model -------------------------------------------

# check for need of random slope
yB1_sl1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = nbinom2                      # Negative binomial family
)
# --> warning: singular convergence
yB1_sl2 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = nbinom2                      # Negative binomial family
)
yB1_sl3 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (0 + rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = nbinom2                      # Negative binomial family
)
yB1_sl4 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (0 + rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = nbinom2                      # Negative binomial family
)

AIC(yB1, yB1_sl1, yB1_sl2, yB1_sl3, yB1_sl4) %>% arrange(AIC)
# random slope is not better --> stay with model 


## 3 Model fitting -------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1")))


# fit the model for reducing interactions,
# keep all main terms


### a backward selection (drop1) ----
# with drop1
drop1(yB1)
drop_model <- update(yB1, . ~ . - land.use.hist:rest.age.std)
drop1(drop_model)
# AIC is equal --> drop anyway because it simplifies model
drop_model <- update(drop_model, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)
# all interactions are gone. stop here


### b all subsets regression (MuMIn) ----
# with MuMIn package

# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# # (for some reason using glmmTMB ends up with not only removing interactions)
# yB1a <- glmer.nb(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist 
#                  + (1|region) + (1|hydrology), 
#                  data = data_model_20y)
# 
# full_model <- yB1a
# model_dredge <- dredge(full_model, beta = "none", evaluate = T, trace = 2,
#                        fixed = c("rest.meth", "land.use.hist", "rest.age.std"),
#                        # m.lim =c(0,5),
#                        rank = AIC) # when do you use AICc? Use AICc when n/k≤40 (n= sample size, k= no. of parameters)
# top_model <- get.models(model_dredge, subset = 1)[[1]]
# summary(top_model)
# 
# options(na.action = "na.omit") # set back to default


### c comparison fitting procedure ----

# AIC(full_model, drop_model, top_model)
# summary(top_model)
# # target.hill.0 ~ land.use.hist + rest.age.std + rest.meth +
# # (1 | region) + (1 | hydrology) + land.use.hist:rest.meth
# summary(drop_model)
# # target.hill.0 ~ rest.meth + land.use.hist + rest.age.std + (1 | region) + (1 | hydrology)
# 
# # MuMIn and drop1 end up different models, but they have the same AIC
# # --> use simpler model without interaction 

yB1_drop <- drop_model


## 4 Model validation drop model -----------------------------------------------


#' As part of the model validation, we need to:
#'  -Plot residuals versus fitted values.
#'  -Plot residuals versus each covariate in the model.
#'  -Plot residuals versus each covariate NOT in the model.
#'  -Plot residuals versus time (if relevant).
#'  -Plot residuals versus spatial coordinates (if relevant).


#' Get residuals and fitted values of model
data_model_20y$E_yB1_drop <- resid(yB1_drop, type = "pearson")   #' Observed eveness minus fitted values.
data_model_20y$F_yB1_drop <- fitted(yB1_drop)  #' Fitted values contain the random effects.


### a Plot residuals vs fitted values ---------------------------------------------


ggplot(data = data_model_20y, aes(x = F_yB1_drop, y = E_yB1_drop)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' Is there a pattern in here?
#' negative fitted values?
#' smaller variance in higher values?



### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okeyish

#' Is there a pattern in here? To answer this question, fit a smoother on 
#' the residuals:
E_yB1_drop <- resid(yB1_drop, type = "pearson")  #' Observed eveness minus fitted values.
F_yB1_drop <- fitted(yB1_drop) #' Fitted values contain the random effects.

T_yB1_drop <- gam(E_yB1_drop ~ s(rest.age.std), data = data_model_20y)
summary(T_yB1_drop)
# smoother not significant

#' Plot the smoother that was applied on the residuals.
plot(T_yB1_drop, ylim = range(E_yB1_drop))
points(x = data_model_20y$rest.age.std,
       y = E_yB1_drop)
abline(h = 0, col = 2, lty = 2)
# no pattern, smoother falls together with line, no GAMM needed

#' Plot residuals versus rest.age for each region:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop, col = region)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' slopes are different for each region
#' It seems that we may need a random intercept and slope model.

#' Plot residuals versus rest.age for each hydrology:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop, col = hydrology)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' slopes are different for each hydrology
#' It seems that we may need a random intercept and slope model.


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus land.use.hist
ggplot(data = data_model_20y, aes(x = land.use.hist, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus region
ggplot(data = data_model_20y, aes(x = region, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus hydrology
ggplot(data = data_model_20y, aes(x = hydrology, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine


### c Plot residuals vs covariates not in the model ------------------------------
# --> no other covariates


### d Model validation with DHARMa ---------------------------------------------


# ### a Plot residuals vs fitted values
# simulation_output <- simulateResiduals(yB1_drop, plot = TRUE)
# # quantile deviations detected
# 
# 
# ### b Plot residuals vs covariates in the model 
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$rest.meth) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$land.use.hist) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$rest.age.std) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$hydrology)
# # within-group deviations from uniformity significant
# plotResiduals(simulation_output$scaledResiduals, data_model_20y$region)
# # within-group deviations from uniformity significant


### e Check collinearity part 2 ------------------------------------------------

# Remove VIF > 3 or > 10
# Zuur et al. 2010 Methods Ecol Evol DOI: 10.1111/j.2041-210X.2009.00001.x

check_collinearity(yB1_drop)
# ok


# use glmer.nb() model
# car::vif(top_model)
# ok


### f check for adjustments of model -------------------------------------------

# check for need of random slope
yB1_drop_sl1 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = nbinom2                      # Negative binomial family
)
# --> warning: singular convergence
yB1_drop_sl2 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = nbinom2                      # Negative binomial family
)
# --> warning: singular convergence
yB1_drop_sl3 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (0 + rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = nbinom2                      # Negative binomial family
)
yB1_drop_sl4 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (0 + rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = nbinom2                      # Negative binomial family
)

AIC(yB1_drop, yB1_drop_sl1, yB1_drop_sl2, yB1_drop_sl3, yB1_drop_sl4) %>% arrange(AIC)
# random slope is not better --> stay with model 



## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1", "yB1_drop")))


### a Summary ----

# the model we use is:
# target.hill.0 ~ rest.meth + land.use.hist + rest.age.std + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models",  "vegetation" "model_plants_restfact_targethill0_yB1_final.Rdata"))

data_model_targethill0_20y <- data_all %>%
  filter(rest.age <= 20)


restfact_targethill0_20y <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                     + (1 |region) + (1 |hydrology),
                     data = data_model_targethill0_20y,
                     family = nbinom2                      # Negative binomial family
)
summary(restfact_targethill0_20y)
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             3.17691    0.16983  18.706  < 2e-16 ***
# rest.methmga            0.11232    0.14849   0.756  0.44940    
# rest.methres            0.27563    0.12119   2.274  0.02295 *  
# rest.methdih            0.40230    0.12642   3.182  0.00146 ** 
# land.use.histgrassland  0.10943    0.07056   1.551  0.12092    
# rest.age.std            0.01181    0.05412   0.218  0.82725  

# --> restoration age is not significant

performance(restfact_targethill0_20y)
#' ICC = 0.455
#' The covariates explain 15 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 54 % of the variation. (R2 cond.)


### b report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_targethill0_20y, effects = "fixed") %>%
#   select(term, estimate, std.error, statistic, p.value) %>%
#   mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
#   mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
#   mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                              .default = p.value))
# 
# # Display the summary as a table
# library(kableExtra)
# model_summary %>%
#   kbl(caption = "GLMM Model Results") %>%
#   kable_styling(full_width = F)
# 
# library(gt)
# model_summary %>%
#   gt() %>%
#   tab_header(title = "GLMM Model Results")
# 
# model_summary %>% 
#   write_csv(
#     here(
#       "outputs", "tables", "vegetation", "model_summary_restfact_targethill0_20y.csv"
#     )
#   )


### c save final model ----

save(restfact_targethill0_20y, data_model_targethill0_20y,
     file = here("outputs", "models",
                 "model_methods_target_hill0_20y.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E - FULL MODEL - NO RESTORATION AGE #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------


# include all restoration sites
data_model <- data_all %>%
  filter(!is.na(rest.meth))

# check NA
colSums(is.na(data_model)) 
# NA in no other variable than rest.age


#  mu_ij = Intercept + rest.meth_ij * land.use.hist_ij + region_i + hydrology_i

# beyound optimal model
## consider interactions between explanatory variables


# B1 crossed random intercept model with region and hydrology as random factors
B1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = poisson
)

#' Check overispersion for the poisson GLMM
check_overdispersion(B1)
# overdispersion detected --> change to neg-binomial


# B1 crossed random intercept model with region and hydrology as random factors
B1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = nbinom2                      # Negative binomial family
)


summary(B1)

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

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus land.use.hist
ggplot(data = data_model, aes(x = land.use.hist, y = E_B1)) +
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

#' Plot residuals versus restoration age
ggplot(data = data_model, aes(x = rest.age.std, y = E_B1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' stronger variance in lower values!


### d Model validation with DHARMa ---------------------------------------------


# ### a Plot residuals vs fitted values
# simulation_output <- simulateResiduals(B1, plot = TRUE)
# # quantile deviations detected
# 
# 
# ### b Plot residuals vs covariates in the model
# plotResiduals(simulation_output$scaledResiduals, data_model$rest.meth) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model$land.use.hist) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model$hydrology)
# # within-group deviations from uniformity significant
# plotResiduals(simulation_output$scaledResiduals, data_model$region)
# # within-group deviations from uniformity significant
# 
# ### c Plot residuals vs covariates not in the model
# plotResiduals(simulation_output$scaledResiduals, data_model$rest.age.std) # ok



### e check for adjustments of model -------------------------------------------

# no adjustments needed. model is fine




## 3 Model fitting -------------------------------------------------------------

# fit the model for reducing interactions,
# keep all main terms


### a backward selection (drop1) ----
# with drop1
drop1(B1)
drop_model <- update(B1, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)
# all interactions are gone. stop here


### b all subsets regression (MuMIn) ----
# with MuMIn package

# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# B1a <- glmer.nb(target.hill.0 ~ rest.meth * land.use.hist 
#                 + (1|region) + (1|hydrology), 
#                 data = data_model)
# 
# full_model <- B1a
# model_dredge <- dredge(full_model, beta = "none", evaluate = T, trace = 2,
#                        fixed = c("rest.meth", "land.use.hist"),
#                        # m.lim =c(0,5),
#                        rank = AIC) # when do you use AICc? Use AICc when n/k≤40 (n= sample size, k= no. of parameters)
# top_model <- get.models(model_dredge, subset = 1)[[1]]
# summary(top_model)
# 
# options(na.action = "na.omit") # set back to default


### c comparison fitting procedure ----

# AIC(full_model, drop_model, top_model)
# # df      AIC
# # full_model 11 869.7192
# # drop_model  8 867.4274
# # top_model   8 867.4711
# summary(top_model)
# # target.hill.0 ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# summary(drop_model)
# # target.hill.0 ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# 
# # MuMIn and drop1 end up with the same model -> super

B1_drop <- drop_model


## 4 Model validation drop model -----------------------------------------------


#' As part of the model validation, we need to:
#'  -Plot residuals versus fitted values.
#'  -Plot residuals versus each covariate in the model.
#'  -Plot residuals versus each covariate NOT in the model.
#'  -Plot residuals versus time (if relevant).
#'  -Plot residuals versus spatial coordinates (if relevant).


#' Get residuals and fitted values of model
data_model$E_B1_drop <- resid(B1_drop, type = "pearson")   #' Observed eveness minus fitted values.
data_model$F_B1_drop <- fitted(B1_drop)  #' Fitted values contain the random effects.


### a Plot residuals vs fitted values ---------------------------------------------


ggplot(data = data_model, aes(x = F_B1_drop, y = E_B1_drop)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' Is there a pattern in here?
#' negative fitted values?
#' maybe smaller variance in residuals with larger values?



### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus land.use.hist
ggplot(data = data_model, aes(x = land.use.hist, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus region
ggplot(data = data_model, aes(x = region, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

#' Plot the residuals versus hydrology
ggplot(data = data_model, aes(x = hydrology, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine


### c Plot residuals vs covariates not in the model ------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model, aes(x = rest.age.std, y = E_B1_drop)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' stronger variance in lower values!


### d Model validation with DHARMa ---------------------------------------------


# ### a Plot residuals vs fitted values
# simulation_output <- simulateResiduals(B1_drop, plot = TRUE)
# # ok
# 
# 
# ### b Plot residuals vs covariates in the model
# plotResiduals(simulation_output$scaledResiduals, data_model$rest.meth) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model$land.use.hist) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model$hydrology)
# # within-group deviations from uniformity significant
# plotResiduals(simulation_output$scaledResiduals, data_model$region)
# # within-group deviations from uniformity significant
# 
# ### c Plot residuals vs covariates not in the model
# plotResiduals(simulation_output$scaledResiduals, data_model$rest.age.std) # ok



### e Check collinearity part 2 ------------------------------------------------

# Remove VIF > 3 or > 10
# Zuur et al. 2010 Methods Ecol Evol DOI: 10.1111/j.2041-210X.2009.00001.x

check_collinearity(B1_drop)
# ok

# use glmer.nb() model
# car::vif(top_model)
# ok



### f check for adjustments of model -------------------------------------------

# no adjustments needed. model is fine




## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model", "B1", "B1_drop")))


### a Summary ----

# the model we use is:
# target.hill.0 ~ rest.meth + land.use.hist + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0_B1_final.Rdata"))

data_model_targethill0 <- data_all %>%
  filter(!is.na(rest.meth))

restfact_targethill0 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist
                    + (1 |region) + (1 |hydrology),
                    data = data_model_targethill0,
                    family = nbinom2                      # Negative binomial family
)
summary(restfact_targethill0)
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             3.22505    0.15154  21.282  < 2e-16 ***
# rest.methmga            0.10544    0.09686   1.089  0.27633    
# rest.methres            0.22939    0.07202   3.185  0.00145 ** 
# rest.methdih            0.32295    0.07568   4.267 1.98e-05 ***
# land.use.histgrassland  0.10190    0.05359   1.901  0.05724 . 

# --> rest.meth is significant, land.use.hist on the edge

performance(restfact_targethill0)
#' ICC = 0.517
#' The covariates explain 13 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 58 % of the variation. (R2 cond.)


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.rest.meth <- emmeans(restfact_targethill0, ~ rest.meth)
emm.rest.meth
summary(emm.rest.meth, infer = F, type = "response")
pairs(emm.rest.meth, adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast  estimate     SE  df z.ratio p.value
# cus - mga  -0.1054 0.0969 Inf  -1.089  0.6965
# cus - res  -0.2294 0.0720 Inf  -3.185  0.0079
# cus - dih  -0.3229 0.0757 Inf  -4.267  0.0001
# mga - res  -0.1239 0.0841 Inf  -1.473  0.4539
# mga - dih  -0.2175 0.0805 Inf  -2.702  0.0348
# res - dih  -0.0936 0.0624 Inf  -1.499  0.4381
# 
# Results are averaged over the levels of: land.use.hist 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

pairs(regrid(emm.rest.meth), adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast  estimate   SE  df z.ratio p.value
# cus - mga    -2.94 2.73 Inf  -1.078  0.7033
# cus - res    -6.82 2.29 Inf  -2.980  0.0153
# cus - dih   -10.09 2.64 Inf  -3.819  0.0008
# mga - res    -3.88 2.67 Inf  -1.456  0.4646
# mga - dih    -7.15 2.76 Inf  -2.591  0.0471
# res - dih    -3.27 2.21 Inf  -1.478  0.4511
# 
# Results are averaged over the levels of: land.use.hist 
# P value adjustment: tukey method for comparing a family of 4 estimates  

## --> different p-values, because of calculating them after back-transformation


# --> res is significantly different from cus (in tot.hill.0 it is on the edge)
# --> res is not significantly different from dih (in tot.hill.0 it is on the edge)


### c report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_targethill0, effects = "fixed") %>%
#   select(term, estimate, std.error, statistic, p.value) %>%
#   mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
#   mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
#   mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
#                              .default = p.value))
# 
# # Display the summary as a table
# library(kableExtra)
# model_summary %>%
#   kbl(caption = "GLMM Model Results") %>%
#   kable_styling(full_width = F)
# 
# library(gt)
# model_summary %>%
#   gt() %>%
#   tab_header(title = "GLMM Model Results")
# 
# model_summary %>% 
#   write_csv(
#     here(
#       "outputs", "tables", "vegetation", "model_summary_restfact_targethill0.csv"
#     )
#   )

### summary statistics ###
data_model_target %>%
  group_by(rest.meth) %>%
  get_summary_stats(target.hill.0, type = "full")


### d save final model ----


save(restfact_targethill0, data_model_targethill0,
     file = here("outputs", "models",
                 "model_methods_target_hill0_full.Rdata"))





## end script
  
  
  
  
  
  
  
  
  
  
  
  glmm_rest_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_model)
summary(glmm_rest_age)

glmm_rest_5 <- glmer.nb(target.hill.0 ~
                          rest.meth
                        * rest.age
                        + (1|region) + (1|hydrology)
                        , data = data_model
)
summary(glmm_rest_5)


range(data_model$rest.age)

# Vorhersagen und Visualisierung
newdat1 <- expand.grid(rest.age = seq(0, 40, by = 0.01),
                       region = unique(data_model$region),
                       hydrology = unique(data_model$hydrology),
                       rest.meth = unique(data_model$rest.meth)
)
# expand.grid kombiniert alle angegebenen Werte einer Variable mit allen
# Werten der anderen. Hier also 101 Werte für soil-prep x 3 Werte region + 3 Werte
# hydrology = 909 Zeilen in der Tabelle


# Vorhersagen mit random effect
newdat1$pred_hill0_glmm <- predict(glmm_rest_5, newdata = newdat1,
                                   type = "response")

# Vorhersagen nur mit fixed effect
newdat1$pred_hill0_glmm_fixed <- predict(glmm_rest_5, newdata = newdat1, type = "response", re.form = NA)

# Finale Grafik ---
ggplot(data_model, aes(rest.age, target.hill.0, color = region)) +
  geom_point() +
  theme_bw() +
  # facet_wrap(~hydrology) +
  facet_grid(hydrology ~ rest.meth) +
  geom_line(aes(rest.age, pred_hill0_glmm), data = newdat1) +
  geom_line(aes(rest.age, pred_hill0_glmm_fixed),
            linewidth = 1, linetype = 2, data = newdat1, color = "black")
ggsave("outputs/figures/plants_species_diversity/model_tothill0_restmeth-restage.jpg")

# Die farbigen Linien beinhalten die random intercepts für Region und Hydrologie
# Die schwarze Linie ist zeigt nur den fixed effect der soil-preparation





glmm_rest_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_model)
summary(glmm_rest_age)

range(data_model$rest.age)

# Vorhersagen und Visualisierung
newdat1 <- expand.grid(rest.age = seq(0, 40, by = 0.01),
                       region = unique(data_model$region),
                       hydrology = unique(data_model$hydrology))
# expand.grid kombiniert alle angegebenen Werte einer Variable mit allen
# Werten der anderen. Hier also 101 Werte für soil-prep x 3 Werte region + 3 Werte
# hydrology = 909 Zeilen in der Tabelle


# Vorhersagen mit random effect
newdat1$pred_hill0_glmm <- predict(glmm_rest_age, newdata = newdat1,
                                   type = "response")

# Vorhersagen nur mit fixed effect
newdat1$pred_hill0_glmm_fixed <- predict(glmm_rest_age, newdata = newdat1, type = "response", re.form = NA)

# Finale Grafik ---
ggplot(data_model, aes(rest.age, target.hill.0, color = region)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~hydrology) +
  geom_line(aes(rest.age, pred_hill0_glmm), data = newdat1) +
  geom_line(aes(rest.age, pred_hill0_glmm_fixed),
            linewidth = 1, linetype = 2, data = newdat1, color = "black")

# Die farbigen Linien beinhalten die random intercepts für Region und Hydrologie
# Die schwarze Linie ist zeigt nur den fixed effect der soil-preparation
































# Code snippets ----

## GAMM ----

# Load mgcv package
library(mgcv)

# Convert the GLMM to a GAMM
B1_gamm <- gamm(
  target.hill.0 ~ rest.meth + land.use.hist + 
    s(rest.age.std, by = land.use.hist),   # Smooth for rest.age.std by land.use.hist
  random = list(region = ~ 1, hydrology = ~ 1),        # Random effects
  # family = nb(),                                       # Negative binomial family
  data = data_model
)
summary(B1_gamm$gam)
summary(B1_gamm$lme)

# get normalized residuals
E4 <- resid(B1_gamm$lme, type = "normalized")

# plot residuals vs covariates in the model
plot(x = data_model$rest.age.std, y = E4)
abline(h = 0, lty = 2)

boxplot(E4 ~ data_model$region, data = data_model)
abline(h = 0, lty = 2)

boxplot(E4 ~ data_model$hydrology, data = data_model)
abline(h = 0, lty = 2)



## MUMIn package

## all subsets regression
# with MuMIn package
options(na.action = "na.fail") # Required for dredge to run

full_model <- B1
model_dredge <- dredge(full_model, beta = "none", evaluate = T, trace = 2,
                       fixed = c("rest.meth", "land.use.hist", "rest.age.std"),
                       # m.lim =c(0,5),
                       rank = AICc) # when do you use AICc? Use AICc when n/k≤40 (n= sample size, k= no. of parameters)
top_model <- get.models(model_dredge, subset = 1)[[1]]
summary(top_model)


# get top models
top_most_models <- get.models(model_dredge, subset = delta < 2)
top_most_models
# second best model (by AICc) is target.hill.0 ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# AIC 821.2433
# which is very very close to top model!
# Rule of thumb often used : models with AIC.delta <=2 are equally supported by the data

# get second best model
B1_top2 <- get.models(model_dredge, subset = 2)[[1]]


# Average model
av_model <- model.avg(model_dredge, subset = delta < 2)
summary(av_model)

options(na.action = "na.omit") # set back to default


B1a <- glm(target.hill.0 ~ rest.meth, data = data_model, family = "poisson")

AICc(B1a, B1_top)







## performance package

check_model(top_model)



## lmerTest

## backward stepwise regression ####
bw_model <- step(full_model, direction = "backward")






## subset models ---------------------------------------------------------------

glmm_dih_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_dih)
glmm_mga_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_mga)
glmm_res_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_res)
glmm_cus_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_cus)
summary(glmm_dih_age)
summary(glmm_mga_age)
summary(glmm_res_age)
summary(glmm_cus_age)

plot(data_all$target.hill.0 ~ data_all$rest.age)

glmm_dry_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_dry)
glmm_fresh_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_fresh)
glmm_moist_age <- glmer.nb(target.hill.0 ~ rest.age + (1|region), data = data_moist)
summary(glmm_dry_age)
summary(glmm_fresh_age)
summary(glmm_moist_age)


