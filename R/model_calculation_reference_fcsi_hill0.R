#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: FCSi (Index of Favourable Conservation Status)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATIION ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(tidyverse)
library(here)
library(performance)
library(emmeans)
library(rstatix)
library(glmmTMB)
library(ggbeeswarm)
library(mgcv)

### Start ###
rm(list = ls())

### Load data ###
data_all <- read_csv(
  here("data", "processed", "data_processed.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## a Missing values ------------------------------------------------------------
colSums(is.na(data_all)) 


## b Outliers ------------------------------------------------------------------


# Outliers: check with Cleveland dotplot
dotchart(data_all$fcsi.hill.0, ylab = "Order of the data")

# rstatix test for outliers
data_all %>% 
  dplyr::select(id.site, fcsi.hill.0) %>% 
  identify_outliers(fcsi.hill.0)


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
ggplot(data_all, aes(x = site.type, y = fcsi.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = region, y = fcsi.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = hydrology, y = fcsi.hill.0)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")



## f distribution --------------------------------------------------------------

library(lattice)
library(fitdistrplus)
library(logspline)

histogram(data_all$fcsi.hill.0)

x <- data_all$fcsi.hill.0
descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.gamma <- fitdist(x, "gamma")

plot(fit.norm)
plot(fit.gamma)

fit.norm$aic
fit.gamma$aic

detach(package:fitdistrplus)


## g Interactions --------------------------------------------------------------

# --> no numerical covariates


# check categorical coviariates
library(MASS)

## site.type vs. X

# region
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$fcsi.hill.0)
int_model <- glm(fcsi.hill.0 ~ region * site.type, data = data_all, family = gaussian)
anova(int_model)

# hydrology
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$fcsi.hill.0)
int_model <- glm(fcsi.hill.0 ~ hydrology * site.type, data = data_all, family = gaussian)
anova(int_model)


## region vs. X

# hydrology
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$fcsi.hill.0)
int_model <- glm(fcsi.hill.0 ~ hydrology * region, data = data_all, family = gaussian)
anova(int_model)


detach(package:MASS)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - FULL MODEL ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------

data_model <- data_all

B1 <- glmmTMB(fcsi.hill.0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = gaussian
)
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

data_model_fcsihill0 <- data_all

restref_fcsihill0 <- glmmTMB(fcsi.hill.0 ~ site.type + (1|region) + (1|hydrology), 
                            data = data_model_fcsihill0,
                            family = gaussian
)

summary(restref_fcsihill0)
performance(restref_fcsihill0)


## 2 Post-hoc test -------------------------------------------------------------

# estimated marginal means (EMMs)
# Tukey-adjusted pariwise comparisons
emm.site.type <- emmeans(restref_fcsihill0, "site.type")
summary(emm.site.type, infer = F, type = "response")
pairs(emm.site.type, adjust = "tukey")


## 3 Save final model ----------------------------------------------------------

save(restref_fcsihill0, data_model_fcsihill0,
     file = here("outputs", "models",
                 "model_reference_fcsi_hill0.Rdata"))






## end script








