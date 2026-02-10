#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites ####
# Response variable: Characteristic Hill-Shannon (q1)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke
# 2025



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


## 1 Missing values ------------------------------------------------------------

colSums(is.na(data_all)) 


## 2 Outliers ------------------------------------------------------------------

# Outliers: check with Cleveland dotplot
dotchart(data_all$target.hill.1, ylab = "Order of the data")

# rstatix test for outliers
data_all %>% 
  dplyr::select(id.site, target.hill.1) %>% 
  identify_outliers(target.hill.1)


## 3 inspect categorical covariates --------------------------------------------

table(data_all$site.type)
table(data_all$hydrology)
table(data_all$region)
table(data_all$site.type, data_all$hydrology)
table(data_all$site.type, data_all$region)


## 4 Check collinearity --------------------------------------------------------

# between continuous covariates

# no numerical variable in model data


## 5 Relationships -------------------------------------------------------------

#' Plot response variable versus each covariate.
ggplot(data_all, aes(x = site.type, y = target.hill.1)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = region, y = target.hill.1)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = hydrology, y = target.hill.1)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")


## 6 distribution --------------------------------------------------------------

library(lattice)
library(fitdistrplus)
library(logspline)

histogram(data_all$target.hill.1)

x <- data_all$target.hill.1
descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.gamma <- fitdist(x, "gamma")

plot(fit.norm)
plot(fit.gamma)

fit.norm$aic
fit.gamma$aic

detach(package:fitdistrplus)

# --> gamma distribution


## 7 Interactions --------------------------------------------------------------

# --> no numerical covariates


# check categorical coviariates
library(MASS)

## site.type vs. X

# region
int_model <- glm(target.hill.1 ~ region * site.type, data = data_all, family = Gamma(link="log"))
anova(int_model)
# no interaction

# hydrology
int_model <- glm(target.hill.1 ~ hydrology * site.type, data = data_all, family = Gamma(link="log"))
anova(int_model)
# no interaction


## region vs. X

# hydrology
int_model <- glm(target.hill.1 ~ hydrology * region, data = data_all, family = Gamma(link="log"))
anova(int_model)
# no interaction


detach(package:MASS)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - FULL MODEL ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rm(list = setdiff(ls(), c("data_all")))


## 1 Model formulation ---------------------------------------------------------

data_model <- data_all

B1 <- glmmTMB(target.hill.1 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = Gamma(link="log")
)
summary(B1)
performance(B1)


## 2 Model validation full model -----------------------------------------------

#' Get residuals and fitted values of model
data_model$E_B1 <- resid(B1, type = "pearson")   
data_model$F_B1 <- fitted(B1)


### a Plot residuals vs fitted values ------------------------------------------


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


### c Plot residuals vs covariates not in the model ----------------------------

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

data_model_targethill1 <- data_all

restref_targethill1 <- glmmTMB(target.hill.1 ~ site.type + (1|region) + (1|hydrology), 
                            data = data_model_targethill1,
                            family = Gamma(link="log")
)

summary(restref_targethill1)
performance(restref_targethill1)


## 2 Post-hoc test -------------------------------------------------------------

# estimated marginal means (EMMs)
# Tukey-adjusted pariwise comparisons
emm.site.type <- emmeans(restref_targethill1, "site.type")
summary(emm.site.type, infer = F, type = "response")
pairs(emm.site.type, adjust = "tukey")


## 3 Save final model ----------------------------------------------------------

# save(restref_targethill1, data_model_targethill1,
#      file = here("outputs", "models",
#                  "model_reference_target_hill1.Rdata"))


## end script