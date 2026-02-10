#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors ####
# Response variable: Characteristic Species Richness (Hill 0)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke
# 2025



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATIION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(tidyverse)
library(here)
library(glmmTMB)
library(performance) # visual check of model assumptions
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(rstatix)
library(mgcv)

### Start ###
rm(list = ls())

### Load data ###
data_all <- read_csv(
  here("data", "processed", "data_processed.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION #########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Missing values ------------------------------------------------------------

colSums(is.na(data_all)) 


## 2 Outliers ------------------------------------------------------------------

# Outliers: check with Cleveland dotplot
dotchart(data_all$target.hill.0, ylab = "Order of the data")
dotchart(data_all$rest.age.std, ylab = "Order of the data")


# rstatix test for outliers
data_all %>% 
  select(id.site, target.hill.0) %>% 
  identify_outliers(target.hill.0)
data_all %>% 
  select(id.site, rest.age.std) %>% 
  identify_outliers(rest.age.std)


## 3 inspect categorical covariates --------------------------------------------

table(data_all$rest.meth)
table(data_all$land.use.hist)
table(data_all$region)
table(data_all$hydrology)
table(data_all$rest.meth, data_all$land.use.hist)
table(data_all$rest.meth, data_all$region)
table(data_all$rest.meth, data_all$hydrology)
table(data_all$land.use.hist, data_all$region)
table(data_all$land.use.hist, data_all$hydrology)
table(data_all$region, data_all$hydrology)


## 4 Check collinearity part 1 -------------------------------------------------

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


## 5 Relationships -------------------------------------------------------------

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


## 6 distribution --------------------------------------------------------------

# --> count data: poisson or non-binomial


## 7 Interactions --------------------------------------------------------------

library(MASS)


## rest.meth vs. X

# land.use.hist
int_model <- glm(target.hill.0 ~ rest.meth * land.use.hist, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * land.use.hist, data = data_all)
anova(int_model)
# no interaction

# rest.age
int_model <- glm(target.hill.0 ~ rest.meth * rest.age.std, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * rest.age.std, data = data_all)
anova(int_model)
# no interaction

# region
int_model <- glm(target.hill.0 ~ rest.meth * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * region, data = data_all)
anova(int_model)
# interaction

# hydrology
int_model <- glm(target.hill.0 ~ rest.meth * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.meth * hydrology, data = data_all)
anova(int_model)
# interaction


## land.use.hist vs. X

# rest.age
int_model <- glm(target.hill.0 ~ land.use.hist * rest.age.std, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * rest.age.std, data = data_all)
anova(int_model)
# no interaction 

# region
int_model <- glm(target.hill.0 ~ land.use.hist * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * region, data = data_all)
anova(int_model)
# no interaction

# hydrology
int_model <- glm(target.hill.0 ~ land.use.hist * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ land.use.hist * hydrology, data = data_all)
anova(int_model)
# no interaction


## rest.age.std vs. X

# region
int_model <- glm(target.hill.0 ~ rest.age.std * region, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.age.std * region, data = data_all)
anova(int_model)
# interaction

# hydrology
int_model <- glm(target.hill.0 ~ rest.age.std * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ rest.age.std * hydrology, data = data_all)
anova(int_model)
# interaction


## region vs. X

# hydrology
int_model <- glm(target.hill.0 ~ region * hydrology, data = data_all, family = "poisson")
check_overdispersion(int_model)
int_model <- glm.nb(target.hill.0 ~ region * hydrology, data = data_all)
anova(int_model)
# interaction


detach(package:MASS)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - MODEL <20 YEARS (APPENDIX) ##############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rm(list = setdiff(ls(), c("data_all")))


## 0 Remove colinearity --------------------------------------------------------

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
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
# no collinearity

data_model_20y %>% 
  anova_test(rest.age ~ rest.meth)
# not significant, no collinearity


## 1 Model formulation ---------------------------------------------------------

yB1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = poisson
)
check_overdispersion(yB1)

# change to negative binomial family due to overdispersion
yB1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = nbinom2
)
summary(yB1)
performance(yB1)


## 2 Model validation full model -----------------------------------------------

#' Get residuals and fitted values of model
data_model_20y$E_yB1 <- resid(yB1, type = "pearson")
data_model_20y$F_yB1 <- fitted(yB1)


### a Plot residuals vs fitted values ---------------------------------------------

ggplot(data = data_model_20y, aes(x = F_yB1, y = E_yB1)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' fit a smoother on the residuals:
E_yB1 <- resid(yB1, type = "pearson")
F_yB1 <- fitted(yB1)
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
#' check need for random intercept and slope model

#' Plot residuals versus rest.age for each hydrology:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1, col = hydrology)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' check need for random intercept and slope model


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus land.use.hist
ggplot(data = data_model_20y, aes(x = land.use.hist, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus region
ggplot(data = data_model_20y, aes(x = region, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus hydrology
ggplot(data = data_model_20y, aes(x = hydrology, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### c Plot residuals vs covariates not in the model ------------------------------
# --> no other covariates


### d check for adjustments of model -------------------------------------------

# check for need of random slope
yB1_sl1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = nbinom2    
)
# --> warning: singular convergence
yB1_sl2 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = nbinom2    
)
yB1_sl3 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (0 + rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = nbinom2    
)
yB1_sl4 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (0 + rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = nbinom2    
)

AIC(yB1, yB1_sl1, yB1_sl2, yB1_sl3, yB1_sl4) %>% arrange(AIC)
# random slope is not better --> stay with model 


## 3 Model fitting -------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1")))


# fit the model for reducing interactions,
# keep all main terms


### a backward selection (drop1) ----
drop1(yB1)
drop_model <- update(yB1, . ~ . - land.use.hist:rest.age.std)
drop1(drop_model)
# AIC is equal --> drop because it simplifies model
drop_model <- update(drop_model, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)

### b drop model ----
yB1_drop <- drop_model


## 4 Model validation drop model -----------------------------------------------

#' Get residuals and fitted values of model
data_model_20y$E_yB1_drop <- resid(yB1_drop, type = "pearson")
data_model_20y$F_yB1_drop <- fitted(yB1_drop)


### a Plot residuals vs fitted values ---------------------------------------------

ggplot(data = data_model_20y, aes(x = F_yB1_drop, y = E_yB1_drop)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okeyish

#' fit a smoother on the residuals:
E_yB1_drop <- resid(yB1_drop, type = "pearson")
F_yB1_drop <- fitted(yB1_drop)
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
#' check need for random intercept and slope model

#' Plot residuals versus rest.age for each hydrology:
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop, col = hydrology)) +
  geom_point(shape = 1, size = 1) +
  labs(x = "values", y = "Residuals") + 
  geom_smooth(method = "glm", se = FALSE)
#' check need for random intercept and slope model


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus land.use.hist
ggplot(data = data_model_20y, aes(x = land.use.hist, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus region
ggplot(data = data_model_20y, aes(x = region, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus hydrology
ggplot(data = data_model_20y, aes(x = hydrology, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### c Plot residuals vs covariates not in the model ------------------------------
# --> no other covariates


### d Check collinearity part 2 ------------------------------------------------

# Remove VIF > 3 or > 10
# Zuur et al. 2010 Methods Ecol Evol DOI: 10.1111/j.2041-210X.2009.00001.x

check_collinearity(yB1_drop)
# ok


### e check for adjustments of model -------------------------------------------

# check for need of random slope
yB1_drop_sl1 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = nbinom2    
)
# --> warning: singular convergence
yB1_drop_sl2 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = nbinom2    
)
# --> warning: singular convergence
yB1_drop_sl3 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (0 + rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = nbinom2    
)
yB1_drop_sl4 <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (0 + rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = nbinom2    
)

AIC(yB1_drop, yB1_drop_sl1, yB1_drop_sl2, yB1_drop_sl3, yB1_drop_sl4) %>% arrange(AIC)
# random slope is not better --> stay with model 


## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all")))

### a Summary ----

data_model_targethill0_20y <- data_all %>%
  filter(rest.age <= 20)


restfact_targethill0_20y <- glmmTMB(target.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                     + (1 |region) + (1 |hydrology),
                     data = data_model_targethill0_20y,
                     family = nbinom2    
)
summary(restfact_targethill0_20y)
# --> restoration age is not significant


### b save final model ----

# save(restfact_targethill0_20y, data_model_targethill0_20y,
#      file = here("outputs", "models",
#                  "model_methods_target_hill0_20y.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D - FULL MODEL - NO RESTORATION AGE #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


rm(list = setdiff(ls(), c("data_all")))


## 1 Model formulation ---------------------------------------------------------

# include all restoration sites
data_model <- data_all %>%
  filter(!is.na(rest.meth))

# check NA
colSums(is.na(data_model)) 
# NA in no other variable than rest.age


B1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = poisson
)

check_overdispersion(B1)


# change to negative binomial family due to overdispersion
B1 <- glmmTMB(target.hill.0 ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = nbinom2    
)
summary(B1)


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

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus land.use.hist
ggplot(data = data_model, aes(x = land.use.hist, y = E_B1)) +
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

#' Plot residuals versus restoration age
ggplot(data = data_model, aes(x = rest.age.std, y = E_B1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### d check for adjustments of model -------------------------------------------

# no adjustments needed. model is fine


## 3 Model fitting -------------------------------------------------------------

# fit the model for reducing interactions,
# keep all main terms


### a backward selection (drop1) ----
drop1(B1)
drop_model <- update(B1, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)


### b drop model ----
B1_drop <- drop_model


## 4 Model validation drop model -----------------------------------------------

#' Get residuals and fitted values of model
data_model$E_B1_drop <- resid(B1_drop, type = "pearson")   
data_model$F_B1_drop <- fitted(B1_drop) 


### a Plot residuals vs fitted values ---------------------------------------------

ggplot(data = data_model, aes(x = F_B1_drop, y = E_B1_drop)) +
  geom_point(size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = FALSE) +  
  labs(x = "Fitted values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus land.use.hist
ggplot(data = data_model, aes(x = land.use.hist, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus region
ggplot(data = data_model, aes(x = region, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 

#' Plot the residuals versus hydrology
ggplot(data = data_model, aes(x = hydrology, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### c Plot residuals vs covariates not in the model ------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model, aes(x = rest.age.std, y = E_B1_drop)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 


### d Check collinearity part 2 ------------------------------------------------

# Remove VIF > 3 or > 10
# Zuur et al. 2010 Methods Ecol Evol DOI: 10.1111/j.2041-210X.2009.00001.x

check_collinearity(B1_drop)
# ok


### e check for adjustments of model -------------------------------------------

# no adjustments needed. model is fine


## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all")))


### a Summary ----

data_model_targethill0 <- data_all %>%
  filter(!is.na(rest.meth))

restfact_targethill0 <- glmmTMB(
  target.hill.0 ~ rest.meth + land.use.hist
  + (1 |region) + (1 |hydrology),
  data = data_model_targethill0,
  family = nbinom2 
)
summary(restfact_targethill0)
performance(restfact_targethill0)

### b Post-hoc test ----

# estimated marginal means (EMMs)
# Tukey-adjusted pariwise comparisons
emm.rest.meth <- emmeans(restfact_targethill0, ~ rest.meth)
summary(emm.rest.meth, infer = F, type = "response")
pairs(emm.rest.meth, adjust = "tukey") 


### c save final model ----

# save(restfact_targethill0, data_model_targethill0,
#      file = here("outputs", "models",
#                  "model_methods_target_hill0_full.Rdata"))


## end script