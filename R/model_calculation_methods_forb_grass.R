#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Response variable: Forb-grass-ratio
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
    id.site,
    site.type, hydrology, region, rest.meth, land.use.hist, rest.age,
    site.forb.legu.grass.ratio) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih"),
  )




## transform input data --------------------------------------------------------

# standardize numeric explanatory variables
data <- sites %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)), .names = "{col}.std"))



## set model data --------------------------------------------------------------


# join diversity data
data_all <- data %>%
  dplyr::select(id.site, site.forb.legu.grass.ratio, hydrology, region,
                rest.meth, rest.age.std, land.use.hist) %>% 
  rename(fg.ratio = site.forb.legu.grass.ratio)



rm(list = setdiff(ls(), c("data_all", "sites")))


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
  select(id.site, tot.hill.0) %>% 
  identify_outliers(tot.hill.0)
data_all %>% 
  select(id.site, rest.age.std) %>% 
  identify_outliers(rest.age.std)
sites %>% 
  filter(id.site %in% c("S_JAU", "S_NOZ", "M_JER", "S_MCH"))
sort(sites$site.cover.legumes)
# ok



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

ggplot(data_all, aes(x = rest.meth, y = fg.ratio)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = land.use.hist, y = fg.ratio)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent") 
ggplot(data_all, aes(x = rest.age, y = fg.ratio)) +
  geom_point() + geom_smooth(method = "glm") 
ggplot(data_all, aes(x = region, y = fg.ratio)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent") 
ggplot(data_all, aes(x = hydrology, y = fg.ratio)) +
  geom_jitter(color = "grey") + geom_boxplot(fill = "transparent")


## f distribution --------------------------------------------------------------

library(lattice)
library(fitdistrplus)
library(logspline)

histogram(data_all$fg.ratio)
histogram(log(data_all$fg.ratio))

x <- data_all$fg.ratio
descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.exp <- fitdist(x, "exp")
fit.gamma <- fitdist(x, "gamma")
fit.lnorm <- fitdist(x, "lnorm")
print(fit.gamma)

df_summary <- data_all %>%
  group_by(site.type) %>%
  summarize(mean_ratio = mean(fg.ratio), 
            var_ratio = var(fg.ratio))
ggplot(df_summary, aes(x = mean_ratio, y = var_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Variance vs. Mean of Forb-Grass Ratio")

plot(fit.norm)
plot(fit.exp)
plot(fit.gamma)
plot(fit.lnorm)

fit.norm$aic
fit.exp$aic
fit.gamma$aic
fit.lnorm$aic

detach(package:fitdistrplus)



## g Interactions --------------------------------------------------------------

# need to be checked: 2025-03-12

# check for possible interactions between covariates with coplot
coplot(fg.ratio ~ rest.age | region * hydrology,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

coplot(fg.ratio ~ rest.age | rest.meth * land.use.hist,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# same but different
# Grafische Daten Exploration ---
ggplot(data_all, aes(rest.age, fg.ratio, color = rest.meth)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~land.use.hist)
ggplot(data_all, aes(rest.age, fg.ratio, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~region)


## rest.meth vs. X
interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$land.use.hist,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(land.use.hist, rest.meth), y = fg.ratio)) + geom_boxplot()
# interaction


ggplot(data_all, aes(rest.age, fg.ratio, color = rest.meth)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# no clear interaction


interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$region,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(region, rest.meth), y = fg.ratio)) + geom_boxplot()
# interaction

interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$hydrology,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(rest.meth, hydrology), y = fg.ratio)) + geom_boxplot()
# interaction


## land.use.hist vs. X

ggplot(data_all, aes(rest.age, fg.ratio, color = land.use.hist)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# no clear interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$region,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(region, land.use.hist), y = fg.ratio)) + geom_boxplot()
# interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$hydrology,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(hydrology, land.use.hist), y = fg.ratio)) + geom_boxplot()
# no clear interaction



## rest.age vs. X

ggplot(data_all, aes(rest.age, fg.ratio, color = region)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction between rest.age and region --> consider random slope

ggplot(data_all, aes(rest.age, fg.ratio, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction between rest.age and hydrology --> consider random slope

## region vs. hydrology

interaction.plot(x.factor = data_all$hydrology, trace.factor = data_all$region,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(region, hydrology), y = fg.ratio)) + geom_boxplot()
# interaction



# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

## rest.meth vs. X

int_model <- glm(fg.ratio ~ rest.meth * land.use.hist, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction between rest.meth and land.use.hist

int_model <- glm(fg.ratio ~ rest.meth * rest.age, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction

int_model <- glm(fg.ratio ~ rest.meth * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction

int_model <- glm(fg.ratio ~ rest.meth * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction


## land.use.hist vs. X

int_model <- glm(fg.ratio ~ land.use.hist * rest.age, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction 

int_model <- glm(fg.ratio ~ land.use.hist * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

int_model <- glm(fg.ratio ~ land.use.hist * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction


## rest.age vs. X

int_model <- glm(fg.ratio ~ rest.age * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction between rest.age and region

int_model <- glm(fg.ratio ~ rest.age * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction between rest.age and hydrology


## region vs. hydrology

int_model <- glm(fg.ratio ~ region * hydrology, data = data_all, family = Gamma(link = "log"))
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



## i conclusions--------------------------------------------------------------

#' extreme outliers with fgratio > 5
#' missing values in rest.age and land.use.hist
#' detected collinearity between restoration age and restoration method --> consider removing rest.age
#' interaction between:
#' - hydrology and rest.meth, rest.age (-> consider random slope)
#' interaction unclear between:
#' - region and land.use.hist, rest.age (-> consider random slope)
#' - region and hydrology
#' no interaction between:
#' - rest.meth and land.use.hist, rest.age, region
#' - land.use.hist and rest.age, hydrology




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))


# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 
# only 108 sites left
# NA in rest.age


R1 <- glmmTMB(fg.ratio ~ 1 + (1|region), data = data_model, family = Gamma(link="log"))
R2 <- glmmTMB(fg.ratio ~ 1 + (1|hydrology), data = data_model, family = Gamma(link="log"))
R3 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|hydrology), data = data_model, family = Gamma(link="log"))
R4 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R5 <- glmmTMB(fg.ratio ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R6 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R7 <- glmmTMB(fg.ratio ~ 1 + (rest.age.std|region) + (1|hydrology), data = data_model, family = Gamma(link="log"))
# -> singularity
R8 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (rest.age.std|hydrology), data = data_model, family = Gamma(link="log"))
R9 <- glmmTMB(fg.ratio ~ 1 + (rest.age.std|region) + (rest.age.std|hydrology), data = data_model, family = Gamma(link="log"))
# -> singularity
Rnull <- glm(fg.ratio ~ 1, data = data_model, family = Gamma(link="log")) # right family??


AIC(R1, R2, R3, R4, R5, R6, R7, R8, R9, Rnull) %>% 
  arrange(AIC)
# --> model 3 is the same good as model 2
# --> go on to build model with model 3 as in other response variables

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(fg.ratio ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later


# # how strong is hydrology?
# glmm_hydr <- glmer.nb(fg.ratio ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(fg.ratio ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(fg.ratio ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(fg.ratio ~ rest.meth * hydrology + (1|region), data = data_model)
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


#' To account for colinearity between restoration age and restoration method we
#' would have to exclude one of the factors from analysis. First we remove colinerarity
#' and check for influence of restoration age in a reduced dataset with sites 
#' younger than 20 years.


## subset data to account for colinearity in age and restoration method
# exclude all sites > 20 years
data_model_20y <- data_all %>% 
  filter(rest.age <= 20) 
# 91 sites

# check NA
colSums(is.na(data_model_20y)) 
# no NA


## check colinearity
ggplot(data_model_20y, aes(x = rest.meth, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
# no colinearity

data_model_20y %>% 
  anova_test(rest.age ~ rest.meth)
# not significant, no colinearity



# beyound optimal model
#  mu_ij = Intercept + rest.meth_ij * land.use.hist_ij + rest.age_ij * land.use.hist_ij + region_i + hydrology_i

# consider interactions between explanatory variables

# B1_y crossed random intercept model with region and hydrology as random factors
yB1 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = Gamma(link="log")
)
check_overdispersion(yB1)
#' No overdispersion detected.

summary(yB1)
# restoration age is not significant



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
#' --> variance in residuals higher with medium fitted values?



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
#' slopes are not different for each hydrology


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okayish

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
yB1_sl1 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = Gamma(link="log")
)
# --> warning: singular convergence
yB1_sl2 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = Gamma(link="log")
)
# --> warning: singular convergence
yB1_sl3 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (0 + rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = Gamma(link="log")
)
yB1_sl4 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (0 + rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = Gamma(link="log")
)

AIC(yB1, yB1_sl1, yB1_sl2, yB1_sl3, yB1_sl4) %>% arrange(AIC)
# random slope is not better (less than delta 2) --> stay with model 


## 3 Model fitting -------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1")))


# fit the model for reducing interactions,
# keep all main terms

### a backward selection (drop1) ----
# with drop1
drop1(yB1)
drop_model <- update(yB1, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)
drop_model <- update(drop_model, . ~ . - land.use.hist:rest.age.std)
drop1(drop_model)
# all interactions are gone. stop here

drop1(drop_model, test = "Chisq")


### b all subsets regression (MuMIn) ----
# with MuMIn package
# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# # (for some reason using glmmTMB ends up with not only removing interactions)
# yB1a <- glmer.nb(fg.ratio ~ rest.meth * land.use.hist + rest.age.std * land.use.hist 
#                  + (1|region) + (1|hydrology), 
#                    data = data_model_20y)
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
# summary(drop_model)
# MuMIn and drop1 end up with the same model -> super

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
#' higher variance in medium values?



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
#' slopes are not really different for each region


#' Plot the residuals versus rest.meth
ggplot(data = data_model_20y, aes(x = rest.meth, y = E_yB1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' cus is not good

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
yB1_drop_sl1 <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist + rest.age.std
                        + (rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = Gamma(link="log")
)
# Model convergence problem
yB1_drop_sl2 <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = Gamma(link="log")
)
# Model convergence problem
yB1_drop_sl3 <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist + rest.age.std
                        + (0 + rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = Gamma(link="log")
)
yB1_drop_sl4 <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (0 + rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = Gamma(link="log")
)

AIC(yB1_drop, yB1_drop_sl1, yB1_drop_sl2, yB1_drop_sl3, yB1_drop_sl4) %>% arrange(AIC)
# sl3 model is slightly better, but < delta 2 --> stay with model 



## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1", "yB1_drop")))


### a Summary ----

# the model we use is:
# fg.ratio ~ rest.meth + land.use.hist + rest.age.std + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation","model_plants_restfact_hill0_yB1_final.Rdata"))


data_model_fgratio_20y <- data_all %>%
  filter(rest.age <= 20) 


restfact_fgratio_20y <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist + rest.age.std
                                 + (1 |region) + (1 |hydrology),
                                 data = data_model_fgratio_20y,
                                 family = Gamma(link="log")
)
summary(restfact_fgratio_20y)
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)            -0.35668    0.45137  -0.790    0.429
# rest.methmga            0.17020    0.50763   0.335    0.737
# rest.methres            0.49859    0.38705   1.288    0.198
# rest.methdih            0.45233    0.38606   1.172    0.241
# land.use.histgrassland -0.27896    0.23297  -1.197    0.231
# rest.age.std            0.08399    0.16619   0.505    0.613

# --> restoration age is not significant

performance(restfact_fgratio_20y)
#' ICC = 0.359
#' The covariates explain 7 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 40 % of the variation. (R2 cond.)


### b report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_fgratio_20y, effects = "fixed") %>%
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
#       "outputs", "tables", "vegetation", "model_summary_restfact_fgratio_20y.csv"
#     )
#   )


### c save final model ----

save(restfact_fgratio_20y, data_model_fgratio_20y,
     file = here("outputs", "models",
                 "model_methods_forb_grass_20y.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E - FULL MODEL - NO RESTORATION AGE #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

## 1 Model formulation ---------------------------------------------------------


# include all restoration sites
data_model <- data_all %>%
  filter(!is.na(rest.meth))
# 121 sites

# check NA
colSums(is.na(data_model)) 
# NA in no other variable than rest.age


#  mu_ij = Intercept + rest.meth_ij * land.use.hist_ij + region_i + hydrology_i

# beyound optimal model
## consider interactions between explanatory variables


# B1 crossed random intercept model with region and hydrology as random factors
B1 <- glmmTMB(fg.ratio ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = Gamma(link="log")
)

#' Check overispersion
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
#' maybe higher variance in residuals with medium values?



### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okayish

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
#' stronger variance in lower values?


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
# B1a <- glmer.nb(fg.ratio ~ rest.meth * land.use.hist 
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
# # fg.ratio ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# 
# options(na.action = "na.omit") # set back to default


### c comparison fitting procedure ----

# AIC(full_model, drop_model, top_model)
# summary(top_model)
# summary(drop_model)
# MuMIn and drop1 end up with the same model -> super

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
#' maybe higher variance in residuals with smaller values?



### b Plot residuals vs covariates in the model --------------------------------

#' Plot the residuals versus rest.meth
ggplot(data = data_model, aes(x = rest.meth, y = E_B1_drop)) +
  geom_boxplot() + 
  labs(x = "factor", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks okayish

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
#' looks okeyish


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
# fg.ratio ~ rest.meth + land.use.hist + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio.Rdata"))

data_model_fgratio <- data_all %>%
  filter(!is.na(rest.meth)) 


restfact_fgratio <- glmmTMB(fg.ratio ~ rest.meth + land.use.hist
                             + (1 |region) + (1 |hydrology),
                             data = data_model_fgratio,
                            family = Gamma(link="log")
)
summary(restfact_fgratio)$sigma
# sigma < 1 --> gamma distribution is okay, and not exponential (sigma ~= 1)
summary(restfact_fgratio)
# Estimate Std. Error z value Pr(>|z|)  
# (Intercept)            -0.29562    0.31436  -0.940   0.3470  
# rest.methmga           -0.16148    0.30262  -0.534   0.5936  
# rest.methres            0.47800    0.23112   2.068   0.0386 *
#   rest.methdih            0.36559    0.22882   1.598   0.1101  
# land.use.histgrassland -0.02449    0.17565  -0.139   0.8891  

# --> rest.meth is significant, land.use.hist not

check_overdispersion(restfact_fgratio)
# no overdispersion detected

performance(restfact_fgratio)
#' ICC = 0.294
#' The covariates explain 9 % of the variation (R2 marg.)
#' The covariates and random effects explain 36 % of the variation. (R2 cond.)


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.rest.meth <- emmeans(restfact_fgratio, ~ rest.meth)
emm.rest.meth
summary(emm.rest.meth, infer = F, type = "response")
# rest.meth response    SE  df
# cus          0.735 0.239 Inf
# mga          0.625 0.201 Inf
# res          1.185 0.362 Inf
# dih          1.059 0.304 Inf
# 
# Results are averaged over the levels of: land.use.hist 
pairs(emm.rest.meth, adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast  estimate    SE  df z.ratio p.value
# cus - mga    0.161 0.303 Inf   0.534  0.9509
# cus - res   -0.478 0.231 Inf  -2.068  0.1636
# cus - dih   -0.366 0.229 Inf  -1.598  0.3799
# mga - res   -0.639 0.268 Inf  -2.390  0.0790
# mga - dih   -0.527 0.244 Inf  -2.158  0.1351
# res - dih    0.112 0.200 Inf   0.561  0.9435
# 
# Results are averaged over the levels of: land.use.hist 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

pairs(regrid(emm.rest.meth), adjust = "tukey") # regrid() for back-transformation from log-scale

## --> different p-values, because of calculating them after back-transformation



### c report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_fgratio, effects = "fixed") %>%
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
#       "outputs", "tables", "vegetation", "model_summary_restfact_fgratio.csv"
#     )
#   )

### summary statistics ###
data_model_tot %>%
  group_by(rest.meth) %>%
  get_summary_stats(fg.ratio, type = "full")

### d save final model ----


save(restfact_fgratio, data_model_fgratio,
     file = here("outputs", "models", "vegetation",
                 "model_methods_forb_grass_full.Rdata"))






## end script











glmm_rest_age <- glmer.nb(fg.ratio ~ rest.age + (1|region) + (1|hydrology), data = data_model)
summary(glmm_rest_age)

glmm_rest_5 <- glmer.nb(fg.ratio ~
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
ggplot(data_model, aes(rest.age, fg.ratio, color = region)) +
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





glmm_rest_age <- glmer.nb(fg.ratio ~ rest.age + (1|region) + (1|hydrology), data = data_model)
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
ggplot(data_model, aes(rest.age, fg.ratio, color = region)) +
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
  fg.ratio ~ rest.meth + land.use.hist + 
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
# second best model (by AICc) is fg.ratio ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# AIC 821.2433
# which is very very close to top model!
# Rule of thumb often used : models with AIC.delta <=2 are equally supported by the data

# get second best model
B1_top2 <- get.models(model_dredge, subset = 2)[[1]]


# Average model
av_model <- model.avg(model_dredge, subset = delta < 2)
summary(av_model)

options(na.action = "na.omit") # set back to default


B1a <- glm(fg.ratio ~ rest.meth, data = data_model, family = "poisson")

AICc(B1a, B1_top)







## performance package

check_model(top_model)



## lmerTest

## backward stepwise regression ####
bw_model <- step(full_model, direction = "backward")






## subset models ---------------------------------------------------------------

glmm_dih_age <- glmer.nb(fg.ratio ~ rest.age + (1|region) + (1|hydrology), data = data_dih)
glmm_mga_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_mga)
glmm_res_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_res)
glmm_cus_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_cus)
summary(glmm_dih_age)
summary(glmm_mga_age)
summary(glmm_res_age)
summary(glmm_cus_age)

plot(data_all$fg.ratio ~ data_all$rest.age)

glmm_dry_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_dry)
glmm_fresh_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_fresh)
glmm_moist_age <- glmer.nb(fg.ratio ~ rest.age + (1|region), data = data_moist)
summary(glmm_dry_age)
summary(glmm_fresh_age)
summary(glmm_moist_age)


