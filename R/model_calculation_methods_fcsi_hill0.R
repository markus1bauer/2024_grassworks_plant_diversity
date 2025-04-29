#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Response variable: FCSi (Index of Favourable Conservation Status)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATION ###############################################################
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
    id.site, site.type, hydrology, region, rest.meth, land.use.hist, rest.age,
    longitude, latitude
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

# standardise explanatory variable (only numerical variables)
# so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”)
# --> z-standardization
# It ensures that the estimated coefficients are all on the same scale, 
# making it easier to compare effect sizes.

# Standardizing the numeric explanatory variables
data <- sites %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.)), .names = "{col}.std"))

# # Verify scaling
# summary(data)




## set model data --------------------------------------------------------------

# join diversity data
data_all <- data %>%
  left_join(diversity, by = "id.site")


data_all <- data_all %>%
  dplyr::select(
    id.site,
    hydrology,
    region,
    rest.meth,
    rest.age,
    rest.age.std,
    land.use.hist,
    longitude,
    latitude,
    fcsi.hill.0
  )









rm(list = setdiff(ls(), c("data_all")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Protocol of data exploration (Steps 1-8)
# used from Zuur et al. (2010) Methods Ecol Evol 
#[DOI: 10.1111/2041-210X.12577](https://doi.org/10.1111/2041-210X.12577)


## a Missing values ------------------------------------------------------------
colSums(is.na(data_all)) 
# 66 NAs in rest.meth, rest.age and land.use.hist due to reference sites
# 13 other NAs in rest.age
vis_dat(data_all)
gg_miss_var(data_all)  




## b Outliers, zero-inflation, transformations? (Step 1, 3, 4) -----------------

# data_all %>%
#   count(location_construction_year)
# plot1 <- ggplot(data_all, aes(x = exposition, y = y)) +
#   geom_quasirandom()
# (plot2 <- ggplot(data_all, aes(x = fcsi.hill.0)) +
#   geom_histogram(binwidth = 1))
# plot3 <- ggplot(data_all, aes(x = fcsi.hill.0)) +
#   geom_density()
# plot4 <- ggplot(data_all, aes(x = log(fcsi.hill.0))) +
#   geom_density()
# (plot1 + plot2) / (plot3 + plot4)





# Outliers: check with Cleveland dotplot

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
# outlier are points far right or left in plot
# doesn't look like there are outliers

sort(data_all$fcsi.hill.0)

## wie umgehen mit ouliers?? Messfehler: unrealistische WErte
# 
# dotchart(sites_data$fcsi.hill.2,
#          ylab = "Order of the data")
# 
# sites_data %>% 
#   filter(fcsi.hill.2 > 5) %>% 
#   select(id.site)

# another test for outliers
library(rstatix)
data_all %>% 
  dplyr::select(id.site, fcsi.hill.0) %>% 
  identify_outliers(fcsi.hill.0)
data_all %>% 
  dplyr::select(id.site, rest.age) %>% 
  identify_outliers(rest.age)
# no extreme outliers


## c inspect categorical covariates -----------------------------------------

table(data_all$rest.meth)
#' Unbalanced...but enough observations per level.

table(data_all$land.use.hist)
#' Unbalanced...but enough observations per level.

table(data_all$hydrology)
#' Unbalanced...but enough observations per level.

table(data_all$region)
#' Balanced.

#' Was each previous land use measured in every restoration method?
table(data_all$rest.meth, data_all$land.use.hist)
histogram( ~ rest.meth | land.use.hist, data_all)
#' Unbalanced, do we have enough observations per combination for an 
#' interaction term?

#' Was each hydrology measured in every restoration method?
table(data_all$rest.meth, data_all$hydrology)
histogram( ~ rest.meth | hydrology, data_all)
#' Unbalanced, do we have enough observations per combination for an 
#' interaction term?

#' Was each region measured in every restoration method?
table(data_all$rest.meth, data_all$region)
histogram( ~ rest.meth | region, data_all)
#' Unbalanced, no mga in south

#' Was each hydrology measured in every previous land use?
table(data_all$land.use.hist, data_all$hydrology)
histogram( ~ land.use.hist | hydrology, data_all)
#' Unbalanced, but enough observations per combination for an interaction term

#' Was each region measured in every previous land use?
table(data_all$land.use.hist, data_all$region)
histogram( ~ land.use.hist | region, data_all)
#' Unbalanced, but fine

#' Was each region measured in every hydrology?
table(data_all$hydrology, data_all$region)
histogram( ~ hydrology | region, data_all)
#' Unbalanced, but fine



## d Check collinearity part 1 ----------------------------------------

# between continuous covariates

# only one numerical variable in model data --> no need to check

# Exclude r > 0.7
#  Dormann et al. 2013 Ecography [DOI: 10.1111/j.1600-0587.2012.07348.x](https://doi.org/10.1111/j.1600-0587.2012.07348.x)

# data_all %>%
#   select(where(is.numeric) %>%
#   GGally::ggpairs(
#     lower = list(continuous = "smooth_loess")
#   ) +
#   theme(strip.text = element_text(size = 7))
#   )
#
# # exclude variable
# data_all <- data_all %>%
#   select(-biotope_area)

# between continuous variables and factors

ggplot(data_all, aes(x = rest.meth, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
data_all %>% 
  anova_test(rest.age ~ rest.meth)
#' boxplots are not next to each other, anova is significant:
#' --> restoration method is collinear with age
ggplot(data_all, aes(x = land.use.hist, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = hydrology, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
ggplot(data_all, aes(x = region, y = rest.age)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
#' region is unbalanced, but fine I guess.
#' The rest is fine.



# source("HighstatLibV14.R")
# MyVar <- c("rest.meth", "rest.age", "land.use.hist", "hydrology", "region")
# corvif(data_all[,MyVar])


## e Relationships --------------------------------------------------------------

#' Plot response variable versus each covariate.


ggplot(data_all, aes(x = rest.meth, y = fcsi.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Restoration method")
ggplot(data_all, aes(x = land.use.hist, y = fcsi.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Previous landuse")
ggplot(data_all, aes(x = rest.age, y = fcsi.hill.0)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Age of Restoration")
ggplot(data_all, aes(x = region, y = fcsi.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = fcsi.hill.0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# relationship between rest.age and response linear?
# difference between regions --> use as random factor
# difference between hydrology --> test as random factor

## f distribution --------------------------------------------------------------

library(lattice)
histogram(data_all$fcsi.hill.0)
# normal distribution?

x <- data_all$fcsi.hill.0

library(fitdistrplus)
library(logspline)

descdist(x.std, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.gamma <- fitdist(x, "gamma")


plot(fit.norm)
plot(fit.gamma)

# normal looks okeyish

fit.norm$aic
fit.gamma$aic
# normal distribution has lowest AIC

detach(package:fitdistrplus)



## g Interactions --------------------------------------------------------------

# TODO: not checked yet --> check! (2025-03-12)

# check for possible interactions between covariates with coplot
coplot(fcsi.hill.0 ~ rest.age | region * hydrology,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

coplot(fcsi.hill.0 ~ rest.age | rest.meth * land.use.hist,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# same but different
# Grafische Daten Exploration ---
ggplot(data_all, aes(rest.age, fcsi.hill.0, color = rest.meth)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~land.use.hist)
ggplot(data_all, aes(rest.age, fcsi.hill.0, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~region)


## rest.meth vs. X
interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$land.use.hist,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(land.use.hist, rest.meth), y = fcsi.hill.0)) + geom_boxplot()
# interaction


ggplot(data_all, aes(rest.age, fcsi.hill.0, color = rest.meth)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), se = F)
# interaction


interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$region,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(region, rest.meth), y = fcsi.hill.0)) + geom_boxplot()
# interaction

interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$hydrology,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(hydrology, rest.meth), y = fcsi.hill.0)) + geom_boxplot()
# interaction


## land.use.hist vs. X

ggplot(data_all, aes(rest.age, fcsi.hill.0, color = land.use.hist)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction between rest.age and land.use.hist

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$region,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(region, land.use.hist), y = fcsi.hill.0)) + geom_boxplot()
# no interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$hydrology,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(hydrology, land.use.hist), y = fcsi.hill.0)) + geom_boxplot()
# interaction



## rest.age vs. X

ggplot(data_all, aes(rest.age, fcsi.hill.0, color = region)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), se = F)
# interaction between rest.age and region --> consider random slope

ggplot(data_all, aes(rest.age, fcsi.hill.0, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), se = F)
# interaction between rest.age and hydrology --> consider random slope

## region vs. hydrology

interaction.plot(x.factor = data_all$hydrology, trace.factor = data_all$region,
                 response = data_all$fcsi.hill.0)
ggplot(data_all, aes(x = interaction(region, hydrology), y = fcsi.hill.0)) + geom_boxplot()
# no interaction



# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

## rest.meth vs. X

int_model <- glm(fcsi.hill.0 ~ rest.meth * land.use.hist, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction between rest.meth and land.use.hist

int_model <- glm(fcsi.hill.0 ~ rest.meth * rest.age, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction

int_model <- glm(fcsi.hill.0 ~ rest.meth * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction

int_model <- glm(fcsi.hill.0 ~ rest.meth * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction


## land.use.hist vs. X

int_model <- glm(fcsi.hill.0 ~ land.use.hist * rest.age, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# interaction 

int_model <- glm(fcsi.hill.0 ~ land.use.hist * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

int_model <- glm(fcsi.hill.0 ~ land.use.hist * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction


## rest.age vs. X

int_model <- glm(fcsi.hill.0 ~ rest.age * region, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction between rest.age and region

int_model <- glm(fcsi.hill.0 ~ rest.age * hydrology, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction between rest.age and hydrology


## region vs. hydrology

int_model <- glm(fcsi.hill.0 ~ region * hydrology, data = data_all, family = Gamma(link = "log"))
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

#' missing values in rest.age and land.use.hist
#' detected collinearity between restoration age and restoration method --> consider removing rest.age
#' unclearity if there is a linear relationship between rest.age and fcsi.hill.0
#' interaction between:
#' - rest.meth and rest.age (-> consider random slope)
#' - hydrology and rest.meth, rest.age (-> consider random slope), region
#' interaction unclear between:
#' - region and land.use.hist, rest.age (-> consider random slope)
#' - land.use.hist and rest.age (-> consider random slope)
#' no interaction between:
#' - rest.meth and land.use.hist, rest.age, region
#' - land.use.hist and hydrology




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))


# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 
# only 108 sites left
# NA in rest.age


R1 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|region), data = data_model, family = gaussian)
R2 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|hydrology), data = data_model, family = gaussian)
R3 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = gaussian)
R4 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = gaussian)
R5 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = gaussian)
R6 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = gaussian)
R7 <- glmmTMB(fcsi.hill.0 ~ 1 + (rest.age.std|region) + (1|hydrology), data = data_model, family = gaussian)
R8 <- glmmTMB(fcsi.hill.0 ~ 1 + (1|region) + (rest.age.std|hydrology), data = data_model, family = gaussian)
# -> singularity
R9 <- glmmTMB(fcsi.hill.0 ~ 1 + (rest.age.std|region) + (rest.age.std|hydrology), data = data_model, family = gaussian)
# -> singularity
Rnull <- glm(fcsi.hill.0 ~ 1, data = data_model, family = gaussian) # right family??


AIC(R1, R2, R3, R4, R5, R6, R7, R8, R9, Rnull) %>% 
  arrange(AIC)
# model 7 is similar to model 3 (< delta 2), model 3 is easier
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(fcsi.hill.0 ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later


# # how strong is hydrology?
# glmm_hydr <- glmer.nb(fcsi.hill.0 ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(fcsi.hill.0 ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(fcsi.hill.0 ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(fcsi.hill.0 ~ rest.meth * hydrology + (1|region), data = data_model)
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
yB1 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = gaussian)

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
#' --> outlier with < -1

data_model_20y %>% 
  filter(E_yB1 < -1)
# N_DAM - maybe exclude from analysis?

## --> 2025-03-12 not excluded, calculated model with N_DAM site as well ##
#' 
#' ## 3 Model adjustment -----------------------------------------------
#' 
#' # exclude N_DAM
#' data_model_20y <- data_model_20y %>% 
#'   filter(id.site != "N_DAM")
#' 
#' # B1_y crossed random intercept model with region and hydrology as random factors
#' yB1 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
#'                + (1|region) + (1|hydrology), 
#'                data = data_model_20y,
#'                family = gaussian
#' )
#' check_overdispersion(yB1)
#' #' No overdispersion detected.
#' 
#' summary(yB1)
#' # restoration age is not significant
#' 
#' 
#' ## 4 Model validation full model -----------------------------------------------
#' 
#' ### a Plot residuals vs fitted values ---------------------------------------------
#' 
#' ggplot(data = data_model_20y, aes(x = F_yB1, y = E_yB1)) +
#'   geom_point(size = 0.8, alpha = 0.5) + 
#'   geom_smooth(method = "gam", se = FALSE) +  
#'   labs(x = "Fitted values", y = "Residuals") + 
#'   geom_hline(yintercept = 0, lty = 2) 
#' #' Is there a pattern in here?
#' #' negative fitted values?
#' #' looks good


### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' --> looks okayish

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
#' slopes are similar for each region

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
yB1_sl1 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = gaussian
)
# --> warning: singular convergence
yB1_sl2 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = gaussian
)
# --> warning: singular convergence
yB1_sl3 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (0 + rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = gaussian
)
yB1_sl4 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (0 + rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = gaussian
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
drop_model <- update(drop_model, . ~ . - rest.meth:land.use.hist)
drop1(drop_model)
# all interactions are gone. stop here

drop1(drop_model, test = "Chisq")


### b all subsets regression (MuMIn) ----
# with MuMIn package
# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# # (for some reason using glmmTMB ends up with not only removing interactions)
# yB1a <- glmer.nb(fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std * land.use.hist 
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
#' looks good



### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1_drop)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks ok

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
#' slopes are similar for each region

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
# moderate collinearity
# keep them in for now


# use glmer.nb() model
# car::vif(top_model)
# ok


### f check for adjustments of model -------------------------------------------

# check for need of random slope
yB1_drop_sl1 <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = gaussian
)
yB1_drop_sl2 <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = gaussian
)
yB1_drop_sl3 <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (0 + rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = gaussian
)
yB1_drop_sl4 <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (0 + rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = gaussian
)

AIC(yB1_drop, yB1_drop_sl1, yB1_drop_sl2, yB1_drop_sl3, yB1_drop_sl4) %>% arrange(AIC)
# random slope is not better --> stay with model 



## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1", "yB1_drop")))


### a Summary ----

# the model we use is:
# fcsi.hill.0 ~ rest.meth * land.use.hist + rest.age.std + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation","model_plants_restfact_hill0_yB1_final.Rdata"))


data_model_fcsihill0_20y <- data_all %>%
  filter(rest.age <= 20) #%>% 
  # filter(id.site != "N_DAM")

restfact_fcsihill0_20y <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist + rest.age.std
                                 + (1 |region) + (1 |hydrology),
                                 data = data_model_fcsihill0_20y,
                                 family = gaussian
)
summary(restfact_fcsihill0_20y)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             3.17921    0.17643  18.020  < 2e-16 ***
#   rest.methmga            0.04866    0.15359   0.317  0.75140    
# rest.methres            0.21997    0.12445   1.768  0.07713 .  
# rest.methdih            0.36867    0.13110   2.812  0.00492 ** 
#   land.use.histgrassland  0.13967    0.07912   1.765  0.07752 .  
# rest.age.std            0.01485    0.06047   0.246  0.80603   

# --> restoration age is not significant

performance(restfact_fcsihill0_20y)
#' ICC = 0.412
#' The covariates explain 13 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 49 % of the variation. (R2 cond.)


### b report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_fcsihill0_20y, effects = "fixed") %>%
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
#       "outputs", "tables", "vegetation", "model_summary_restfact_fcsihill0_20y.csv"
#     )
#   )


### c save final model ----

save(restfact_fcsihill0_20y, data_model_fcsihill0_20y,
     file = here("outputs", "models",
                 "model_methods_fcsi_hill0_20y.Rdata"))



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
B1 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = gaussian
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
#' maybe higher variance in residuals with lower values?

data_model %>% 
  filter(E_B1 < -1)
# N_DAM - maybe exclude from analysis?

## --> 2025-03-12 not excluded, calculated model with N_DAM site as well ##

#' ## 3 Model adjustment -----------------------------------------------
#' 
#' # exclude N_DAM
#' data_model <- data_model %>% 
#'   filter(id.site != "N_DAM")
#' 
#' # B1 crossed random intercept model with region and hydrology as random factors
#' B1 <- glmmTMB(fcsi.hill.0 ~ rest.meth * land.use.hist
#'               + (1|region) + (1|hydrology), 
#'               data = data_model,
#'               family = gaussian
#' )
#' 
#' #' Check overispersion
#' check_overdispersion(B1)
#' #' No overdispersion detected.
#' 
#' 
#' 
#' ## 2 Model validation full model -----------------------------------------------
#' 
#' 
#' #' As part of the model validation, we need to:
#' #'  -Plot residuals versus fitted values.
#' #'  -Plot residuals versus each covariate in the model.
#' #'  -Plot residuals versus each covariate NOT in the model.
#' #'  -Plot residuals versus time (if relevant).
#' #'  -Plot residuals versus spatial coordinates (if relevant).
#' 
#' 
#' #' Get residuals and fitted values of model
#' data_model$E_B1 <- resid(B1, type = "pearson")   #' Observed eveness minus fitted values.
#' data_model$F_B1 <- fitted(B1)  #' Fitted values contain the random effects.
#' 
#' 
#' ### a Plot residuals vs fitted values ---------------------------------------------
#' 
#' 
#' ggplot(data = data_model, aes(x = F_B1, y = E_B1)) +
#'   geom_point(size = 0.8, alpha = 0.5) + 
#'   geom_smooth(method = "gam", se = FALSE) +  
#'   labs(x = "Fitted values", y = "Residuals") + 
#'   geom_hline(yintercept = 0, lty = 2) 
#' #' Is there a pattern in here?
#' #' negative fitted values?
#' #' looks ok

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

drop1(drop_model, test = "Chisq")


### b all subsets regression (MuMIn) ----
# with MuMIn package

# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# B1a <- glmer.nb(fcsi.hill.0 ~ rest.meth * land.use.hist 
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
# # fcsi.hill.0 ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
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
#' looks okayish



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
#' maybe higher variance in smaller values


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
# fcsi.hill.0 ~ rest.meth + land.use.hist + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restfact_hill0_B1_final.Rdata"))

data_model_fcsihill0 <- data_all %>%
  filter(!is.na(rest.meth)) #%>% 
  # filter(id.site != "N_DAM")

restfact_fcsihill0 <- glmmTMB(fcsi.hill.0 ~ rest.meth + land.use.hist
                             + (1 |region) + (1 |hydrology),
                             data = data_model_fcsihill0,
                             family = gaussian
)
summary(restfact_fcsihill0)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             3.19054    0.15846  20.135  < 2e-16 ***
#   rest.methmga            0.08641    0.10438   0.828  0.40773    
# rest.methres            0.21419    0.07626   2.809  0.00498 ** 
#   rest.methdih            0.32579    0.08307   3.922 8.78e-05 ***
#   land.use.histgrassland  0.12963    0.06226   2.082  0.03732 * 

# --> rest.meth is significant, land.use.hist is also significant (!)

performance(restfact_fcsihill0)
#' ICC = 0.460
#' The covariates explain 12 % of the variation (R2 marg.)
#' The covariates and random effects explain 52 % of the variation. (R2 cond.)


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.rest.meth <- emmeans(restfact_fcsihill0, ~ rest.meth)
emm.rest.meth
summary(emm.rest.meth, infer = F, type = "response")
# rest.meth emmean    SE  df
# cus         3.26 0.161 113
# mga         3.34 0.160 113
# res         3.47 0.155 113
# dih         3.58 0.153 113
# 
# Results are averaged over the levels of: land.use.hist 
pairs(emm.rest.meth, adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast  estimate     SE  df t.ratio p.value
# cus - mga  -0.0864 0.1044 113  -0.828  0.8411
# cus - res  -0.2142 0.0763 113  -2.809  0.0295
# cus - dih  -0.3258 0.0831 113  -3.922  0.0009
# mga - res  -0.1278 0.0936 113  -1.366  0.5235
# mga - dih  -0.2394 0.0896 113  -2.671  0.0425
# res - dih  -0.1116 0.0718 113  -1.555  0.4088
# 
# Results are averaged over the levels of: land.use.hist 
# P value adjustment: tukey method for comparing a family of 4 estimates  

pairs(regrid(emm.rest.meth), adjust = "tukey") # regrid() for back-transformation from log-scale

## --> different p-values, because of calculating them after back-transformation



### c report final model ----

# # Tidy up the model summary
# model_summary <- broom.mixed::tidy(restfact_fcsihill0, effects = "fixed") %>%
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
#       "outputs", "tables", "vegetation", "model_summary_restfact_fcsihill0.csv"
#     )
#   )

### summary statistics ###
data_model_fcsihill0 %>%
  group_by(rest.meth) %>%
  get_summary_stats(fcsi.hill.0, type = "full")

### d save final model ----


save(restfact_fcsihill0, data_model_fcsihill0,
     file = here("outputs", "models",
                 "model_methods_fcsi_hill0_full.Rdata"))






## end script











glmm_rest_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_model)
summary(glmm_rest_age)

glmm_rest_5 <- glmer.nb(fcsi.hill.0 ~
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
ggplot(data_model, aes(rest.age, fcsi.hill.0, color = region)) +
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





glmm_rest_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_model)
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
ggplot(data_model, aes(rest.age, fcsi.hill.0, color = region)) +
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
  fcsi.hill.0 ~ rest.meth + land.use.hist + 
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
# second best model (by AICc) is fcsi.hill.0 ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# AIC 821.2433
# which is very very close to top model!
# Rule of thumb often used : models with AIC.delta <=2 are equally supported by the data

# get second best model
B1_top2 <- get.models(model_dredge, subset = 2)[[1]]


# Average model
av_model <- model.avg(model_dredge, subset = delta < 2)
summary(av_model)

options(na.action = "na.omit") # set back to default


B1a <- glm(fcsi.hill.0 ~ rest.meth, data = data_model, family = Gamma(link = "log"))

AICc(B1a, B1_top)







## performance package

check_model(top_model)



## lmerTest

## backward stepwise regression ####
bw_model <- step(full_model, direction = "backward")






## subset models ---------------------------------------------------------------

glmm_dih_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region) + (1|hydrology), data = data_dih)
glmm_mga_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_mga)
glmm_res_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_res)
glmm_cus_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_cus)
summary(glmm_dih_age)
summary(glmm_mga_age)
summary(glmm_res_age)
summary(glmm_cus_age)

plot(data_all$fcsi.hill.0 ~ data_all$rest.age)

glmm_dry_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_dry)
glmm_fresh_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_fresh)
glmm_moist_age <- glmer.nb(fcsi.hill.0 ~ rest.age + (1|region), data = data_moist)
summary(glmm_dry_age)
summary(glmm_fresh_age)
summary(glmm_moist_age)


