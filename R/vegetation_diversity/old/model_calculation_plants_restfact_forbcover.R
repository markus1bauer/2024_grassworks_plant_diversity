#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Response variable: relative cover of forbs (forb cover)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATIION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(lme4)
library(lmerTest)
library(DHARMa)
library(MuMIn) # automated model selection
library(performance) # visual check of model assumptions
library(ggbeeswarm)
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(rstatix)
library(glmmTMB)
library(mgcv)
library(broom.mixed) # Tidy up the model summary
library(scales) # label_number()



### Start ###
rm(list = ls())
setwd(here())

### Functions ###

# glm.nb (negative-binominal distribution) due to overdispersion
# library(MASS) # for glm.nb -> detach afterwards because it masks "select"

# test for overdispersion
# With function from Ben Bolker
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion 
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
# overdisp_fun(R_3a)
# significant result shows overdispersion
# --> use negative binomial-distribution



## load data -------------------------------------------------------------------

### site environment data ####

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20241114.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site,
    site.type, hydrology, region, rest.meth, land.use.hist, rest.age,
    longitude, latitude
  ) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih"),
  )


forb_index <- read_csv(
  here("data", "processed", "data_processed_plants_forbindex_20241129.csv"),
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

data <- data %>% 
  left_join(forb_index, by = "id.site")


data_all <- data %>%
  dplyr::select(
    id.site,
    forb.rel.cover.site,
    hydrology,
    region,
    rest.meth,
    rest.age,
    rest.age.std,
    land.use.hist,
    longitude,
    latitude
  ) %>% 
  rename(forb.cover = forb.rel.cover.site)

# remove zeros for beta distribution (doesn't allow zeros)
data_all <- data_all %>% 
  mutate(forb.cover = forb.cover + 0.001)



rm(list = setdiff(ls(), c("data_all", "overdisp_fun", "sites")))


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
# (plot2 <- ggplot(data_all, aes(x = forb.cover)) +
#   geom_histogram(binwidth = 1))
# plot3 <- ggplot(data_all, aes(x = forb.cover)) +
#   geom_density()
# plot4 <- ggplot(data_all, aes(x = log(forb.cover))) +
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

sort(data_all$forb.cover)

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
  select(id.site, forb.cover) %>% 
  identify_outliers(forb.cover)
# no outliers
data_all %>% 
  select(id.site, forb.cover, rest.meth, land.use.hist) %>% 
  group_by(rest.meth) %>% 
  identify_outliers(forb.cover)
data_all %>% 
  select(id.site, forb.cover, rest.meth, land.use.hist) %>% 
  group_by(land.use.hist) %>% 
  identify_outliers(forb.cover)

data_all %>% 
  select(id.site, rest.age) %>% 
  identify_outliers(rest.age)
# no extreme outliers
data_all %>% 
  select(id.site, rest.age, rest.meth) %>% 
  group_by(rest.meth) %>% 
  identify_outliers(rest.age)
data_all %>% 
  select(id.site, rest.age, land.use.hist) %>% 
  group_by(land.use.hist) %>% 
  identify_outliers(rest.age)




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


ggplot(data_all, aes(x = rest.meth, y = forb.cover)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Restoration method")
ggplot(data_all, aes(x = land.use.hist, y = forb.cover)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Previous landuse")
ggplot(data_all, aes(x = rest.age, y = forb.cover)) +
  geom_point() + geom_smooth(method = "glm") +
  labs(title = "Age of Restoration")
ggplot(data_all, aes(x = region, y = forb.cover)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = forb.cover)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# relationship between rest.age and response linear?
# difference between regions --> use as random factor
# difference between hydrology --> test as random factor


## f distribution --------------------------------------------------------------
library(lattice)
histogram(data_all$forb.cover)
# normal?

x <- data_all$forb.cover

library(fitdistrplus)
library(logspline)

descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.beta <- fitdist(x, "beta")


plot(fit.beta)
plot(fit.norm)
# beta and normal looks okeyish

fit.norm$aic
fit.gamma$aic
fit.beta$aic
# beta has lowes AIC

library(gamlss)
library(gamlss.dist)
library(gamlss.add)

fit <- fitDist(x, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)

summary(fit)


## f Interactions --------------------------------------------------------------

# check for possible interactions between covariates with coplot
coplot(forb.cover ~ rest.age | region * hydrology,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

coplot(forb.cover ~ rest.age | rest.meth * land.use.hist,
       data = data_all,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# same but different
# Grafische Daten Exploration ---
ggplot(data_all, aes(rest.age, forb.cover, color = rest.meth)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~land.use.hist)
ggplot(data_all, aes(rest.age, forb.cover, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~region)


# TODO check interactions

## rest.meth vs. X
interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$land.use.hist,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(land.use.hist, rest.meth), y = forb.cover)) + geom_boxplot()
# interaction


ggplot(data_all, aes(rest.age, forb.cover, color = rest.meth)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# no interaction


interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$region,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(region, rest.meth), y = forb.cover)) + geom_boxplot()
# interaction

interaction.plot(x.factor = data_all$rest.meth, trace.factor = data_all$hydrology,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(rest.meth, hydrology), y = forb.cover)) + geom_boxplot()
# interaction


## land.use.hist vs. X

ggplot(data_all, aes(rest.age, forb.cover, color = land.use.hist)) +
  geom_point() + theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# no interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$region,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(region, land.use.hist), y = forb.cover)) + geom_boxplot()
# interaction

interaction.plot(x.factor = data_all$land.use.hist, trace.factor = data_all$hydrology,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(hydrology, land.use.hist), y = forb.cover)) + geom_boxplot()
# interaction



## rest.age vs. X

ggplot(data_all, aes(rest.age, forb.cover, color = region)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction between rest.age and region --> consider random slope

ggplot(data_all, aes(rest.age, forb.cover, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F)
# interaction between rest.age and hydrology --> consider random slope

## region vs. hydrology

interaction.plot(x.factor = data_all$hydrology, trace.factor = data_all$region,
                 response = data_all$forb.cover)
ggplot(data_all, aes(x = interaction(region, hydrology), y = forb.cover)) + geom_boxplot()
# interaction



# test interactions between covariates
# interaction term significant --> interaction

## rest.meth vs. X
int_model <- glm(forb.cover ~ rest.meth * land.use.hist, data = data_all, family = beta_family(link = "logit"))
anova(int_model)
# --> no interaction between rest.meth and land.use.hist

int_model <- glm(forb.cover ~ rest.meth * rest.age, data = data_all, family = beta_family(link = "logit"))
anova(int_model)
# --> no interaction

int_model <- glm(forb.cover ~ rest.meth * region, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# --> no interaction

int_model <- glm(forb.cover ~ rest.meth * hydrology, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction


## land.use.hist vs. X

int_model <- glm(forb.cover ~ land.use.hist * rest.age, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# no interaction 

int_model <- glm(forb.cover ~ land.use.hist * region, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

int_model <- glm(forb.cover ~ land.use.hist * hydrology, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# interaction


## rest.age vs. X

int_model <- glm(forb.cover ~ rest.age * region, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# interaction between rest.age and region

int_model <- glm(forb.cover ~ rest.age * hydrology, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# --> interaction between rest.age and hydrology


## region vs. hydrology

int_model <- glm(forb.cover ~ region * hydrology, data = data_all, family = beta_family(link = "logit"))
check_overdispersion(int_model)
anova(int_model)
# no interaction


## g Spatial dependency --------------------------------------------------------------

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



## h conclusions--------------------------------------------------------------

#' missing values in rest.age and land.use.hist
#' detected collinearity between restoration age and restoration method --> consider removing rest.age
#' unclearity if there is a linear relationship between rest.age and forb.cover
#' interaction between:
#' - hydrology and rest.meth, land.use.hist, rest.age (-> consider random slope)
#' - region and rest.age
#' interaction unclear between:
#' - rest.meth and land.use.hist, region
#' - region and hydrology
#' no interaction between:
#' - rest.meth and rest.age
#' - land.use.hist and rest.age, region




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all", "overdisp_fun")))


# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 
# only 108 sites left
# NA in rest.age


R1 <- glmmTMB(forb.cover ~ 1 + (1|region), data = data_model, family = beta_family(link = "logit"))
R2 <- glmmTMB(forb.cover ~ 1 + (1|hydrology), data = data_model, family = beta_family(link = "logit"))
R3 <- glmmTMB(forb.cover ~ 1 + (1|region) + (1|hydrology), data = data_model, family = beta_family(link = "logit"))
R4 <- glmmTMB(forb.cover ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = beta_family(link = "logit"))
R5 <- glmmTMB(forb.cover ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = beta_family(link = "logit"))
R6 <- glmmTMB(forb.cover ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = beta_family(link = "logit"))
R7 <- glmmTMB(forb.cover ~ 1 + (rest.age.std|region) + (1|hydrology), data = data_model, family = beta_family(link = "logit"))
# -> Model convergence problem
R8 <- glmmTMB(forb.cover ~ 1 + (1|region) + (rest.age.std|hydrology), data = data_model, family = beta_family(link = "logit"))
# -> Model convergence problem
R9 <- glmmTMB(forb.cover ~ 1 + (rest.age.std|region) + (rest.age.std|hydrology), data = data_model, family = beta_family(link = "logit"))
# -> Model convergence problem
Rnull <- glm(forb.cover ~ 1, data = data_model, family = beta_family(link = "logit")) # right family??


AIC(R1, R2, R3, R4, R5, R6,
    # R7, R8, R9, 
    Rnull
) %>% 
  arrange(AIC)
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(forb.cover ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later


# # how strong is hydrology?
# glmm_hydr <- glmer.nb(forb.cover ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(forb.cover ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(forb.cover ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(forb.cover ~ rest.meth * hydrology + (1|region), data = data_model)
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
yB1 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
               + (1|region) + (1|hydrology), 
               data = data_model_20y,
               family = beta_family(link = "logit")
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
#' looks great



### b Plot residuals vs covariates in the model --------------------------------

#' Plot residuals versus restoration age
ggplot(data = data_model_20y, aes(x = rest.age.std, y = E_yB1)) +
  geom_point( size = 0.8, alpha = 0.5) + 
  geom_smooth(method = "gam", se = TRUE) +  
  labs(x = "values", y = "Residuals") + 
  geom_hline(yintercept = 0, lty = 2) 
#' looks fine

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
yB1_sl1 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = beta_family(link = "logit")
)
# --> Model convergence problem
yB1_sl2 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = beta_family(link = "logit")
)
# --> Model convergence problem
yB1_sl3 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (0 + rest.age.std|region) + (1|hydrology), 
                   data = data_model_20y,
                   family = beta_family(link = "logit")
)
yB1_sl4 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist
                   + (1|region) + (0 + rest.age.std|hydrology), 
                   data = data_model_20y,
                   family = beta_family(link = "logit")
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
# yB1a <- glmer.nb(forb.cover ~ rest.meth * land.use.hist + rest.age.std * land.use.hist 
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
#' looks fine



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
yB1_drop_sl1 <- glmmTMB(forb.cover ~ rest.meth + land.use.hist + rest.age.std
                        + (rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = beta_family(link = "logit")
)
# Model convergence problem
yB1_drop_sl2 <- glmmTMB(forb.cover ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = beta_family(link = "logit")
)
# Model convergence problem
yB1_drop_sl3 <- glmmTMB(forb.cover ~ rest.meth + land.use.hist + rest.age.std
                        + (0 + rest.age.std|region) + (1|hydrology), 
                        data = data_model_20y,
                        family = beta_family(link = "logit")
)
yB1_drop_sl4 <- glmmTMB(forb.cover ~ rest.meth + land.use.hist + rest.age.std
                        + (1|region) + (0 + rest.age.std|hydrology), 
                        data = data_model_20y,
                        family = beta_family(link = "logit")
)

AIC(yB1_drop, yB1_drop_sl1, yB1_drop_sl2, yB1_drop_sl3, yB1_drop_sl4) %>% 
  arrange(AIC)
# random slope is not better --> stay with model 



## 5 Final model ---------------------------------------------------------------

rm(list = setdiff(ls(), c("data_all", "data_model_20y", "yB1", "yB1_drop")))


### a Summary ----

# the model we use is:
# forb.cover ~ rest.meth + land.use.hist + rest.age.std + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation","model_plants_restfact_hill0_yB1_final.Rdata"))


data_model_forbcover_20y <- data_all %>%
  filter(rest.age <= 20)

restfact_forbcover_20y <- glmmTMB(forb.cover ~ rest.meth + land.use.hist + rest.age.std
                                  + (1 |region) + (1 |hydrology),
                                  data = data_model_forbcover_20y,
                                  family = beta_family(link = "logit")
)
summary(restfact_forbcover_20y)
# Estimate Std. Error z value Pr(>|z|)   
# (Intercept)            -1.28423    0.42936  -2.991  0.00278 **
#   rest.methmga           -0.01146    0.40366  -0.028  0.97735   
# rest.methres            0.20777    0.32061   0.648  0.51696   
# rest.methdih            0.57109    0.33210   1.720  0.08550 . 
# land.use.histgrassland  0.13465    0.19640   0.686  0.49297   
# rest.age.std           -0.19997    0.15017  -1.332  0.18298  

# --> restoration age is not significant

performance(restfact_forbcover_20y)
#' ICC = 0.567
#' The covariates explain 11 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 61 % of the variation. (R2 cond.)


### b report final model ----

# Tidy up the model summary
model_summary <- broom.mixed::tidy(restfact_forbcover_20y, effects = "fixed") %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
  mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
  mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
                             .default = p.value))

# Display the summary as a table
library(kableExtra)
model_summary %>%
  kbl(caption = "GLMM Model Results") %>%
  kable_styling(full_width = F)

library(gt)
model_summary %>%
  gt() %>%
  tab_header(title = "GLMM Model Results")

model_summary %>% 
  write_csv(
    here(
      "outputs", "tables", "vegetation", "model_summary_restfact_forbcover_20y.csv"
    )
  )


### c save final model ----

save(restfact_forbcover_20y, data_model_forbcover_20y,
     file = here("outputs", "models", "vegetation",
                 "model_plants_restfact_forbcover_20y.Rdata"))



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
B1 <- glmmTMB(forb.cover ~ rest.meth * land.use.hist
              + (1|region) + (1|hydrology), 
              data = data_model,
              family = beta_family(link = "logit")
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
#' looks okay



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
#' looks okay


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
drop1(drop_model, test = "Chisq")
# all interactions are gone. stop here


### b all subsets regression (MuMIn) ----
# with MuMIn package

# options(na.action = "na.fail") # Required for dredge to run
# 
# # use glmer.nb() for dredge
# B1a <- glmer.nb(forb.cover ~ rest.meth * land.use.hist 
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
# # forb.cover ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
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
#' looks okay



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
#' looks ok


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
# forb.cover ~ rest.meth + land.use.hist + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restfact_hill0_B1_final.Rdata"))

data_model_forbcover <- data_all %>%
  filter(!is.na(rest.meth))

restfact_forbcover <- glmmTMB(forb.cover ~ rest.meth + land.use.hist
                              + (1 |region) + (1 |hydrology),
                              data = data_model_forbcover,
                              family = beta_family(link = "logit")
)
summary(restfact_forbcover)
# Estimate Std. Error z value Pr(>|z|)   
# (Intercept)             -1.0024     0.3121  -3.212  0.00132 **
#   rest.methmga            -0.2606     0.2846  -0.916  0.35988   
# rest.methres             0.1085     0.1993   0.544  0.58627   
# rest.methdih             0.3405     0.2218   1.535  0.12469   
# land.use.histgrassland   0.2926     0.1635   1.790  0.07343 . 

# --> rest.meth is NOT significant, land.use.hist on the edge



### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.rest.meth <- emmeans(restfact_forbcover, ~ rest.meth)
emm.rest.meth
summary(emm.rest.meth, infer = F, type = "response")
pairs(emm.rest.meth, adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast  estimate    SE  df z.ratio p.value
# cus - mga    0.261 0.285 Inf   0.916  0.7965
# cus - res   -0.108 0.199 Inf  -0.544  0.9481
# cus - dih   -0.341 0.222 Inf  -1.535  0.4162
# mga - res   -0.369 0.244 Inf  -1.513  0.4295
# mga - dih   -0.601 0.237 Inf  -2.541  0.0538
# res - dih   -0.232 0.183 Inf  -1.271  0.5817
# 
# Results are averaged over the levels of: land.use.hist 
# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 4 estimates 

pairs(regrid(emm.rest.meth), adjust = "tukey") # regrid() for back-transformation from log-scale


## --> different p-values, because of calculating them after back-transformation



### c report final model ----

# Tidy up the model summary
model_summary <- broom.mixed::tidy(restfact_forbcover, effects = "fixed") %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  mutate(across(c(estimate, std.error, statistic), round, 3)) %>%   # Rounding values
  mutate(p.value = label_number(accuracy = 0.0001)(p.value)) %>% 
  mutate(p.value = case_when(p.value < 0.001 ~ "< 0.001",
                             .default = p.value))

# Display the summary as a table
library(kableExtra)
model_summary %>%
  kbl(caption = "GLMM Model Results") %>%
  kable_styling(full_width = F)

library(gt)
model_summary %>%
  gt() %>%
  tab_header(title = "GLMM Model Results")

model_summary %>% 
  write_csv(
    here(
      "outputs", "tables", "vegetation", "model_summary_restfact_forbcover.csv"
    )
  )

### summary statistics ###
data_model_forbcover %>%
  group_by(rest.meth) %>%
  get_summary_stats(forb.cover, type = "full")

### d save final model ----


save(restfact_forbcover, data_model_forbcover,
     file = here("outputs", "models", "vegetation",
                 "model_plants_restfact_forbcover.Rdata"))






## end script











glmm_rest_age <- glmer.nb(forb.cover ~ rest.age + (1|region) + (1|hydrology), data = data_model)
summary(glmm_rest_age)

glmm_rest_5 <- glmer.nb(forb.cover ~
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
ggplot(data_model, aes(rest.age, forb.cover, color = region)) +
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





glmm_rest_age <- glmer.nb(forb.cover ~ rest.age + (1|region) + (1|hydrology), data = data_model)
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
ggplot(data_model, aes(rest.age, forb.cover, color = region)) +
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
  forb.cover ~ rest.meth + land.use.hist + 
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
# second best model (by AICc) is forb.cover ~ land.use.hist + rest.meth + (1 | region) + (1 | hydrology)
# AIC 821.2433
# which is very very close to top model!
# Rule of thumb often used : models with AIC.delta <=2 are equally supported by the data

# get second best model
B1_top2 <- get.models(model_dredge, subset = 2)[[1]]


# Average model
av_model <- model.avg(model_dredge, subset = delta < 2)
summary(av_model)

options(na.action = "na.omit") # set back to default


B1a <- glm(forb.cover ~ rest.meth, data = data_model, family = "poisson")

AICc(B1a, B1_top)







## performance package

check_model(top_model)



## lmerTest

## backward stepwise regression ####
bw_model <- step(full_model, direction = "backward")






## subset models ---------------------------------------------------------------

glmm_dih_age <- glmer.nb(forb.cover ~ rest.age + (1|region) + (1|hydrology), data = data_dih)
glmm_mga_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_mga)
glmm_res_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_res)
glmm_cus_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_cus)
summary(glmm_dih_age)
summary(glmm_mga_age)
summary(glmm_res_age)
summary(glmm_cus_age)

plot(data_all$forb.cover ~ data_all$rest.age)

glmm_dry_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_dry)
glmm_fresh_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_fresh)
glmm_moist_age <- glmer.nb(forb.cover ~ rest.age + (1|region), data = data_moist)
summary(glmm_dry_age)
summary(glmm_fresh_age)
summary(glmm_moist_age)


