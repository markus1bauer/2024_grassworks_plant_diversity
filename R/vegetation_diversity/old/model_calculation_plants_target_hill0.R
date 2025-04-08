#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Response variable: Characteristic Species Richness
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
#library(vegan)
#library(nlme)
library(lme4)
library(lmerTest)
library(DHARMa)
library(MuMIn) # automated model selection
library(performance) # visual check of model assumptions


### Start ###
rm(list = ls())

## load data -------------------------------------------------------------------

### site environment data ####

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20240813.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, rest.meth, region, rest.age, 
    land.use.hist, hydrology, hydrology.cont, site.cwm.abu.oek.f, lui,
    mngm.type, obs.year,
    c.n, c.perc, n.perc, ph.value, tic.perc, toc.perc, ends_with(c("B10", "B30")),
    clay.perc, silt.perc, sand.perc,
  ) %>%
  distinct() %>% 
  # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
  filter(!is.na(id.site))

sites$obs.year <- as.factor(sites$obs.year)

# todo: resolve infinite numbers in c.n.B30 and c.n in sites N_RET and N_HOR
# quick and dirty:
sites[sites$id.site == "N_RET", "c.n"] <- NA
sites[sites$id.site == "N_HOR", "c.n"] <- NA
sites[sites$id.site == "N_RET", "c.n.B30"] <- NA
sites[sites$id.site == "N_HOR", "c.n.B30"] <- NA

# choose one restoration method for "dih&seed"
# --> dih because it is most likely that dih has the greatest influence
# --> Line said for M_SKF: ReS
sites %>% 
  filter(rest.meth == "dih&seed")
sites[sites$id.site == "S_FHZ", "rest.meth"] <- "dih"
sites[sites$id.site == "S_GIG", "rest.meth"] <- "dih"
sites[sites$id.site == "M_SKF", "rest.meth"] <- "res"
sites[sites$id.site == "S_WTZ", "rest.meth"] <- "dih"





### diversity data ####

diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20240814.csv"),
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
# 
# colSums(is.na(data))
# # M_WDG is missing
# # lui 11 rows are missing (not M_WDG)
# # start.rest 88 are missing (66 are reference sites)
# # some other: resolve later



## set model data --------------------------------------------------------------


# join diversity data
data_all <- data %>%
  left_join(diversity, by = "id.site")


# 
# # only restored sites
# data_rest <- data_all %>% 
#   filter(site.type == "restored") %>% # one restored site has NA
#   # remove dih&seed sites (check if ok!)
#   filter(rest.meth != "dih&seed") # 4 sites with rest.meth = "dih&seed"



data_model <- data_all %>%
  dplyr::select(
    target.hill.0,
    hydrology,
    region,
    rest.meth,
    rest.age.std,
    land.use.hist,
    # site.cwm.abu.oek.f.std,
    lui.std,
    # mngm.type,
    c.n.std,
    c.perc.std,
    # n.perc.std,
    ph.value.std,
    # tic.perc.std,
    # toc.perc.std,
    # clay.perc.std,
    # silt.perc.std,
    sand.perc.std,
    site.type
  )



# # eliminate all rows with any missing values
data_model <- na.omit(data_model) 
# %>% dplyr::filter(rest.meth != "dih&seed")
# only 98 sites left: how to keep more rows?
# solve c.n NA issue

data_all %>% 
  count(rest.meth, region)



plot(data_all$region, data_all$ph.value)
plot(data_all$hydrology, data_all$ph.value)


rm(list = setdiff(ls(), c("data_model", "data_rest")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Data Exploration ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# see Zuur et al. 2010

## 1 - Outliers ----------------------------------------------------------------

# check with Cleveland dotplot

library(lattice)
Z <- data_all %>% 
  select(ends_with("std"), starts_with(c("tot", "target", "fcsi")))


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
# outlier in response variables: fcsi.hill.1, fcsi.hill.2
# outlier in explanatory variables: 
#   toc.perc.std, n.perc.std, c.perc.std.
#   toc.perc.B10.std, n.perc.B10.std, c.n.B10.std, c.perc.B10.std
#   toc.perc.B30.std, n.perc.B30.std, c.perc.B30.std
#   lui.std
# c.n.B30.std and c.n.std missing !
# maybe in rest.age ?

sort(data_all$n.perc)

## wie umgehen mit ouliers?? Messfehler: unrealistische WErte

dotchart(sites_data$fcsi.hill.2,
         ylab = "Order of the data")

sites_data %>% 
  filter(fcsi.hill.2 > 5) %>% 
  select(id.site)

# another test for outliers
library(rstatix)
sites %>% 
  select(id.site, lui) %>% 
  identify_outliers(lui)

# LUI: N_TUT, N_EIC --> uncertainty in survey data in M_EIC -> LUI maybe incorrect
# toc: N_PFE
# n: N_PFE, S_MAS
# c.n.B10.std: N_EIC, S_TRN
# c: N_PFE
# fcsi.hill.1: N_JAS, N_NLG
# fcsi.hill.2: M_KOT, N_NLG


## 2 - Collinearity -------------------------------------------------------------

library(Hmisc)

soil_data <- data_all %>% 
  select(ends_with(c("B10.std", "B30.std")))
corr_result <- rcorr(as.matrix(soil_data), type = c("spearman"))
corr_result
library(PerformanceAnalytics)
chart.Correlation(soil_data, histogram=T, method = c ("spearman"), pch=19)
# soil data from B10 and B30 are highly correlated --> use only mean values!

# all numeric variables (incl. mean soil variables)
data <- data_all %>%
  select(ends_with("std"), -ends_with(c("B10.std", "B30.std"))) 


corr_result <- rcorr(as.matrix(data), type = c("spearman"))
corr_result
# looking in spearmann correlation for selecting variables (correlation < +|- 0.6)


chart.Correlation(data, histogram=T, method = c ("spearman"), pch=19)

# correlations >= 0.6
# sand with clay and silt --> take sand
# cwm f value with hydrology.cont --> skip hydrology.cont
# c/n moderately with n (-0.55), not with c --> keep
# c with toc and moderately with tic (0.53) (but tic not with toc) and n --> take c
# ph with tic --> take ph




## 3 - Interactions -------------------------------------------------------------

# check for possible interactions between covariates with coplot
coplot(target.hill.0 ~ rest.age | region * hydrology,
       data = sites_data,
       ylab = "Total species richness",
       xlab = "",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Modelling #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all", "data_rest")))




## different models -------------------------------------------------------------

# glm.nb (negative-binominal distribution) due to overdispersion
library(MASS) # for glm.nb -> detach afterwards because it masks "select"

# todo: test this again!

### random factor ####

glmm_random_1 <- glmer.nb(target.hill.0 ~ 1 + (1|region), data = data_model)
glmm_random_2 <- glmer.nb(target.hill.0 ~ 1 + (1|hydrology), data = data_model)
glmm_random_3 <- glmer.nb(target.hill.0 ~ 1 + (1|region) + (1|hydrology), data = data_model)
glmm_random_4 <- glmer.nb(target.hill.0 ~ 1 + (1|region) + (1|region:hydrology), data = data_model)
glmm_random_5 <- glmer.nb(target.hill.0 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model)
# boundary (singular) fit: see help('isSingular')
glmm_random_6 <- glmer.nb(target.hill.0 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model)
glm_null <- glm(target.hill.0 ~ 1, data = data_model)


AIC(glmm_random_1, glmm_random_2, glmm_random_3, glmm_random_4, glmm_random_5,
    glmm_random_6,
    glm_null) %>% 
  arrange(AIC)
# df      AIC
# glmm_random_3  4 762.8285
# glmm_random_4  4 763.2123
# glmm_random_6  5 763.4610
# glmm_random_5  4 765.8226
# glmm_random_1  3 772.4974
# glmm_random_2  3 796.3223
# glm_null       2 807.1700
# --> model 3 (region and hydrology as random factors, no interaction) has lowest AIC (1492.334)
# --> go on to build model with this random structure



### only restoration variables ####
## single variables
glmm_rest_meth <- glmer.nb(target.hill.0 ~ rest.meth + (1|region) + (1|hydrology), data = data_model)
glmm_rest_age <- glmer.nb(target.hill.0 ~ rest.age.std + (1|region) + (1|hydrology), data = data_model)
glmm_land_use_hist <- glmer.nb(target.hill.0 ~ land.use.hist + (1|region) + (1|hydrology), data = data_model)

AIC(glmm_rest_meth, glmm_rest_age, glmm_land_use_hist) %>% 
  arrange(AIC)
# --> rest.meth has lowest AIC
summary(glmm_rest_meth)


## combination of restoration variables
glmm_rest <- glmer.nb(target.hill.0 ~ rest.meth + rest.age.std + land.use.hist
                      + (1|region) + (1|hydrology), data = data_model)
glmm_rest_1 <- glmer.nb(target.hill.0 ~ rest.meth + rest.age.std
                        + (1|region) + (1|hydrology), data = data_model)
glmm_rest_2 <- glmer.nb(target.hill.0 ~ rest.meth + land.use.hist
                        + (1|region) + (1|hydrology), data = data_model)
glmm_rest_3 <- glmer.nb(target.hill.0 ~ rest.age.std + land.use.hist
                        + (1|region) + (1|hydrology), data = data_model)

AIC(glmm_rest_meth,
    glmm_rest, glmm_rest_1, glmm_rest_2, glmm_rest_3) %>% 
  arrange(AIC)
# rest.meth + land.use.hist is best model


## interactions
glmm_rest_4 <- glmer.nb(target.hill.0 ~
                          rest.meth
                        : rest.age.std
                        + (1|region) + (1|hydrology)
                        , data = data_model
)
glmm_rest_5 <- glmer.nb(target.hill.0 ~
                          rest.meth
                        * rest.age.std
                        + (1|region) + (1|hydrology)
                        , data = data_model
)
glmm_rest_6 <- glmer.nb(target.hill.0 ~
                          rest.meth
                        : land.use.hist
                        + (1|region) + (1|hydrology)
                        , data = data_model
)
glmm_rest_7 <- glmer.nb(target.hill.0 ~
                          rest.meth
                        * land.use.hist
                        + (1|region) + (1|hydrology)
                        , data = data_model
)
AIC(glmm_rest_meth, glmm_rest_2,
    glmm_rest_4, glmm_rest_5, glmm_rest_6, glmm_rest_7) %>% 
  arrange(AIC)
# interaction between rest.meth and land.use.hist is best model


### only site conditions ####
## single variables
glmm_lui <- glmer.nb(target.hill.0 ~ lui.std + (1|region) + (1|hydrology), data = data_model)
glmm_cn <- glmer.nb(target.hill.0 ~ c.n.std + (1|region) + (1|hydrology), data = data_model)
glmm_c <- glmer.nb(target.hill.0 ~ c.perc.std + (1|region) + (1|hydrology), data = data_model)
glmm_ph <- glmer.nb(target.hill.0 ~ ph.value.std + (1|region) + (1|hydrology), data = data_model)
glmm_sand <- glmer.nb(target.hill.0 ~ sand.perc.std + (1|region) + (1|hydrology), data = data_model)

AIC(glmm_lui, glmm_cn, glmm_c, glmm_ph, glmm_sand) %>% 
  arrange(AIC)
# cn has lowest AIC, then sand, ph, lui, c


## combination of site variables
# all variables
glmm_site <- glmer.nb(target.hill.0 ~
                        lui.std
                      + c.n.std
                      + c.perc.std
                      + ph.value.std
                      + sand.perc.std
                      + (1|region) + (1|hydrology)
                      , data = data_model)


# all subsets regression with site variables (dredge)
options(na.action = "na.fail") # Required for dredge to run
model_dredge <- dredge(glmm_site, beta = "none", evaluate = T, trace = 2,
                       # fixed = "rest.meth",
                       # m.lim =c(0,5),
                       rank = AIC) # when do you use AICc?
top_model <- get.models(model_dredge, subset = 1)[[1]]
summary(top_model)
# top model: target.hill.0 ~ (1 | region) + (1 | hydrology)
# --> null-model
options(na.action = "na.omit") # set back to default



# combination of site variables with lowest single AIC
glmm_site_1 <- glmer.nb(target.hill.0 ~ c.n.std + ph.value.std
                        + (1|region) + (1|hydrology)
                        , data = data_model)
glmm_site_2 <- glmer.nb(target.hill.0 ~ c.n.std + sand.perc.std
                        + (1|region) + (1|hydrology)
                        , data = data_model)
glmm_site_3 <- glmer.nb(target.hill.0 ~ c.n.std + lui.std
                        + (1|region) + (1|hydrology)
                        , data = data_model)

AIC(glmm_cn, glmm_random_3,
    glmm_site_1, glmm_site_2, glmm_site_3) %>% 
  arrange(AIC)
# null-model is best

## comparison only rest and only site
AIC(glmm_cn, glmm_site, glmm_random_3,
    glmm_rest_meth, glmm_rest_2) %>% 
  arrange(AIC)
# only site variables is not better
# --> restoration (and explicitly method) is main factor effecting species richness


### combination of restoration and site conditions ####
glmm_rest_2_cn <- update(glmm_rest_2, .~. + c.n.std)
glmm_rest_2_lui <- update(glmm_rest_2, .~. + lui.std)
glmm_rest_2_c <- update(glmm_rest_2, .~. + c.perc.std)
glmm_rest_2_ph <- update(glmm_rest_2, .~. + ph.value.std)
glmm_rest_2_sand <- update(glmm_rest_2, .~. + sand.perc.std)

AIC(glmm_rest_2_lui, glmm_rest_2_cn, glmm_rest_2_c, glmm_rest_2_ph,
    glmm_rest_2_sand,
    glmm_rest_meth, glmm_rest_2) %>% 
  arrange(AIC)
# rest.meth + land.use.hist is best model
# 2nd best: rest.meth + ph
# rest.meth + lui is even better then rest.meth + cn


# glmm_rest_lui <- update(glmm_rest, .~. + lui.std)
# glmm_rest_cn <- update(glmm_rest, .~. + c.n.std)
# glmm_rest_c <- update(glmm_rest, .~. + c.perc.std)
# glmm_rest_ph <- update(glmm_rest, .~. + ph.value.std)
# glmm_rest_sand <- update(glmm_rest, .~. + sand.perc.std)
# 
# AIC(glmm_rest_lui, glmm_rest_cn, glmm_rest_c, glmm_rest_ph, glmm_rest_sand,
#     glmm_rest) %>% 
#   arrange(AIC)

# combine restored and ph, lui, cn
glmm_rest_2_site_1 <- update(glmm_rest_2, .~. + ph.value.std + lui.std)
glmm_rest_2_site_2 <- update(glmm_rest_2, .~. + ph.value.std + c.n.std)
glmm_rest_2_site_3 <- update(glmm_rest_2, .~. + lui.std + c.n.std)
AIC(glmm_rest_2_site_1, glmm_rest_2_site_2, glmm_rest_2_site_3,
    glmm_rest_2_ph, glmm_rest_2, glmm_rest_2_lui, glmm_rest_2_cn) %>% 
  arrange(AIC)
# rest.meth + land.use.hist is best model
# rest.meth + land.use.hist + ph is 2nd best
summary(glmm_rest_2)



sort(data_model$ph.value)


## Interaction between rest.meth and pH
glmm_rest_meth_ph_int_1 <-  glmer.nb(target.hill.0 ~ rest.meth*ph.value.std
                                     + (1|region) + (1|hydrology)
                                     , data = data_model)
glmm_rest_meth_ph_int_2 <-  glmer.nb(target.hill.0 ~ rest.meth:ph.value.std
                                     + (1|region) + (1|hydrology)
                                     , data = data_model)
summary(glmm_rest_meth_ph_int_1)
summary(glmm_rest_meth_ph_int_2)
AIC(glmm_rest_meth_ph_int_1, glmm_rest_meth_ph_int_2,
    glmm_rest_meth_ph) %>% 
  arrange(AIC)
# interaction doesn't make it better

## model evaluation ------------------------------------------------------------
res <- simulateResiduals(glmm_rest_meth_ph, plot = T)

op <- par(mfrow = c(2,2))
plot(glmm_rest_meth_ph)
par(op)

check_model(glmm_rest_meth_ph)




# -------------------------------------------------------------------------





glm_rest <- glm(
  target.hill.0 ~
    rest.meth
  + rest.age.std
  + land.use.hist
  , data = model_data
)
summary(glm_rest)
# --> overdispersion

glm_rest_1 <- glm.nb(
  target.hill.0 ~
    rest.meth
  + rest.age.std
  + land.use.hist
  , data = model_data,
)
summary(glm_rest_1)
# AIC: 745.13

# interaction
glm_rest_2 <- glm.nb(
  target.hill.0 ~
    rest.meth
  * rest.age.std
  + land.use.hist
  # + rest.meth*rest.age.std
  , data = model_data,
)
summary(glm_rest_2)
drop1(glm_rest_2)


glmm_rest_1 <- glmer.nb(
  target.hill.0 ~
    rest.meth
  + rest.age.std
  + land.use.hist
  + (1|region)
  , data = model_data, family = "poisson",
)



# only site conditions
# glm.nb due to overdispersion in glm
glm_site <- glm.nb(
  target.hill.0 ~
    hydrology
  # + site.cwm.abu.oek.f.std
  + lui.std
  # + mngm.type
  + c.n.std
  + c.perc.std
  # + n.perc.std
  + ph.value.std
  # + tic.perc.std
  # + toc.perc.std
  # + clay.perc.std
  # + silt.perc.std
  + sand.perc.std
  + region
  , data = model_data
)
# summary(glm_site)

glmm_site_1 <- glmer.nb(
  target.hill.0 ~
    hydrology
  # + site.cwm.abu.oek.f.std
  + lui.std
  # + mngm.type
  + c.n.std
  + c.perc.std
  # + n.perc.std
  + ph.value.std
  # + tic.perc.std
  # + toc.perc.std
  # + clay.perc.std
  # + silt.perc.std
  + sand.perc.std
  + (1|region)
  , data = model_data, family = "poisson"
)

glmm_site_2 <- glmer.nb(
  target.hill.0 ~
    # + site.cwm.abu.oek.f.std
    lui.std
  # + mngm.type
  + c.n.std
  + c.perc.std
  # + n.perc.std
  + ph.value.std
  # + tic.perc.std
  # + toc.perc.std
  # + clay.perc.std
  # + silt.perc.std
  + sand.perc.std
  + (1|region)
  + (1|hydrology)
  , data = model_data, family = "poisson"
)

# null model
glm_null <- glm.nb(
  target.hill.0 ~
    1
  , data = model_data
)

AIC(glm_rest_meth, glm_rest_age, glm_land_use_hist, glm_region,
    glmm_rest_meth_1, glmm_rest_meth_2,
    glm_rest_1, glm_rest_2, glm_site,
    glmm_rest_1, glmm_rest_2, glmm_site_1, glmm_site_2,
    glm_null)
# Warnmeldung:
#   In AIC.default(m_rest, m_site) :
#   Modelle sind nicht alle mit der gleichen Datensatzgröße angepasst worden




## fit full model ####
full_model <- glm(
  target.hill.0 ~
    rest.meth
  + rest.age.std
  + land.use.hist
  + hydrology
  # + site.cwm.abu.oek.f.std
  + lui.std
  # + mngm.type
  # + c.n.std
  + c.perc.std
  # + n.perc.std
  + ph.value.std
  # + tic.perc.std
  # + toc.perc.std
  # + clay.perc.std
  # + silt.perc.std
  + sand.perc.std
  + region
  , data = data_model
  # , na.action = na.omit,
)
summary(full_model)
# we see clearly overdispersion in summary
# Residual deviance / df = 7890.2/76 = 103.8186
# change to GLM with negative-binominal distribution

library(MASS) # for glm.nb -> detach afterwards because it masks "select"
full_model <- glmer.nb(
  target.hill.0 ~
    rest.meth
  + rest.age.std
  + land.use.hist
  # + site.cwm.abu.oek.f.std
  + lui.std
  # + mngm.type
  + c.n.std
  + c.perc.std
  # + n.perc.std
  + ph.value.std
  # + tic.perc.std
  # + toc.perc.std
  # + clay.perc.std
  # + silt.perc.std
  + sand.perc.std
  + (1|region) + (1|region:hydrology)
  , data = data_model
  # , na.action = na.omit,
)

# library(car)
vif(full_model)

detach("package:MASS", unload = TRUE)

summary(full_model)
# Mit der overdispersion das ist viel besser
# 95.057 / 76 = 1.25 --> unter 1.5 gilt als OK

# diagnostic plots
check_model(full_model)
# VIF of rest.meth is moderate high (> 5 but < 10)
# but I cannot skip this factor as it is my main explanatory variable! ok??
# observed residual variance does not follow predicted residual variance: ?? what does this mean??
op <- par(mfrow = c(2,2))
plot(full_model)
par(op)
# looks ok to me




## all subsets regression ####
# with MuMIn package
options(na.action = "na.fail") # Required for dredge to run

model_dredge <- dredge(full_model, beta = "none", evaluate = T, trace = 2,
                       # fixed = "rest.meth",
                       # m.lim =c(0,5),
                       rank = AIC) # when do you use AICc?
top_model <- get.models(model_dredge, subset = 1)[[1]]
summary(top_model)
summary(glmm_rest_site_1)
# with CWM F:
# glm.nb(formula = target.hill.0 ~ ph.value.std + region + rest.meth + 
#          site.cwm.abu.oek.f.std + 1)
# AIC: 684.37

# with hydrology:
# glm.nb(formula = target.hill.0 ~ hydrology + land.use.hist + ph.value.std + 
#          region + rest.age.std + rest.meth + sand.perc.std + rest.age.std:rest.meth + 
#          1, data = model_data, init.theta = 42.94222405, link = log)
# AIC: 664.06
print(model_dredge)
op <- par(mfrow = c(2,2))
plot(glmm_rest_site_1)
par(op)

check_model(glmm_rest_site_1)
anova(top_model)

subset(model_dredge, delta < 2)
par(mar = c(3,5,6,4))
plot(model_dredge, labAsExpr = TRUE)

top_most_models <- get.models(model_dredge, subset = delta < 2)
top_most_models

av_model <- model.avg(model_dredge, subset = delta < 2)
summary(av_model)

options(na.action = "na.omit") # set back to default




#______________________________

drop1(full_model, test = "Chi")
summary(full_model)
redu_m1 <- update(full_model, .~. - rest.age.std - land.use.hist - lui.std 
                  -c.perc.std - sand.perc.std - ph.value.std)
drop1(redu_m1, test = "Chi")
drop1(glmm_rest_site_1, test = "Chi")


#______________________________

## backward stepwise regression ####
bw_model <- step(glmm_rest_site_1, direction = "backward")

# with CWM F:
# glm.nb
# Step:  AIC=682.37
# target.hill.0 ~ rest.meth + site.cwm.abu.oek.f.std + ph.value.std + region

# with hydrology:
# Step:  AIC=662.06
# target.hill.0 ~ rest.meth + rest.age.std + land.use.hist + hydrology + 
#   ph.value.std + sand.perc.std + region + rest.meth:rest.age.std


## GLMM ####
full_model_glmm <- glmer(target.hill.0 ~ 
                           rest.meth
                         + rest.age.std
                         + land.use.hist
                         + hydrology
                         # + site.cwm.abu.oek.f.std
                         + lui.std
                         # + mngm.type
                         # + c.n.std
                         + c.perc.std
                         # + n.perc.std
                         + ph.value.std
                         # + tic.perc.std
                         # + toc.perc.std
                         # + clay.perc.std
                         # + silt.perc.std
                         + sand.perc.std
                         + rest.meth * rest.age.std
                         + (1|region),
                         data = model_data, family = "poisson")
# warum poisson? kommt von Felix script
summary(full_model_glmm)

# Auch bei GLMMs sollte die overdispersion getestet werden
# Das geht mit der folgenden Funktion

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

overdisp_fun(full_model_glmm)
# Ein signifikantes Ergebnis deutet auf overdispersion hin,
# also auch hier wieder der Umstieg auf die negative Binomial-verteilung

full_model_glmm2 <- glmer.nb(target.hill.0 ~
                               rest.meth
                             + rest.age.std
                             + land.use.hist
                             # + hydrology
                             # + site.cwm.abu.oek.f.std
                             + lui.std
                             # + mngm.type
                             # + c.n.std
                             + c.perc.std
                             # + n.perc.std
                             + ph.value.std
                             # + tic.perc.std
                             # + toc.perc.std
                             # + clay.perc.std
                             # + silt.perc.std
                             + sand.perc.std
                             + rest.meth * rest.age.std
                             + (1|region)
                             + (1|hydrology),
                             data = model_data)
overdisp_fun(full_model_glmm2)
# Und auch die overdispersion ist nicht mehr vorhanden danke der neuen
# Fehlerverteilung (negativ-binomial)

anova(full_model_glmm, full_model_glmm2)

# all subsets regression
# with MuMIn package
options(na.action = "na.fail") # Required for dredge to run

model_dredge <- dredge(full_model_glmm, beta = "none", evaluate = T, trace = 2,
                       # fixed = "rest.meth",
                       # m.lim =c(0,5),
                       rank = AIC) # when do you use AICc?
top_model <- get.models(model_dredge, subset = 1)[[1]]
summary(top_model)

# target.hill.0 ~ hydrology + land.use.hist + ph.value.std + rest.age.std +  
#   rest.meth + sand.perc.std + (1 | region) + rest.age.std:rest.meth
# AIC 671.2


# backward stepwise regression
bw_model <- step(full_model_glmm, direction = "backward")




# Grafische Daten Exploration ---
ggplot(data, aes(sum.soil.prep, target.hill.0, color = hydrology)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm", method.args = list(family = poisson), se = F) +
  facet_wrap(~region)

# prediction table
newdat1 <- expand.grid(sum.soil.prep = seq(0,1, by = 0.01),
                       region = unique(data$region),
                       hydrology = unique(data$hydrology))
# expand.grid kombiniert alle angegebenen Werte einer Variable mit allen
# Werten der anderen. Hier also 101 Werte für soil-prep x 3 Werte region + 3 Werte
# hydrology = 909 Zeilen in der Tabelle

# Und für all diese Werte mache ich jetzt vorhersagen
newdat1$pred.target.hill.0 <- predict(glm3, newdata = newdat1, type = "response")

ggplot(data, aes(sum.soil.prep, target.hill.0, color = hydrology)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~region) +
  geom_line(aes(sum.soil.prep, pred_hill0), data = newdat1)

summary(glm3)
# Interpretation der Ergebnisse
# Die Artenzahl steigt von Nord nach Süd
# Bei frischen Standorten gibt es keinen oder einen leicht negativen
# Effekt der Bodenvorbereitung ???!!!
# Bei moist und dry gibt es einen deutlichen positiven Effekt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D  ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 - Homogeneity ----------------------------------------------------------------

## 2 - Normality ----------------------------------------------------------------

