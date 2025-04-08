#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: forb-index
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A PREPARATION ###############################################################
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
library(broom.mixed) # Tidy up the model summary
library(scales) # label_number()
library(glmmTMB)

### Start ###
rm(list = ls())


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
    id.site, site.type, hydrology, region, 
    longitude, latitude
  ) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         site.type = fct_relevel(site.type, "negative", "restored", "positive"),
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

# # Standardizing the numeric explanatory variables
# data <- data %>%
#   mutate(across(where(is.character), as.factor)) %>%
#   mutate(across(where(is.numeric), ~ as.numeric(scale(.)), .names = "{col}.std"))

# # Verify scaling
# summary(data)
# 
# colSums(is.na(data))
# data %>% 
#   filter(site.type == "restored") %>% 
#   filter(is.na(rest.meth))
# M_ALT: restoration method still unclear (2024-10-22)



## set model data --------------------------------------------------------------

data <- sites %>% 
  left_join(forb_index, by = "id.site")

data_all <- data %>%
  dplyr::select(
    id.site,
    forb.index.site,
    hydrology,
    region,
    site.type
  ) %>% 
  rename(forb.index = forb.index.site)

# remove zeros for gamma distribution (doesn't allow zeros)
data_all <- data_all %>% 
  mutate(forb.index = forb.index + 0.001)

# remove outliers
## see data exploration (b)
data_all <- data_all %>% 
  filter(forb.index < 20) %>% 
  filter(id.site != "M_NIS")





rm(list = setdiff(ls(), c("data_all", "overdisp_fun", "sites")))


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
# there are outliers

sort(data_all$forb.index)

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
  select(id.site, forb.index) %>% 
  identify_outliers(forb.index)
# extreme outliers with forb.index > 20

# outliers within groups
data_all %>% 
  select(id.site, forb.index, site.type) %>% 
  group_by(site.type) %>% 
  identify_outliers(forb.index)
# M_NIS

# remove outliers
data_all <- data_all %>% 
  filter(forb.index < 20) %>% 
  filter(id.site != "M_NIS")

data_all %>% 
  select(id.site, forb.index) %>% 
  identify_outliers(forb.index)
# no extreme outliers


## c inspect categorical covariates -----------------------------------------

table(data_all$site.type)
#' Unbalanced...but enough observations per level.

table(data_all$hydrology)
#' Unbalanced...but enough observations per level.

table(data_all$region)
#' Balanced.

#' Was each sitetype measured in every hydrology?
table(data_all$site.type, data_all$hydrology)
histogram( ~ sitetype | hydrology, data_all)
#' Unbalanced, do we have enough observations per combination?
#' only 4 observation in dry-negative

#' Was each sitetype measured in every region?
table(data_all$site.type, data_all$region)
histogram( ~ sitetype | region, data_all)
#' Unbalanced, but enough observations




## d Check collinearity part 1 (Step 5) ----------------------------------------

# between continuous covariates

# no numerical variable in model data --> no need to check

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

# ggplot(data_all, aes(x = rest.meth, y = rest.age)) +
#   geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
# ggplot(data_all, aes(x = land.use.hist, y = rest.age)) +
#   geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
# ggplot(data_all, aes(x = hydrology, y = rest.age)) +
#   geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")
# ggplot(data_all, aes(x = region, y = rest.age)) +
#   geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent")

#' restoration method is collinear with Age; the boxplots are not next to each other.
#' region is unbalanced, but fine I guess.
#' The rest is fine.



# source("HighstatLibV14.R")
# MyVar <- c("rest.meth", "rest.age", "land.use.hist", "hydrology", "region")
# corvif(data_all[,MyVar])


## e Relationships --------------------------------------------------------------

#' Plot response variable versus each covariate.


ggplot(data_all, aes(x = site.type, y = forb.index)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Site type")
ggplot(data_all, aes(x = region, y = forb.index)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = forb.index)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# difference between site types
# difference between regions not very clear --> test as random factor
# difference between hydrology --> use as random factor

## f distribution --------------------------------------------------------------
library(lattice)
histogram(data_all$forb.index)
# gamma?

x <- data_all$forb.index

library(fitdistrplus)
library(logspline)

descdist(x, discrete = FALSE)

fit.gamma <- fitdist(x, "gamma")
fit.beta <- fitdist(x, "beta")


plot(fit.lnorm)
plot(fit.gamma)
plot(fit.norm)
# gamma looks okeyish

fit.lnorm$aic
fit.gamma$aic
fit.norm$aic
# gamma has lowes AIC

library(gamlss)
library(gamlss.dist)
library(gamlss.add)

fit <- fitDist(x, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)

summary(fit)

## f Interactions --------------------------------------------------------------

# --> no numerical covariates

# check for possible interactions between covariates with coplot
# coplot(forb.index ~ site.type | region * hydrology,
#        data = data_all,
#        ylab = "Total species richness",
#        xlab = "",
#        panel = function(x, y, ...) {
#          tmp <- lm(y ~ x, na.action = na.omit)
#          abline(tmp)
#          points(x, y) })


# check categorical coviariates

ggplot(data_all, aes(x = interaction(region, site.type), y = forb.index))+ 
  geom_boxplot()
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$forb.index)
# interaction between region and site.type

ggplot(data_all, aes(x = interaction(hydrology, site.type), y = forb.index))+ 
  geom_boxplot()
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$forb.index)
# --> interaction between hydrology and site.type

ggplot(data_all, aes(x = interaction(hydrology, region), y = forb.index))+ 
  geom_boxplot()
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$forb.index)
# interaction between hydrology and region


# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

int_model <- glm(forb.index ~ region * site.type, data = data_all, family = Gamma(link = "log"))
anova(int_model)
# no interaction

int_model <- glm(forb.index ~ hydrology * site.type, data = data_all, family =  Gamma(link = "log"))
anova(int_model)
# no interaction

int_model <- glm(forb.index ~ hydrology * region, data = data_all, family =  Gamma(link = "log"))
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


## h conclusions  --------------------------------------------------------------

#' no missing values
#' no outliers
#' very little observations in negative-dry
#' interactions unclear



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all", "overdisp_fun")))

# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 



### a Random structure ---------------------------------------------------------

R1 <- glmmTMB(forb.index ~ 1 + (1|region), data = data_model, family = Gamma(link="log"))
R2 <- glmmTMB(forb.index ~ 1 + (1|hydrology), data = data_model, family = Gamma(link="log"))
R3 <- glmmTMB(forb.index ~ 1 + (1|region) + (1|hydrology), data = data_model, family = Gamma(link="log"))
R4 <- glmmTMB(forb.index ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R5 <- glmmTMB(forb.index ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
R6 <- glmmTMB(forb.index ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link="log"))
Rnull <- glm(forb.index ~ 1, data = data_model, family = Gamma(link="log")) # right family??


AIC(R1, R2, R3, R4, R5, R6,
    Rnull
) %>% 
  arrange(AIC)
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(forb.index ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later


# # how strong is hydrology?
# glmm_hydr <- glmer.nb(forb.index ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(forb.index ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(forb.index ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(forb.index ~ rest.meth * hydrology + (1|region), data = data_model)
# summary(glmm_hydr)
# summary(glmm_test)
# summary(glmm_test_1)
# summary(glmm_test_2)
# # rest.meth is not significant anymore --> hydrology covers all effects



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
# B1 <- glmer(forb.index ~ site.type
#             + (1|region) + (1|hydrology), data = data_model,
#             family = "poisson")
B1 <- glmmTMB(forb.index ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = Gamma(link="log")
)
check_overdispersion(B1)
#' No overdispersion detected.


summary(B1)




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
#' smaller variance in residuals with larger values?



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



### d Model validation with DHARMa ---------------------------------------------

# ### a Plot residuals vs fitted values
# simulation_output <- simulateResiduals(B1, plot = TRUE)
# # quantile deviations detected
# 
# 
# ### b Plot residuals vs covariates in the model
# plotResiduals(simulation_output$scaledResiduals, data_model$site.type) # ok
# plotResiduals(simulation_output$scaledResiduals, data_model$hydrology)
# # within-group deviations from uniformity significant
# plotResiduals(simulation_output$scaledResiduals, data_model$region)
# # within-group deviations from uniformity significant
# 
# ### c Plot residuals vs covariates not in the model
# # no other covariates


### e check for adjustments of model -------------------------------------------

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
# forb.index ~ site.type + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restref_forbindex.Rdata"))

data_model_forbindex <- data_all

restref_forbindex <- glmmTMB(forb.index ~ site.type + (1|region) + (1|hydrology), 
                           data = data_model_forbindex,
                           family = Gamma(link="log")
)

summary(restref_forbindex)
# Estimate Std. Error z value Pr(>|z|)  
# (Intercept)        0.05662    0.24180   0.234    0.815    
# site.typerestored  1.40656    0.17624   7.981 1.45e-15 ***
#   site.typepositive  1.62048    0.21656   7.483 7.27e-14 ***

# --> site.type is significant

performance(restref_forbindex)
#' ICC = 0.155
#' The covariates explain 32 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 43 % of the variation. (R2 cond.)


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.site.type <- emmeans(restref_forbindex, "site.type")
summary(emm.site.type, infer = F, type = "response")
pairs(emm.site.type, adjust = "tukey", 
      # type = "response"
)
# contrast            estimate    SE  df z.ratio p.value
# negative - restored   -1.407 0.176 Inf  -7.981  <.0001
# negative - positive   -1.620 0.217 Inf  -7.483  <.0001
# restored - positive   -0.214 0.174 Inf  -1.232  0.4339
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 3 estimates 

pairs(regrid(emm.site.type), adjust = "tukey") # regrid() for back-transformation from log-scale

## --> different p-values, because of calculating them after back-transformation


### c report final model ----


# Tidy up the model summary
model_summary <- broom.mixed::tidy(restref_forbindex, effects = "fixed") %>%
  dplyr::select(term, estimate, std.error, statistic, p.value) %>%
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
      "outputs", "tables", "vegetation", "model_summary_restref_forbindex.csv"
    )
  )

### summary statistics ###
data_model_forbindex %>%
  group_by(site.type) %>%
  get_summary_stats(forb.index, type = "full")


### d save final model ----


save(restref_forbindex, data_model_forbindex,
     file = here("outputs", "models", "vegetation",
                 "model_plants_restref_forbindex.Rdata"))






## end script








