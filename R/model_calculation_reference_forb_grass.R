#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: forb-grass-ratio
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
    id.site,
    site.type, hydrology, region,
    site.cover.grass, site.cover.forbs, site.cover.legumes,
    site.forb.legu.grass.ratio, 
    longitude, latitude
  ) %>%
  distinct() %>% 
  mutate(region = fct_relevel(region, "north", "centre", "south"),
         hydrology = fct_relevel(hydrology, "dry", "fresh", "moist"),
         site.type = fct_relevel(site.type, "negative", "restored", "positive"),
  )






## transform input data --------------------------------------------------------

# standardise explanatory variable (only numerical variables)

# --> no numerical explanatory variabls



## set model data --------------------------------------------------------------

data_all <- sites %>%
  dplyr::select(
    id.site,
    site.forb.legu.grass.ratio,
    hydrology,
    region,
    site.type
  ) %>% 
  rename(fg.ratio = site.forb.legu.grass.ratio)




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
dotchart(data_all$fg.ratio,
         ylab = "Order of the data")
# outlier are points far right or left in plot
# doesn't look like there are outliers


## wie umgehen mit ouliers?? Messfehler: unrealistische WErte
# sort(data_all$fg.ratio)


# another test for outliers
data_all %>% 
  select(id.site, fg.ratio) %>% 
  identify_outliers(fg.ratio)
# extreme outliers with fg.ratio >= 3.96

sites %>% 
  filter(id.site %in% c("S_JAU", "S_NOZ", "M_JER", "S_MCH"))
sort(sites$site.cover.legumes)
# nothing strange with outliers --> keep them



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


ggplot(data_all, aes(x = site.type, y = fg.ratio)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Site type")
ggplot(data_all, aes(x = region, y = fg.ratio)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = fg.ratio)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# difference between site types
# difference between regions --> use as random factor
# difference between hydrology --> use as random factor


## f distribution --------------------------------------------------------------

library(lattice)
histogram(data_all$fg.ratio)
# always positive, right skewed with long tail
# gamma distribution? exponential?
histogram(log(data_all$fg.ratio))
# doesn't look like a normal distribution --> not lognormal

x <- data_all$fg.ratio

library(fitdistrplus)
library(logspline)

descdist(x, discrete = FALSE)

fit.norm <- fitdist(x, "norm")
fit.exp <- fitdist(x, "exp")
fit.gamma <- fitdist(x, "gamma")
fit.lnorm <- fitdist(x, "lnorm")
print(fit.gamma)
# shape is close to 1 --> exponential distribution could be possible

# If the variance increases with the mean, 
# this suggests a Gamma (or exponential) distribution might be a good fit
df_summary <- data_all %>%
  group_by(site.type) %>%
  summarize(mean_ratio = mean(fg.ratio), 
            var_ratio = var(fg.ratio))
ggplot(df_summary, aes(x = mean_ratio, y = var_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Variance vs. Mean of Forb-Grass Ratio")
# increases
plot(fit.norm)
plot(fit.exp)
plot(fit.gamma)
plot(fit.lnorm)
# gamma and exponential looks okay

fit.norm$aic
fit.exp$aic
fit.gamma$aic
fit.lnorm$aic
# exponential has lowest AIC, but < delta 2 to gamma
# exponential is a special case of gamma, glmmTMB uses gamma with link log for both
# --> use gamma distribution and check shape = 1 if necessary

detach(package:fitdistrplus)


## g Interactions --------------------------------------------------------------

# need to be checked again (2025-03-05)

# --> no numerical covariates


# check categorical coviariates

## site.type vs. X

# region
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(region, site.type), y = fg.ratio))+ 
  geom_boxplot()
# --> no clear interaction

# hydrology
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(hydrology, site.type), y = fg.ratio))+ 
  geom_boxplot()
# --> interaction


## region vs. X

# hydrology
interaction.plot(x.factor = data_all$region, trace.factor = data_all$hydrology,
                 response = data_all$fg.ratio)
ggplot(data_all, aes(x = interaction(hydrology, region), y = fg.ratio))+ 
  geom_boxplot()
# --> no interaction


# test interactions between covariates
# interaction term significant --> interaction
library(MASS)

## site.type vs. X

# region
int_model <- glm(fg.ratio ~ region * site.type, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

# hydrology
int_model <- glm(fg.ratio ~ hydrology * site.type, data = data_all, family = Gamma(link = "log"))
check_overdispersion(int_model)
anova(int_model)
# no interaction

## region vs. X

# hydrology
int_model <- glm(fg.ratio ~ hydrology * region, data = data_all, family = Gamma(link = "log"))
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
#' no interaction between:
#'  - site.type and region
#'  - region and hydrology
#' interaction unclear between:
#'  - site.type and hydrology


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all")))

# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 



### a Random structure ---------------------------------------------------------

R1 <- glmmTMB(fg.ratio ~ 1 + (1|region), data = data_model, family = Gamma(link = "log"))
R2 <- glmmTMB(fg.ratio ~ 1 + (1|hydrology), data = data_model, family = Gamma(link = "log"))
R3 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|hydrology), data = data_model, family = Gamma(link = "log"))
R4 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|region:hydrology), data = data_model, family = Gamma(link = "log"))
R5 <- glmmTMB(fg.ratio ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link = "log"))
R6 <- glmmTMB(fg.ratio ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model, family = Gamma(link = "log"))
Rnull <- glm(fg.ratio ~ 1, data = data_model, family = Gamma(link = "log")) # right family??


AIC(R1, R2, R3, R4, R5, R6, Rnull) %>% 
  arrange(AIC)
# --> model 2 is the same good as model 3, and simpler. but we keep going on with model 3
# because of study design structure
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
B1 <- glmmTMB(fg.ratio ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = Gamma(link = "log")
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
# looks ok


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
# fg.ratio ~ site.type + (1 |region) + (1 |hydrology)

# load(file = here("outputs", "models", "vegetation", "model_plants_restref_fgratio.Rdata"))

data_model_fgratio <- data_all

restref_fgratio <- glmmTMB(fg.ratio ~ site.type + (1|region) + (1|hydrology), 
                               data = data_model_fgratio,
                               family = Gamma(link = "log")
)

summary(restref_fgratio)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -1.2526     0.2522  -4.968 6.78e-07 ***
#   site.typerestored   1.2516     0.1717   7.291 3.07e-13 ***
#   site.typepositive   1.2255     0.2105   5.821 5.86e-09 ***

# --> site.type is significant

performance(restref_fgratio)
#' ICC = 0.186
#' The covariates explain 26 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 40 % of the variation. (R2 cond.)

### summary statistics ###
data_model_fgratio %>%
  group_by(site.type) %>%
  get_summary_stats(fg.ratio, type = "full")


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.site.type <- emmeans(restref_fgratio, "site.type")
summary(emm.site.type, infer = F, type = "response")
# site.type response     SE  df
# negative     0.286 0.0721 Inf
# restored     0.999 0.2178 Inf
# positive     0.973 0.2425 Inf
pairs(emm.site.type, adjust = "tukey", 
      # type = "response"
)
# contrast            estimate    SE  df z.ratio p.value
# negative - restored   -1.252 0.172 Inf  -7.291  <.0001
# negative - positive   -1.226 0.211 Inf  -5.821  <.0001
# restored - positive    0.026 0.171 Inf   0.152  0.9873
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 3 estimates  

pairs(regrid(emm.site.type), adjust = "tukey") # regrid() for back-transformation from log-scale
## --> different p-values, because of calculating them after back-transformation

# --> restored are significantly different from neg ref but not from pos ref


### c report final model ----

# see skript show_table_model_results_restref.R



### d save final model ----


save(restref_fgratio, data_model_fgratio,
     file = here("outputs", "models",
                 "model_reference_forb_grass.Rdata"))






## end script








