#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 1: Restoration vs. Reference sites
# Response variable: ratio tot.hill q1/q0
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


### diversity data ####

diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20241119.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))

diversity <- diversity %>% 
  mutate(ratio.q1.q0 = tot.hill.1 / tot.hill.0)
         



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
# data %>% 
#   filter(site.type == "restored") %>% 
#   filter(is.na(rest.meth))
# M_ALT: restoration method still unclear (2024-10-22)



## set model data --------------------------------------------------------------


# join diversity data
data_all <- data %>%
  left_join(diversity, by = "id.site")


data_all <- data_all %>%
  dplyr::select(
    id.site,
    ratio.q1.q0,
    hydrology,
    region,
    site.type
  )







rm(list = setdiff(ls(), c("data_all", "overdisp_fun")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - DATA EXPLORATION ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

range(data_all$ratio.q1.q0)

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
# doesn't look like there are outliers

sort(data_all$ratio.q1.q0)

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
  select(id.site, ratio.q1.q0) %>% 
  identify_outliers(ratio.q1.q0)
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


ggplot(data_all, aes(x = site.type, y = ratio.q1.q0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Site type")
ggplot(data_all, aes(x = region, y = ratio.q1.q0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Region")
ggplot(data_all, aes(x = hydrology, y = ratio.q1.q0)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Hydrology")
# difference between site types
# difference between regions --> use as random factor
# no difference between hydrology --> check if random factor necessary



## f Interactions --------------------------------------------------------------

# --> no numerical covariates

# check for possible interactions between covariates with coplot
# coplot(ratio.q1.q0 ~ site.type | region * hydrology,
#        data = data_all,
#        ylab = "Total species richness",
#        xlab = "",
#        panel = function(x, y, ...) {
#          tmp <- lm(y ~ x, na.action = na.omit)
#          abline(tmp)
#          points(x, y) })


# check categorical coviariates

ggplot(data_all, aes(x = interaction(region, site.type), y = ratio.q1.q0))+ 
  geom_boxplot()
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$region,
                 response = data_all$ratio.q1.q0)
# --> unclear if interaction between region and site.type

ggplot(data_all, aes(x = interaction(hydrology, site.type), y = ratio.q1.q0))+ 
  geom_boxplot()
interaction.plot(x.factor = data_all$site.type, trace.factor = data_all$hydrology,
                 response = data_all$ratio.q1.q0)
# --> interaction between hydrology and site.type



# test interactions between covariates
# interaction term significant --> interaction

int_model <- glm(ratio.q1.q0 ~ region * site.type, data = data_all, family = "poisson")
anova(int_model)
# no interaction

int_model <- glm(ratio.q1.q0 ~ hydrology * site.type, data = data_all, family = "poisson")
anova(int_model)
# no interaction




## g Spatial dependency --------------------------------------------------------------


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
#' no interaction between site.type and region
#' interaction between site.type and hydrology unclear



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - RANDOM STRUCTURE ########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list = setdiff(ls(), c("data_all", "overdisp_fun")))

# # eliminate all rows with any missing values
data_model <- na.omit(data_all) 



### a Random structure ---------------------------------------------------------

R1 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region), data = data_model)
R2 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|hydrology), data = data_model)
R3 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region) + (1|hydrology), data = data_model)
R4 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region) + (1|region:hydrology), data = data_model)
R5 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|hydrology) + (1|region:hydrology), data = data_model)
# -> singularity
R6 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region) + (1|hydrology) + (1|region:hydrology), data = data_model)
# R7 <- glmer.nb(ratio.q1.q0 ~ 1 + (rest.age.std|region) + (1|hydrology), data = data_model)
# R8 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region) + (rest.age.std|hydrology), data = data_model)
# -> singularity
# R9 <- glmer.nb(ratio.q1.q0 ~ 1 + (rest.age.std|region) + (rest.age.std|hydrology), data = data_model)
# -> singularity
Rnull <- glm(ratio.q1.q0 ~ 1, data = data_model, family = "poisson") # right family??


AICc(R1, R2, R3, R4,
     # R5,
     R6, #R7,
     # R8, R9,
     Rnull) %>% 
  arrange(AICc)
# --> model 3 (region and hydrology as random factors, no interaction, no random slope)
# --> go on to build model with this random structure

# # is hydrology better as fixed factor?
# R10 <- glmer.nb(ratio.q1.q0 ~ 1 + (1|region) + hydrology, data = data_model)
# AICc(R3, R10)
# # --> we should consider this (R10 has lower AICc)
# # --> test this with full model later

# test if use of binomial-distribution is ok (due to overdispersion)
R3a <- glmer(ratio.q1.q0 ~ 1 + (1|region) + (1|hydrology), data = data_model, family = "poisson")
overdisp_fun(R3a)
# significant result shows overdispersion
# --> using of negative binomial-distribution is neccessary

# # how strong is hydrology?
# glmm_hydr <- glmer.nb(ratio.q1.q0 ~ hydrology + (1|region), data = data_model)
# glmm_test <- glmer.nb(ratio.q1.q0 ~ rest.meth + hydrology + (1|region), data = data_model)
# glmm_test_1 <- glmer.nb(ratio.q1.q0 ~ rest.meth : hydrology + (1|region), data = data_model)
# glmm_test_2 <- glmer.nb(ratio.q1.q0 ~ rest.meth * hydrology + (1|region), data = data_model)
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
# B1 <- glmer(ratio.q1.q0 ~ site.type
#             + (1|region) + (1|hydrology), data = data_model,
#             family = "poisson")
B1 <- glmmTMB(ratio.q1.q0 ~ site.type + (1|region) + (1|hydrology), 
              data = data_model,
              family = beta_family(link = "logit")
)
check_overdispersion(B1)
# significant -> overdispersion

# # change to non-binomial
# # B1 <- glmer.nb(ratio.q1.q0 ~ site.type
# #                + (1|region) + (1|hydrology), data = data_model)
# B1 <- glmmTMB(ratio.q1.q0 ~ site.type + (1|region) + (1|hydrology), 
#               data = data_model,
#               family = nbinom2                      # Negative binomial family
# )

summary(B1)
performance(B1)
#' ICC = 0.365
#' The covariates explain 28 % of the variation in richness (R2 marg.)
#' The covariates and random effects explain 54 % of the variation. (R2 cond.)


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
# ratio.q1.q0 ~ site.type + (1 |region) + (1 |hydrology)

data_model <- data_all

B1_final <- glmmTMB(ratio.q1.q0 ~ site.type + (1|region) + (1|hydrology), 
                    data = data_model,
                    family = beta_family(link = "logit"))

summary(B1_final)
#                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        3.13399    0.13443  23.312   <2e-16 ***
# site.typerestored  0.53811    0.06168   8.724   <2e-16 ***
# site.typepositive  0.62580    0.07465   8.383   <2e-16 ***

# --> site.type is significant


### b Post-hoc test ----

# Tukey-adjusted pariwise comparisons
# generate estimated marginal means (EMMs) and then apply a Tukey correction to 
# pairwise comparisons
emm.site.type <- emmeans(B1_final, "site.type")
summary(emm.site.type, infer = F, type = "response")
pairs(emm.site.type, adjust = "tukey", 
      # type = "response"
)
# contrast            estimate     SE  df z.ratio p.value
# negative - restored  -0.5381 0.0617 Inf  -8.724  <.0001
# negative - positive  -0.6258 0.0746 Inf  -8.383  <.0001
# restored - positive  -0.0877 0.0557 Inf  -1.574  0.2570
# 
# Results are given on the log (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 3 estimates 

# pairs(regrid(emm.site.type), adjust = "tukey") # regrid() for back-transformation from log-scale
# contrast            estimate   SE  df z.ratio p.value
# negative - restored    -16.4 2.58 Inf  -6.335  <.0001
# negative - positive    -20.0 3.46 Inf  -5.772  <.0001
# restored - positive     -3.6 2.39 Inf  -1.510  0.2861
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

## --> different p-values, because of calculating them after back-transformation


### c report final model ----


# Tidy up the model summary
model_summary <- broom.mixed::tidy(B1_final, effects = "fixed") %>%
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
      "outputs", "tables", "vegetation", "model_summary_restref_hill0_B1_final.csv"
    )
  )

### d save final model ----


save(B1_final, file = here("outputs", "models", "vegetation", "model_plants_restref_q1q0ratio_B1_final.Rdata"))

# load(file = here("outputs", "models", "vegetation", "model_plants_restref_q1q0ratio_B1_final.Rdata"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E - VISUALIZATION  ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### summary statistics ###
data_model %>%
  group_by(site.type) %>%
  get_summary_stats(ratio.q1.q0, type = "full")

# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.site.type <- emmeans(B1_final, ~ site.type)
emm.df <- summary(emm.site.type, infer = F, type = "response") # type ="response" for back-transformation from log-scale

# compact letter display
cld_results <- multcomp::cld(emm.site.type, adjust = "tukey", Letters = letters)

# plot with EMMs and SE
ggplot(data_model, aes(x = site.type, y = ratio.q1.q0)) +
  geom_violin(fill = "lightblue") +
  # geom_quasirandom() +
  # geom_jitter2(width = 0.05, alpha = 0.5) +
  # geom_line(data = means, aes(y = Mean, group = 1), size = 1) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 1,
    color = "black"
  ) +
  geom_text(data = cld_results, aes(x = site.type, y = 0.75, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4) +
  # stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "red") +
  theme_minimal()

ggsave("outputs/figures/plants_species_diversity/model_restref_q1q0ratio_sitetype_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")




## end script








