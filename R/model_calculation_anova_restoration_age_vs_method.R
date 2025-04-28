#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Collinearity restoration method and restoration age - ANOVA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A - PREPARATIION ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)


### Start ###
rm(list = ls())



## load data -------------------------------------------------------------------

### site environment data ####

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20241114.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, rest.meth, rest.age,
    ) %>%
  distinct() %>% 
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih")
  )


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
sites %>% 
  select(id.site, rest.meth, rest.age) %>% 
  group_by(rest.meth) %>%
  identify_outliers(rest.age)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(rest.age ~ rest.meth, data = sites, na.action = na.omit)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, as all the points fall approximately along the reference line,
# we can assume normality. This conclusion is supported by the Shapiro-Wilk test.
# The p-value is not significant (p = 0.191), so we can assume normality.

# Check normality assumption by groups
sites %>%
  filter(!is.na(rest.meth)) %>% 
  group_by(rest.meth) %>%
  shapiro_test(rest.age)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(sites, "rest.age", facet.by = "rest.meth")
# looks okayish

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# unlear if homogeneity is true

# Levene's test #
sites %>% levene_test(rest.age ~ rest.meth)
# From the output above, we can see that the p-value is < 0.05, which is 
# significant. This means that, there is significant difference between 
# variances across groups. Therefore, we can not assume the homogeneity of variances
# in the different treatment groups.

# --> maybe change to glm ?

# model_glm  <- glm(rest.age ~ rest.meth, data = sites, na.action = na.omit)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

res_aov <- anova_test(model)
# res_aov <- sites %>% 
#   anova_test(rest.age ~ rest.meth)
res_aov
# ANOVA Table (type II tests)
# 
#       Effect DFn DFd      F        p p<.05   ges
# 1 rest.meth   3 104 19.129 5.93e-10     * 0.356

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(3, 104) = 19.129, p < 0.001, eta2[g] = 0.356.

sites_20y <- sites %>% 
  filter(rest.age <= 20)
model_20y  <- lm(rest.age ~ rest.meth, data = sites_20y, na.action = na.omit)
res_aov_20y <- anova_test(model_20y)
# res_aov_20y <- sites_20y %>% 
#   anova_test(rest.age ~ rest.meth)
res_aov_20y
# ANOVA Table (type II tests)
# 
#       Effect DFn DFd    F     p p<.05   ges
# 1 rest.meth   3  87 2.13 0.102       0.068

# From the above ANOVA table, it can be seen that there are no significant 
# differences between groups (p > 0.001), 
# F(3, 87) = 2.13, p = 0.102, eta2[g] = 0.068.



