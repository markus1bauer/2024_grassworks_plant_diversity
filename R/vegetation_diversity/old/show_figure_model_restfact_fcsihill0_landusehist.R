#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Plot figure: Model restfact_fcsihill0 - Previous landuse
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(emmeans) # calculate estimated marginal means and post-hoc Tukey



### Start ###
rm(list = ls())


## load data -------------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - Plot  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.land.use.hist <- emmeans(restfact_fcsihill0, ~ land.use.hist)
emm.df <- summary(emm.rest.land.use.hist, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# land.use.hist emmean    SE  df
# arable land     3.36 0.142 112
# grassland       3.47 0.145 112
# 
# Results are averaged over the levels of: rest.meth 

# compact letter display
cld_results <- multcomp::cld(emm.land.use.hist, adjust = "tukey", Letters = letters)

# plot with EMMs and SE
ggplot(data_model_fcsihill0, aes(x = land.use.hist, y = fcsi.hill.0)) +
  geom_violin(fill = "grey", trim = T) +
  # geom_quasirandom() +
  # geom_jitter2(width = 0.05, alpha = 0.5) +
  # geom_line(data = means, aes(y = Mean, group = 1), size = 1) +
  geom_pointrange(
    data = emm.df,
    aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE),
    size = 1,
    color = "black"
  ) +
  geom_text(data = cld_results, aes(x = land.use.hist, y = 4.5, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4) +
  # stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "red") +
  theme_minimal()



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave("outputs/figures/plants_species_diversity/model_restfact_fcsihill0_landusehist_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")

