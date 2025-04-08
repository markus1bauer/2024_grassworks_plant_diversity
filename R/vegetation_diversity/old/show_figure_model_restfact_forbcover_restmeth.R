#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Plot figure: Model restfact_forbcover - Restoration method
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

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_forbcover.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - Plot  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_forbcover, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response     SE  df
# cus          0.298 0.0674 Inf
# mga          0.247 0.0585 Inf
# res          0.321 0.0642 Inf
# dih          0.374 0.0668 Inf
# 
# Results are averaged over the levels of: land.use.hist 

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)

# plot with EMMs and SE
ggplot(data_model_forbcover, aes(x = rest.meth, y = forb.cover)) +
  geom_violin(fill = "grey", trim = F) +
  # geom_quasirandom() +
  # geom_jitter2(width = 0.05, alpha = 0.5) +
  # geom_line(data = means, aes(y = Mean, group = 1), size = 1) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 1,
    color = "black"
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 1, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4) +
  # stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "red") +
  theme_minimal()



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave("outputs/figures/plants_species_diversity/model_restfact_forbcover_restmeth_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")

