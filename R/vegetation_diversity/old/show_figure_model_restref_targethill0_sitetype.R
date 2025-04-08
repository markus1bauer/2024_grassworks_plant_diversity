#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Plot figure: Model restref_targethill0 B1 - Site type
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

load(file = here("outputs", "models", "vegetation", "model_plants_restref_targethill0_B1_final.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - Plot  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.site.type <- emmeans(B1_final, ~ site.type)
emm.df <- summary(emm.site.type, infer = F, type = "response") # type ="response" for back-transformation from log-scale
# site.type response   SE  df
# negative      17.7 2.61 Inf
# restored      32.4 4.49 Inf
# positive      37.0 5.33 Inf

# compact letter display
cld_results <- multcomp::cld(emm.site.type, adjust = "tukey", Letters = letters)

# plot with EMMs and SE
plot2 <- ggplot(data_model, aes(x = site.type, y = target.hill.0)) +
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
  geom_text(data = cld_results, aes(x = site.type, y = 100, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4) +
  labs(x = NULL, y = "Characteristic species richness (16 m²)") +
  scale_y_continuous(limits = c(0, 100)) +
  # stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "red") +
  theme_minimal()
plot2

library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)
ggsave("outputs/figures/plants_species_diversity/model_restref_sitetype_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave("outputs/figures/plants_species_diversity/model_restref_targethill0_sitetype_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")



