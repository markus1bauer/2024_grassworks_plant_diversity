#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# # Plot figure: Model restfact_targethill0 - Restoration method (german, boxplot)
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

### Colours ###
my_col <- c("#4E79A7", "#F28E2B" ,"#E15759", "#59A14F")


### Labels ###
gglayer_labs <- list(
  scale_x_discrete(
    labels = c(
      "Cultivar\nSeed Mixture",
      "Management\nAdaptation",
      "Regional\nSeed Mixture",
      "Direct\nHarvesting"),
  )
)

### Fonts ###
windowsFonts("Franklin Gothic Book" = windowsFont("Franklin Gothic Book"))
windowsFonts("Trebuchet MS" = windowsFont("Trebuchet MS"))

### Theme ###
gglayer_theme <- list(
  theme_minimal(
    base_size = 10, 
    base_family = "Franklin Gothic Book"
  ),
  theme(
    legend.position = "none",
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = rel(1)),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = rel(1.2)),
    axis.line = element_line(color = "grey"),
    panel.grid.major =  element_line(linetype = 1),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      family = "Trebuchet MS",
      size = rel(1.2),
      face = "bold"
    ),
    plot.caption = element_text(size = rel(0.5))
  ),
  scale_fill_manual(values = my_col)
)


## load data -------------------------------------------------------------------

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0.Rdata"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - Plot  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_targethill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response   SE  df
# cus           26.5 4.06 Inf
# mga           29.4 4.46 Inf
# res           33.3 4.90 Inf
# dih           36.6 5.29 Inf

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)


(p4 <- ggplot(data_model_target, aes(x = rest.meth, y = target.hill.0, fill = rest.meth)) +
  gglayer_theme +
  scale_x_discrete(
      labels = c(
        "Regelsaatgut",
        "Management",
        "Regiosaatgut",
        "Direkternte"),
    ) +
  geom_boxplot() +
  # geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
  # geom_pointrange(
  #   data = emm.df,
  #   aes(y = response, ymin = response - SE, ymax = response + SE),
  #   size = 0.7, linewidth = 0.8
  # ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 71, label = .group),          # Add compact letters
            vjust = 0, hjust = 0.7, size = rel(4),
            family = "Franklin Gothic Book", color = "black", 
            # fontface = "bold"
  ) +
  labs(
    title = "Grünland-typische Pflanzenarten",
    y = "Artenzahl",
    # caption = "Tukey adjusted pair-wise comparisons"
    # tag = "A"
  )) +
  coord_cartesian(ylim = c(0, 72))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave("outputs/figures/plants_species_diversity/targethill0_restmeth_boxplot_de.jpg",
       dpi = 300, width = 12, height = 12, units = "cm")
