#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Plot figure: Model restfact_tothill0 B1 - Restoration method
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

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B - Plot  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_tothill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response   SE  df
# cus           32.3 4.46 Inf
# mga           36.4 4.96 Inf
# res           38.6 5.08 Inf
# dih           45.0 5.79 Inf

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters)

# sample size
sample_size = data_model_tot %>%
  group_by(rest.meth) %>%
  summarize(num = n())



# plot with EMMs and SE
# palette("Tableau 10")
# palette()  
my_col <- c("#4E79A7", "#F28E2B" ,"#E15759", "#59A14F")
# 
# my_theme <- theme_minimal()

gglayer_labs <- list(
  scale_x_discrete(
    labels = c(
    "Cultivar\nSeed Mixture",
    "Regional\nSeed Mixture",
    "Management\nAdaptation",
    "Direct\nHarvesting"),
    )
)

library(systemfonts)
fonts <- system_fonts()
# See all monospace fonts
fonts <- system_fonts()
fonts[fonts$monospace, ]

fonts %>% 
  filter(width == "condensed",
         italic == F,
         monospace == F,
         weight == "normal") %>% 
  select(name, family, style, weight) %>% 
  arrange(family) %>% 
  print(n= 25)

windowsFonts("Arial Narrow" = windowsFont("Arial Narrow"))
windowsFonts("Franklin Gothic Book" = windowsFont("Franklin Gothic Book"))
windowsFonts("Trebuchet MS" = windowsFont("Trebuchet MS"))


gglayer_theme <- list(
  theme_minimal(
    base_size = 20, 
    # base_family = "Arial Narrow"
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
    # plot.margin = margin(t = rel(1.5)),
    plot.caption = element_text(size = rel(0.5))
  ),
  scale_fill_manual(values = my_col)
)



plot_A <- ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth,)) +
  gglayer_theme +
  gglayer_labs +
  geom_violin(trim = T, width = 0.8, linewidth = 0.7) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 1.5, linewidth = 1
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 85, label = .group),          # Add compact letters
            vjust = 0, hjust = 0.7, size = rel(7),
            family = "Franklin Gothic Book", color = "black", 
            # fontface = "bold"
            ) +
  labs(
    title = "(A) Total Species Richness",
    y = "Species Number (per 16 m²)",
    # caption = "Tukey adjusted pair-wise comparisons"
    # tag = "A"
  )


ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth,)) +
  geom_violin(trim = T) +
  # scale_fill_manual(values = my_col) +
  # geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # geom_quasirandom() +
  # geom_jitter2(width = 0.05, alpha = 0.5) +
  # geom_line(data = means, aes(y = Mean, group = 1), size = 1) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 1,
    color = "black"
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 80, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4) +
  # stat_summary(fun = "mean", geom = "point", shape = 18, size = 3, color = "red") +
  # scale_x_discrete(labels = c(
  #   "Cultivar Seed Mixture",
  #   "Regional Seed Mixture",
  #   "Management Adaptation",
  #   "Direct Harvesting"))+
  
  # labs(x = "Restoratin method", y = "Effective number of species")+
  theme_mb()
  # theme_minimal()
p1
p1 + theme_light()

p2 <- ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth,)) +
  geom_violin(trim = T)
library(ggthemes)
p2 + theme_light() 

library(ggpubr)
p2 + theme_pubclean() 

p2 + scale_colour_colorblind()

dsamp <- diamonds[sample(nrow(diamonds), 1000), ]
p <- qplot(carat, price, data=dsamp, colour=clarity) + theme_igray()
p + scale_colour_colorblind()
p
palette("default")
p1 + 
  scale_x_discrete(labels = c(
    "Cultivar Seed Mixture",
    "Regional Seed Mixture",
    "Management Adaptation",
    "Direct Harvesting"))
  

p <- ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth,)) +
  scale_fill_manual(values = my_col)
p


p1 +
  geom_violin(trim = T) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 1,
    color = "black"
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 80, label = .group),          # Add compact letters
            vjust = -0.5, color = "black", size = 4)

library(hrbrthemes)
ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth,)) +
  scale_fill_manual(values = my_col) +
  theme_ipsum()

  

plot_tothill0

library(nationalparkcolors)
names(park_palettes)
pal <- park_palette("Badlands")


ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth)) +
  geom_violin(trim = T) +
  scale_fill_manual(values = my_col)
ggplot(data_model_tot, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth)) +
  geom_violin(trim = T) +
  scale_fill_manual(values = pal)



library(ggstatsplot)
plt <- ggbetweenstats(data = data_model_tot, x = rest.meth, y = tot.hill.0)
plt

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C - Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave("outputs/figures/plants_species_diversity/model_restfact_tothill0_restmeth_emm_se.jpg",
       dpi = 300, width = 16, height = 14, units = "cm")

