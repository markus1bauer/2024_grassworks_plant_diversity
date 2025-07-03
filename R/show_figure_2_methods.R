#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure 2 ####
# Restoration methods
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke
# 2025



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preparation #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(tidyverse)
library(here)
library(emmeans)
library(patchwork)
library(glmmTMB)
library(multcompView)


### Start ###
rm(list = ls())

### Colours ###
my_col <- c("#4E79A7", "#F28E2B" ,"#E15759", "#59A14F")

### Labels ###
gglayer_labs <- list(
  scale_x_discrete(
    labels = c(
      "Cultivar\nseed mixture",
      "Management\nadaptation",
      "Regional\nseed mixture",
      "Direct\nharvesting"),
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
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(),
    panel.grid = element_blank(),
    plot.title = element_text(
      family = "Trebuchet MS",
      size = rel(1.2),
      face = "bold"
    ),
    plot.caption = element_text(size = rel(0.5)),
    plot.tag = element_text(size = rel(1.5), face = "bold")
  ),
  scale_fill_manual(values = my_col),
  guides(y = guide_axis(minor.ticks = TRUE))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 1 - Total Species Richness #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_total_hill0_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_tothill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

# sample size
sample_size = data_model_tothill0 %>%
  group_by(rest.meth) %>%
  summarize(num = n())

## plot with EMMs and SE
(p1 <- ggplot(
  data_model_tothill0, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth)
  ) +
  gglayer_theme +
  gglayer_labs +
  geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 0.5, linewidth = 0.8
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 85, label = .group),
            vjust = 0, hjust = 0.4, size = rel(4),
            family = "Franklin Gothic Book", color = "black"
  ) +
  labs(
    title = "Total species richness",
    y = "no. of species"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 2 - Characteristic Species Richness ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_target_hill0_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_targethill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(p2 <- ggplot(
  data_model_targethill0,
  aes(x = rest.meth, y = target.hill.0, fill = rest.meth)
  ) +
    gglayer_theme +
    gglayer_labs +
    geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.5, linewidth = 0.8
    ) +
    geom_text(
      data = cld_results, aes(x = rest.meth, y = 85, label = .group),
      vjust = 0, hjust = 0.4, size = rel(4),
      family = "Franklin Gothic Book", color = "black"
    ) +
    labs(
      title = "Characteristic species richness",
      y = "no. of species"
    ) +
    coord_cartesian(ylim = c(0, 85))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Total Hill-Shannon #################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_total_hill1_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_tothill1, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(p3 <- ggplot(
  data_model_tothill1, aes(x = rest.meth, y = tot.hill.1, fill = rest.meth)
  ) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = response, ymin = response - SE, ymax = response + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(
     data = cld_results, aes(x = rest.meth, y = 36, label = .group),
     vjust = 0, hjust = 0.4, size = rel(4),
     family = "Franklin Gothic Book", color = "black"
   ) +
   labs(
     title = "Total Hill-Shannon",
     y = "effective no. of species"
   ) +
  coord_cartesian(ylim = c(0, 37))
)


rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 4 - Characteristic Hill-Shannon ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_target_hill1_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_targethill1, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(p4 <- ggplot(
  data_model_targethill1,
  aes(x = rest.meth, y = target.hill.1, fill = rest.meth)
  ) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = response, ymin = response - SE, ymax = response + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(
     data = cld_results, aes(x = rest.meth, y = 36, label = .group),
     vjust = 0, hjust = 0.4, size = rel(4),
     family = "Franklin Gothic Book", color = "black"
   ) +
   labs(
     title = "Characteristic Hill-Shannon",
     y = "effective no. of species"
   ) +
  coord_cartesian(ylim = c(0, 37))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 5 - FCSi ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_fcsi_hill0_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_fcsihill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(p5 <- ggplot(
  data_model_fcsihill0, aes(x = rest.meth, y = fcsi.hill.0, fill = rest.meth)
  ) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(
     data = cld_results, aes(x = rest.meth, y = 4.3, label = .group),
     vjust = 0, hjust = 0.4, size = rel(4),
     family = "Franklin Gothic Book", color = "black"
   ) +
   labs(
     title = "Fav. Cons. Status Index (FCSi)",
     y = "FCSi"
   ) +
  coord_cartesian(ylim = c(1.9, 4.3))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 6 - Forb-Grass Ratio ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data
load(file = here("outputs", "models", "model_methods_forb_grass_full.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.rest.meth <- emmeans(restfact_fgratio, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.rest.meth, adjust = "tukey", Letters = letters
  ) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(p6 <- ggplot(
  data_model_fgratio, aes(x = rest.meth, y = fg.ratio, fill = rest.meth)
  ) +
    gglayer_theme +
    gglayer_labs +
    geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.5, linewidth = 0.8
    ) +
    geom_text(data = cld_results, aes(x = rest.meth, y = 5.3, label = .group),
              vjust = 0, hjust = 0.4, size = rel(4),
              family = "Franklin Gothic Book", color = "black"
    ) +
    labs(
      title = "Forb-grass ratio",
      y = "cover of forbs / cover of grasses"
    ) +
    coord_cartesian(ylim = c(0, 5.5))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint plot  #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



layout_3x2 <- "
A#B
C#D
E#F
"

# insert tags
p1 <- p1 + labs(tag = "(A)")
p2 <- p2 + labs(tag = "(B)")
p3 <- p3 + labs(tag = "(C)")
p4 <- p4 + labs(tag = "(D)")
p5 <- p5 + labs(tag = "(E)")
p6 <- p6 + labs(tag = "(F)")


(plot_3x2 <- p1 + p2 + p3 + p4 + p5 + p6 +
    plot_layout(design = layout_3x2,
                widths = c(1, 0.1, 1)))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ggsave(
#   plot = plot_3x2,
#   here("outputs", "figures", "figure_2_300dpi_24x24cm.tiff"),
#   dpi = 300, width = 24, height = 24, units = "cm"
# )
