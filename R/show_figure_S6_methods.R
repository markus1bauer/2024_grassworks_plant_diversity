#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure S6 ####
# Forb, Graminoids and Legumes cover n restoration methods
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Alina Twerski & Christin Juno Laschke
# 2026-02-11


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
library(cowplot)

### Start ###
rm(list = ls())

### Colours ###
my_col <-  c("#781c6d", "#bc3754" ,"#ed6925", "#fbb61a")


### Labels ###
gglayer_labs <- list(
  scale_x_discrete(
    labels = c(
      "cus",
      "mga",
      "res",
      "dih"),
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
    plot.tag = element_text(size = rel(1.5), face = "bold"),
    strip.text = element_text(size=10)
  ),
  scale_fill_manual(values = my_col),
  guides(y = guide_axis(minor.ticks = TRUE))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# load data #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data_all <- read_csv(
  here("data", "processed", "data_processed.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  filter(!is.na(rest.meth)) %>% 
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih"),
         region = fct_relevel(region, "north", "centre", "south"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 1 - Forb cover #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p1 <- ggplot(
  data_all, aes(x = rest.meth, y = cover.forbs, fill = rest.meth)
) +
  gglayer_theme +
  gglayer_labs +
  geom_boxplot() +
  labs(
    title = "Forb cover",
    y = "Cover [%]"
  ) +
  coord_cartesian(ylim = c(0, 85))
)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 2 - Legume cover ######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p2 <- ggplot(
  data_all, aes(x = rest.meth, y = cover.legumes, fill = rest.meth)
) +
  gglayer_theme +
  gglayer_labs +
  geom_boxplot() +
  
  labs(
    title = "Legume cover",
    y = "Cover [%]"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Graminoid cover #################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p3 <- ggplot(
  data_all, aes(x = rest.meth, y = cover.grass, fill = rest.meth)
) +
  gglayer_theme +
  gglayer_labs +
  geom_boxplot() +
  labs(
    title = "Graminoid cover",
    y = "Cover [%]"
  ) +
  coord_cartesian(ylim = c(0, 110))
)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint plot  #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cover_panel_method <- plot_grid(p1, p2, p3,
                         nrow = 3, ncol = 1, 
                         labels = c("A", "B", "C"))
plot(cover_panel_method)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave(
  plot = cover_panel_method,
  here("outputs", "figures", "figure_S6_300dpi_12x24cm.tiff"),
  dpi = 300, width = 12, height = 24, units = "cm"
)



