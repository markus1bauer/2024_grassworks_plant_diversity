#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure S5 ####
# Forb, Graminoids and Legumes cover in reference and restored sites
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
my_col <- c("#66027e", "#fde725" ,"#21918c")


### Labels ###
gglayer_labs <- list(
  scale_x_discrete(
    labels = c(
      "neg",
      "rest",
      "pos"),
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
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"),
         hydrology = fct_relevel(hydrology, "moist", "fresh", "dry"))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 1 - Forb cover #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p1 <- ggplot(
  data_all, aes(x = site.type, y = cover.forbs, fill = site.type)
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
# Plot 2 - Legumes cover ######################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p2 <- ggplot(
  data_all, aes(x = site.type, y = cover.legumes, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  geom_boxplot() +
  
  labs(
    title = "Legumes cover",
    y = "Cover [%]"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Graminoids cover #################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p3 <- ggplot(
  data_all, aes(x = site.type, y = cover.grass, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  geom_boxplot() +
  labs(
    title = "Graminoids cover",
    y = "Cover [%]"
  ) +
  coord_cartesian(ylim = c(0, 110))
)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint plot  #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


cover_panel_reference <- plot_grid(p1, p2, p3,
                         nrow = 3, ncol = 1, 
                        labels = c("A", "B", "C"))
plot(cover_panel_reference)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave(
  plot = cover_panel_reference,
  here("outputs", "figures", "figure_S5_300dpi_12x24cm.tiff"),
  dpi = 300, width = 12, height = 24, units = "cm"
)



