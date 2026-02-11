#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure S3 ####
# Soil moisture effects in reference and restored sites
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke & Alina Twerski
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

### Start ###
rm(list = ls())

### Colours ###
my_col <- c("#66027e", "#fde725" ,"#21918c")
#my_col <- c("#4E79A7", "#F28E2B" ,"#E15759", "#59A14F")


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
# Plot 1 - Total Species Richness #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p1 <- ggplot(
  data_all, aes(x = site.type, y = tot.hill.0, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Total species richness",
    y = "No. of species"
  ) +
  coord_cartesian(ylim = c(0, 85))
)

# rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
#                           "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 2 - Characteristic Species Richness ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p2 <- ggplot(
  data_all, aes(x = site.type, y = target.hill.0, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  
  labs(
    title = "Characteristic species richness",
    y = "No. of species"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Total Hill-Shannon #################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p3 <- ggplot(
  data_all, aes(x = site.type, y = tot.hill.1, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Total Hill-Shannon",
    y = "ENS"
  ) +
  coord_cartesian(ylim = c(0, 41))
)





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 4 - Characteristic Hill-Shannon ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p4 <- ggplot(
  data_all, aes(x = site.type, y = target.hill.1, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Characteristic Hill-Shannon",
    y = "ENS"
  ) +
  coord_cartesian(ylim = c(0, 41))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 5 - FCSi ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p5 <- ggplot(
  data_all, aes(x = site.type, y = fcsi.hill.0, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Fav. Cons. Status Index (FCSi)",
    y = "FCSi"
  ) +
  coord_cartesian(ylim = c(1.7, 4.3))
)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 6 - Forb-Graminoid Ratio ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p6 <- ggplot(
  data_all, aes(x = site.type, y = fg.ratio, fill = site.type)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Forb-graminoid ratio",
    y = "Cover of forbs / cover of graminoids"
  ) +
  coord_cartesian(ylim = c(0, 5.5))
)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint plot  #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



layout_3x2 <- "
A#B
C#D
E#F
"

# insert tags
p1 <- p1 + labs(tag = "A")
p2 <- p2 + labs(tag = "B")
p3 <- p3 + labs(tag = "C")
p4 <- p4 + labs(tag = "D")
p5 <- p5 + labs(tag = "E")
p6 <- p6 + labs(tag = "F")


(plot_3x2 <- p1 + p2 + p3 + p4 + p5 + p6 +
    plot_layout(design = layout_3x2,
                widths = c(1, 0.15, 1)))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave(
  plot = plot_3x2,
  here("outputs", "figures", "figure_S3_300dpi_24x24cm.tiff"),
  dpi = 300, width = 24, height = 24, units = "cm"
  )

