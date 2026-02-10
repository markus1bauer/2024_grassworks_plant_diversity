#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure S2 ####
# Soil moisture effects in restoration methods
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke
# 2026



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
         hydrology = fct_relevel(hydrology, "moist", "fresh", "dry"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 1 - Total Species Richness #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(p1 <- ggplot(
  data_all, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth)
) +
  gglayer_theme +
  gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
    title = "Total species richness",
    y = "no. of species"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 2 - Characteristic Species Richness ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(p2 <- ggplot(
  data_all, aes(x = rest.meth, y = target.hill.0, fill = rest.meth)
) +
    gglayer_theme +
    gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
      title = "Characteristic species richness",
      y = "no. of species"
    ) +
    coord_cartesian(ylim = c(0, 85))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Total Hill-Shannon #################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(p3 <- ggplot(
  data_all, aes(x = rest.meth, y = tot.hill.1, fill = rest.meth)
) +
    gglayer_theme +
    gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
      title = "Total Hill-Shannon",
      y = "effective no. of species"
    ) +
    coord_cartesian(ylim = c(0, 37))
)






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 4 - Characteristic Hill-Shannon ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(p4 <- ggplot(
  data_all, aes(x = rest.meth, y = target.hill.1, fill = rest.meth)
) +
    gglayer_theme +
    gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
      title = "Characteristic Hill-Shannon",
      y = "effective no. of species"
    ) +
    coord_cartesian(ylim = c(0, 37))
)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 5 - FCSi ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(p5 <- ggplot(
  data_all, aes(x = rest.meth, y = fcsi.hill.0, fill = rest.meth)
) +
    gglayer_theme +
    gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
      title = "Fav. Cons. Status Index (FCSi)",
      y = "FCSi"
    ) +
    coord_cartesian(ylim = c(1.9, 4.3))
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 6 - Forb-Grass Ratio ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(p6 <- ggplot(
  data_all, aes(x = rest.meth, y = fg.ratio, fill = rest.meth)
) +
    gglayer_theme +
    gglayer_labs +
  facet_wrap(~hydrology) +
  geom_boxplot() +
  labs(
      title = "Forb-grass ratio",
      y = "cover of forbs / cover of grasses"
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
p1 <- p1 + labs(tag = "(a)")
p2 <- p2 + labs(tag = "(b)")
p3 <- p3 + labs(tag = "(c)")
p4 <- p4 + labs(tag = "(d)")
p5 <- p5 + labs(tag = "(e)")
p6 <- p6 + labs(tag = "(f)")


(plot_3x2 <- p1 + p2 + p3 + p4 + p5 + p6 +
    plot_layout(design = layout_3x2,
                widths = c(1, 0.1, 1)))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  #######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ggsave(
  plot = plot_3x2,
  here("outputs", "figures", "figure_S4_300dpi_24x24cm.tiff"),
  dpi = 300, width = 24, height = 24, units = "cm"
)
