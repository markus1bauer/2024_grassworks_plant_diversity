#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure 2F ####
# Restoration sites and references
# Forb-graminoid ratio (cover)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Christin Juno Laschke & Markus Bauer
# 2026-02-10



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
rm(list = setdiff(ls(), c(
  "graph_a", "graph_b", "graph_c", "graph_d", "graph_e", "graph_f", "m"
)))

### Theme ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 11, color = "black"),
    axis.text = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.line = element_line(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot ########################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# load data
load(file = here("outputs", "models", "model_reference_forb_grass.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.site.type <- emmeans(restref_fgratio, ~ site.type)
emm.df <- summary(emm.site.type, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.site.type, adjust = "tukey", Letters = letters
) %>% 
  mutate(.group = str_trim(.group))

## plot with EMMs and SE
(graph <- ggplot(
  data_model_fgratio, aes(x = site.type, y = fg.ratio, fill = site.type)
) +
    geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.4, linewidth = 0.8
    ) +
    geom_text(
      data = cld_results, aes(x = site.type, y = 5.3, label = .group),
      vjust = 0, hjust = 0.4, size = rel(4)
    ) +
    labs(
      title = "Forb\u2013graminoid ratio (cover)",
      y = "Ratio"
    ) +
    scale_y_continuous(limits = c(0, 5.5), breaks = seq(0, 10, 1)) +
    scale_x_discrete(
      labels = c("Negative\nreference", "Restored", "Positive\nreference")
    ) +
    scale_fill_manual(values = c("#66027e", "#fde725" ,"#21918c")) +
    theme_mb()
)

# Save ####

ggsave(
  plot = graph,
  here(
    "outputs", "figures", "figure_2_forb-graminoid_ratio_300dpi_8x6cm.tiff"
  ),
  dpi = 300, width = 8, height = 6, units = "cm"
)

graph_f <- graph