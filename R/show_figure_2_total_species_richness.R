#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Show figure 2A ####
# Restoration sites and references
# Total species richness
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



## load data
load(file = here("outputs", "models", "model_reference_total_hill0.Rdata"))

# calculate estimated marginal means (EMMs) and standard error (SE)
emm.site.type <- emmeans(restref_tothill0, ~ site.type)
emm.df <- summary(emm.site.type, infer = FALSE, type = "response")
emm.df

# compact letter display
cld_results <- multcomp::cld(
  emm.site.type, adjust = "tukey", Letters = letters
) %>% 
  mutate(.group = str_trim(.group))

# sample size
sample_size = data_model_tothill0 %>%
  group_by(site.type) %>%
  summarize(num = n())

## plot with EMMs and SE
(graph <- ggplot(
  data_model_tothill0, aes(x = site.type, y = tot.hill.0, fill = site.type)
) +
    geom_violin(trim = TRUE, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.4, linewidth = 0.8
    ) +
    geom_text(
      data = cld_results, aes(x = site.type, y = 84, label = .group),
      vjust = 0, hjust = 0.4, size = rel(4)
    ) +
    labs(
      title = "Total species richness",
      y = "Number of species"
    ) +
    scale_x_discrete(
      labels = c("Negative\nreference", "Restored", "Positive\nreference")
    ) +
    scale_fill_manual(values = c("#66027e", "#fde725" ,"#21918c")) +
    scale_y_continuous(limits = c(0, 88), breaks = seq(0, 100, 10)) +
    theme_mb()
)

# Save ####

ggsave(
  plot = graph,
  here("outputs", "figures", "figure_2_total_richness_300dpi_8x6cm.tiff"),
  dpi = 300, width = 8, height = 6, units = "cm"
)

graph_a <- graph +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )