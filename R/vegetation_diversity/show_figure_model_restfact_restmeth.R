#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Question 2: Restoration factors
# Plot figure: Model restfact_response - Restoration method
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(emmeans) # calculate estimated marginal means and post-hoc Tukey
library(patchwork)
library(glmmTMB)




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
    # panel.grid.major =  element_line(linetype = 1),
    # panel.grid.major.x = element_blank(),
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
# Plot 1 - Total Species Richness #########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill0.Rdata"))



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
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))


# sample size
sample_size = data_model_tot %>%
  group_by(rest.meth) %>%
  summarize(num = n())


## plot with EMMs and SE
(p1 <- ggplot(data_model_tothill0, aes(x = rest.meth, y = tot.hill.0, fill = rest.meth)) +
  gglayer_theme +
  gglayer_labs +
  geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
  geom_pointrange(
    data = emm.df,
    aes(y = response, ymin = response - SE, ymax = response + SE),
    size = 0.5, linewidth = 0.8
  ) +
  geom_text(data = cld_results, aes(x = rest.meth, y = 85, label = .group),          # Add compact letters
            vjust = 0, hjust = 0.4, size = rel(4),
            family = "Franklin Gothic Book", color = "black", 
            # fontface = "bold"
  ) +
  labs(
    title = "Total species richness",
    y = "no. of species",
    # caption = "Tukey adjusted pair-wise comparisons"
    # tag = "A"
  ) +
  coord_cartesian(ylim = c(0, 85))
)



rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 2 - Characteristic Species Richness ################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill0.Rdata"))



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
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))



## plot with EMMs and SE
(p2 <- ggplot(data_model_targethill0, aes(x = rest.meth, y = target.hill.0, fill = rest.meth)) +
    gglayer_theme +
    gglayer_labs +
    geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.5, linewidth = 0.8
    ) +
    geom_text(data = cld_results, aes(x = rest.meth, y = 85, label = .group),          # Add compact letters
              vjust = 0, hjust = 0.4, size = rel(4),
              family = "Franklin Gothic Book", color = "black", 
              # fontface = "bold"
    ) +
    labs(
      title = "Characteristic species richness",
      y = "no. of species",
      # caption = "Tukey adjusted pair-wise comparisons"
      # tag = "A"
    ) +
    coord_cartesian(ylim = c(0, 85))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 3 - Total Hill-Shannon #############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_tothill1.Rdata"))


# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_tothill1, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response   SE  df
# cus           9.76 2.02 Inf
# mga          13.46 2.78 Inf
# res          12.76 2.54 Inf
# dih          14.38 2.83 Inf
# 
# Results are averaged over the levels of: land.use.hist 

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))




## plot with EMMs and SE
(p3 <- ggplot(data_model_tothill1, aes(x = rest.meth, y = tot.hill.1, fill = rest.meth)) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = response, ymin = response - SE, ymax = response + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(data = cld_results, aes(x = rest.meth, y = 36, label = .group),          # Add compact letters
             vjust = 0, hjust = 0.4, size = rel(4),
             family = "Franklin Gothic Book", color = "black", 
             # fontface = "bold"
   ) +
   labs(
     title = "Total Hill-Shannon",
     y = "effective no. of species",
     # caption = "Tukey adjusted pair-wise comparisons"
     # tag = "A"
   ) +
  coord_cartesian(ylim = c(0, 37))
)


rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 4 - Characteristic Hill-Shannon ####################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_targethill1.Rdata"))



# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_targethill1, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response   SE  df
# cus           8.63 1.86 Inf
# mga          11.74 2.52 Inf
# res          11.77 2.43 Inf
# dih          12.35 2.52 Inf
# 
# Results are averaged over the levels of: land.use.hist 

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))




## plot with EMMs and SE
(p4 <- ggplot(data_model_targethill1, aes(x = rest.meth, y = target.hill.1, fill = rest.meth)) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = response, ymin = response - SE, ymax = response + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(data = cld_results, aes(x = rest.meth, y = 36, label = .group),          # Add compact letters
             vjust = 0, hjust = 0.4, size = rel(4),
             family = "Franklin Gothic Book", color = "black",
             # fontface = "bold"
   ) +
   labs(
     title = "Characteristic Hill-Shannon",
     y = "effective no. of species",
     # caption = "Tukey adjusted pair-wise comparisons"
     # tag = "A"
   ) +
  coord_cartesian(ylim = c(0, 37))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))






#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 5 - FCSi ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fcsihill0.Rdata"))



# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_fcsihill0, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth emmean    SE  df
# cus         3.26 0.161 113
# mga         3.34 0.160 113
# res         3.47 0.155 113
# dih         3.58 0.153 113
# 
# Results are averaged over the levels of: land.use.hist 

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))




## plot with EMMs and SE
(p5 <- ggplot(data_model_fcsihill0, aes(x = rest.meth, y = fcsi.hill.0, fill = rest.meth,)) +
   gglayer_theme +
   gglayer_labs +
   geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
   geom_pointrange(
     data = emm.df,
     aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE),
     size = 0.5, linewidth = 0.8
   ) +
   geom_text(data = cld_results, aes(x = rest.meth, y = 4.3, label = .group),          # Add compact letters
             vjust = 0, hjust = 0.4, size = rel(4),
             family = "Franklin Gothic Book", color = "black", 
             # fontface = "bold"
   ) +
   labs(
     title = "Fav. Cons. Status Index (FCSi)",
     y = "FCSi",
     # caption = "Tukey adjusted pair-wise comparisons"
     # tag = "A"
   ) +
  coord_cartesian(ylim = c(1.9, 4.3))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot 6 - Forb-Grass Ratio ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## load data

load(file = here("outputs", "models", "vegetation", "model_plants_restfact_fgratio.Rdata"))



# calculate estimated marginal means (EMMs) and standard error (SE) for each group level
emm.rest.meth <- emmeans(restfact_fgratio, ~ rest.meth)
emm.df <- summary(emm.rest.meth, infer = F, type = "response") # type ="response" for back-transformation from log-scale
emm.df
# rest.meth response    SE  df
# cus          0.735 0.239 Inf
# mga          0.625 0.201 Inf
# res          1.185 0.362 Inf
# dih          1.059 0.304 Inf
# 
# Results are averaged over the levels of: land.use.hist 

# compact letter display
cld_results <- multcomp::cld(emm.rest.meth, adjust = "tukey", Letters = letters) %>% 
  mutate(.group = str_trim(.group))




## plot with EMMs and SE
(p6 <- ggplot(data_model_fgratio, aes(x = rest.meth, y = fg.ratio, fill = rest.meth,)) +
    gglayer_theme +
    gglayer_labs +
    geom_violin(trim = T, width = 0.8, linewidth = 0.5) +
    geom_pointrange(
      data = emm.df,
      aes(y = response, ymin = response - SE, ymax = response + SE),
      size = 0.5, linewidth = 0.8
    ) +
    geom_text(data = cld_results, aes(x = rest.meth, y = 5.3, label = .group),          # Add compact letters
              vjust = 0, hjust = 0.4, size = rel(4),
              family = "Franklin Gothic Book", color = "black", 
              # fontface = "bold"
    ) +
    labs(
      title = "Forb-grass ratio",
      y = "cover of forbs / cover of grasses",
      # caption = "Tukey adjusted pair-wise comparisons"
      # tag = "A"
    ) +
    coord_cartesian(ylim = c(0, 5.5))
)

rm(list = setdiff(ls(), c("my_col", "gglayer_labs", "gglayer_theme",
                          "p1", "p2", "p3", "p4", "p5", "p6")))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Joint plot  #############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# (plot_3x1 <- p1 + p2 + p3)

# (plot_3x3 <- (p1 + p2 + p3) / (p4 + plot_spacer() + plot_spacer() )/ (p7 + plot_spacer() + plot_spacer() ))

# p1_mod <- p1 +
#   coord_cartesian(ylim = 0, 150)

# layout <- "
# ABC
# #EF
# GH#
# "

layout_3x2 <- "
A#B
C#D
E#F
"

# insert tags
p1m <- p1 + labs(tag = "(A)")
p2m <- p2 + labs(tag = "(B)")
p3m <- p3 + labs(tag = "(C)")
p4m <- p4 + labs(tag = "(D)")
p5m <- p5 + labs(tag = "(E)")
p6m <- p6 + labs(tag = "(F)")

# OR

# insert tags in title
p1m <- p1 + labs(title = "(A) Total species richness")
p2m <- p2 + labs(title = "(B) Char. species richness")
p3m <- p3 + labs(title = "(C) Total Hill-Shannon")
p4m <- p4 + labs(title = "(D) Char. Hill-Shannon")
p5m <- p5 + labs(title = "(E) Fav. Cons. Status Index (FCSi)")
p6m <- p6 + labs(title = "(F) Forb-grass ratio")

# remove axis titles and texts
p1m <- p1m + theme(axis.text.x = element_blank())
p2m <- p2m + theme(
  axis.text.x = element_blank(),
  axis.title.y = element_blank())
p3m <- p3m + theme(axis.text.x = element_blank())
p4m <- p4m + theme(
  axis.text.x = element_blank(),
  axis.title.y = element_blank())


(plot_3x2 <- p1m + p2m + p3m + p4m + p5m + p6m +
    plot_layout(design = layout_3x2,
                widths = c(1, 0.1, 1)))

# layout_3x2 <- "
# AB
# CD
# EF
# "
# 
# (plot_3x2 <- p1 + p2 + p4 + p5 + p7 + p8 +
#     plot_layout(design = layout_3x2))


# layout_3x3 <- "
# ABC
# DEF
# GHI
# "
# 
# (plot_3x3 <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
#     plot_layout(design = layout_3x3))

# plot_3x3 <- plot_3x3 & theme(
#   plot.margin = margin(t = rel(20),
#                        r = rel(1),
#                        b = rel(1),
#                        l = rel(1)
#   )
# )
# 
# 
# plot_3x3 <- plot_3x3 & theme(
#   axis.text = element_text(size = rel(1.4)),
#   axis.title.y = element_text(size = rel(1.4)),
#   plot.title = element_text(size = rel(1.6)),
#   plot.margin = margin(t = rel(20),
#                        r = rel(5),
#                        b = rel(1),
#                        l = rel(5)
#                        )
# )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Save  ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 3 x 1 grid
# ggsave(
#   plot = plot_3x1,
#   "outputs/figures/plants_species_diversity/model_restfact_restmeth_emm_se_3x1.jpg",
#   dpi = 300, width = 36, height = 10.5, units = "cm")

# 3 x 2 grid
ggsave(
  plot = plot_3x2,
  "outputs/figures/plants_species_diversity/model_restfact_restmeth_emm_se_3x2.svg",
  dpi = 300, width = 24, height = 24, units = "cm"
  )

# alternative:
# 3 x 2 grid
ggsave(
  plot = plot_3x2,
  "outputs/figures/plants_species_diversity/model_restfact_restmeth_emm_se_3x2_alt.svg",
  dpi = 300, width = 24, height = 24, units = "cm"
)


# 3 x 3 grid
# ggsave(
#   plot = plot_3x3,
#   "outputs/figures/plants_species_diversity/model_restfact_restmeth_emm_se_3x3.jpg",
#   dpi = 300, width = 32, height = 24, units = "cm")


