#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Plot figure: Species Diversity
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)


### Start ###
rm(list = ls())


## load data -------------------------------------------------------------------

### environment data ###

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20240813.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, rest.meth, region, rest.age, 
    land.use.hist, hydrology, hydrology.cont, site.cwm.abu.oek.f, lui,
    mngm.type, obs.year,
    c.n, c.perc, n.perc, ph.value, tic.perc, toc.perc,
    clay.perc, silt.perc, sand.perc,
  ) %>% 
  distinct() %>% 
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"),
         region = fct_relevel(region, "north", "centre", "south"))
  

# todo: resolve infinite numbers in c.n.B30 and c.n in sites N_RET and N_HOR
# quick and dirty:
sites[sites$id.site == "N_RET", "c.n"] <- NA
sites[sites$id.site == "N_HOR", "c.n"] <- NA

# choose one restoration method for "dih&seed"
# --> dih because it is most likely that dih has the greatest influence
# --> Line said for M_SKF: ReS
sites %>% 
  filter(rest.meth == "dih&seed")
sites[sites$id.site == "S_FHZ", "rest.meth"] <- "dih"
sites[sites$id.site == "S_GIG", "rest.meth"] <- "dih"
sites[sites$id.site == "M_SKF", "rest.meth"] <- "res"
sites[sites$id.site == "S_WTZ", "rest.meth"] <- "dih"

sites$rest.meth <- droplevels(as.factor(sites$rest.meth))
str(sites)

sites <- sites %>% 
  # add a variable with restoration method and site type
  mutate(rest.meth.type = if_else(
    site.type == "restored", rest.meth, site.type))
# sites <- read_csv(
#   here("data", "processed", "sites_processed_environment_nms.csv"),
#   col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
#     .default = "?"
#   )
# ) %>%
#   # calculate mean of CWM Ellenberg F per site
#   group_by(id.site) %>% 
#   mutate(site_cwm_abu_oek_f = mean(cwm_abu_oek_f),
#          site_cwm_pres_oek_f = mean(cwm_pres_oek_f)) %>%
#   ungroup() %>%
#   mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"),
#          region = fct_relevel(region, "north", "centre", "south")) %>% 
#   # add a variable with restoration method and site type
#   mutate(
#     rest.meth.type = if_else(
#       site.type == "restored", rest.meth, site.type),
#     # add a variable with management type
#     mngm.type = if_else(
#       mngm.mow == 1 & mngm.graz == 1, "both", if_else(
#         mngm.mow == 1, "mowing", if_else(
#           mngm.graz == 1, "grazing", if_else(
#             mngm.no == 1, "none", NA
#           )
#         )
#       ))) %>%  
#   select(
#     id.site, site.type, rest.meth, rest.meth.type, region, rest.age,
#     land.use.hist, hydrology, site_cwm_abu_oek_f, site_cwm_pres_oek_f, lui,
#     mngm.type
#   ) %>% 
#   distinct() %>% 
#   # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
#   filter(!is.na(id.site))


### Diversity data ###

diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20240814.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))

# add to environment data
sites <- sites %>% 
  left_join(diversity, by = "id.site")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Figures ###################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 1 - Species Richness -----------------------------------------------------------

plot_data <- sites %>% 
  pivot_longer(
    cols = starts_with(c("tot", "target")),
    names_to = "system",
    values_to = "richness"
  ) 

plot_target_data <- plot_data %>% 
  filter(system %in% c(
    "tot.hill.0",
    "target.hill.0"
  ))


get_box_stats <- function(y, upper_limit = stats_height * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y), "\n"
    )
  ))
}

stats_height <- max(plot_target_data$richness, na.rm = T)

palette("Tableau 10")
# palette("Classic Tableau")
# palette("R4")
# palette("Dark2")
# palette("Set1")

### by site type ####
plot_target_data %>% 
  ggplot(aes(x = site.type, y = richness, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  # labs(x = "order q", y = "Effective number of target species")+
  facet_grid(cols=vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette())

  # scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-sitetype.jpg")

#### and region ####
plot_target_data %>% 
  # filter(system == "richness_ffh_calth_btt") %>% 
  ggplot(aes(x = site.type, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette())
  # scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-sitetype_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

### by restoration method ####
plot_target_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(
    # x = fct_reorder(rest.meth, richness, median),
    x = rest.meth,
    y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system),) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette()[c(1:3,5)])

  # scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-rest_meth.jpg")

#### and region ####
plot_target_data %>% 
  # filter(system == "richness_ffh_calth_btt") %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = rest.meth, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-restmeth_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


#### restored (by method) and reference sites ####
plot_target_data %>% 
  filter(rest.meth.type != "dih&seed") %>%
  ggplot(aes(x = fct_relevel(rest.meth.type, "negative", "cus", "mga", "res", "dih", "positive"),
             y = richness, fill = rest.meth.type)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette()[c(1:3,5,10,10)])
# scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-rest_meth+site_type.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

palette("Classic Tableau")

### by region ####
plot_target_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = region, y = richness, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-region.jpg")


### by age of restoration ####
plot_target_data %>% 
  filter(rest.age != "NA") %>% 
  ggplot(aes(x = rest.age, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/richness-age.jpg")

#### and region ####
plot_target_data %>% 
  filter(rest.age != "NA") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = rest.age, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_manual(values = palette())

  # scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         # name = "Region",
  # )
ggsave("outputs/figures/plants_species_diversity/richness-age_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


#### and region ####


plot_target_data %>% 
  filter(rest.age != "NA") %>% 
  filter(system == "tot.hill.0") %>%
  ggplot(aes(x = rest.age, y = richness, colour = region)) +
  geom_point() +
  # facet_grid(cols=vars(system)) + 
  # geom_smooth(formula = , se = F) +
  geom_line(aes(rest.age, pred_hill0_glmm), data = newdat1) +
  labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_manual(values = palette()) +
  theme_bw()

# scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
# name = "Region",
# )
ggsave("outputs/figures/plants_species_diversity/richness-age_region_SERE.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")



#### and method ####

plot_target_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(system == "tot.hill.0") %>%
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "res", "mga", "dih")) %>% 
  ggplot(aes(x = rest.age, y = richness, colour = rest.meth)) +
  geom_point(size = 4) +
  # facet_grid(cols=vars(system)) + 
  # geom_smooth(method = "glm", se = F) +
  labs(x = "Restoration Age", y = "Species Richness (16 m²)") +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(40, 40, 40, 40)) +
  scale_colour_manual(values = palette()[c(1,2,4,3)],
                      name = "Restoration Method",
                      labels = c("Cultivar Seed Mixture",
                                 "Regional Seed Mixture",
                                 "Management Adaptation",
                                 "Direct Harvesting")
  )
ggsave("outputs/figures/plants_species_diversity/richness-age_restmeth_SERE.jpg",
       dpi = 300, width = 22, height = 16, units = "cm")

### by land use history ####
plot_target_data %>% 
  filter(land.use.hist != "NA") %>% 
  ggplot(aes(x = land.use.hist, y = richness, fill = land.use.hist)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-land_use_history.jpg")

#### and method ####
plot_target_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(land.use.hist)) %>% 
  ggplot(aes(x = land.use.hist, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-land_use_hist_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


#  SERE figure
data <- plot_target_data %>%
  filter(system == "tot.hill.0", !is.na(land.use.hist), !is.na(rest.meth)) %>% 
  mutate(rest.meth = fct_relevel(rest.meth.type, "cus", "res", "mga", "dih"))


ggboxplot(data, x = "rest.meth", y = "richness",
          fill = "land.use.hist",
          legend = "right",
) +
  labs(
    y = "Species Richness (16 m²)",
    # title = "Previous land use",
    caption = "n = 121",
  ) +
  xlab(NULL) +
  scale_x_discrete(labels = c("Cultivar Seed Mixture",
                              "Regional Seed Mixture",
                              "Management Adaptation",
                              "Direct Harvesting")) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(40, 40, 40, 40)) +
  rotate_x_text(45) +
  # geom_text(
  #   data = dt,
  #   aes(x = rest.meth.type, y = 80, label = cld), vjust = -0.7) +
  scale_fill_manual(
    values = c("lightgrey", "darkgrey"),
    name = "Previous land use",
    # labels = c("Negative Reference",
    #            "Cultivar Seed Mixture",
    #            "Regional Seed Mixture",
    #            "Management Adaptation",
    #            "Direct Harvesting",
    #            "Positive Reference")
    )
ggsave("outputs/figures/plants_species_diversity/richness-land_use_hist_restmeth_SERE.jpg",
       dpi = 300, width = 22, height = 20, units = "cm")





### by hydrology ####
#### factorial ####
plot_target_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = richness, fill = hydrology)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-hydrology.jpg")

##### and region ####
plot_target_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-hydrology_region.jpg")

##### and method ####
plot_target_data %>% 
  filter(!is.na(rest.meth)) %>%
  ggplot(aes(x = hydrology, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-hydrology_method.jpg")

#### CWM Ellenberg F abundance ####
plot_target_data %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/richness-CWM_abu_oek_f.jpg")

##### and region ####
plot_target_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = hydrology, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/richness-CWM_abu_oek_f_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

##### and method ####
plot_target_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = richness, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/richness-CWM_abu_oek_f_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

# #### CWM Ellenberg F presence ###
# plot_target_data %>% 
#   ggplot(aes(x = site_cwm_pres_oek_f, y = richness, colour = system)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   # labs(x = "Restoration Age", y = "Target Species Richness") +
#   scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
#                          name = "system",
#   )
# ggsave("outputs/figures/plants_species_diversity/richness-CWM_pres_oek_f.jpg")


### by LUI ####
plot_target_data %>% 
  ggplot(aes(x = lui, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/richness-LUI.jpg")

#### and region ####
plot_target_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/richness-LUI_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_target_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = richness, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/richness-LUI_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


### by management ####
plot_target_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = mngm.type)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-mngm.jpg")

#### and region ####
plot_target_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = region)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-mngm_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_target_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/richness-mngm_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


## 2 - FCSi ------------------------------------------------------------------------


plot_data <- sites %>% 
  pivot_longer(
    cols = starts_with("fcsi"),
    names_to = "system",
    values_to = "fcsi"
  ) %>% 
  filter(system %in% c(
    "fcsi.hill.0"
    # , "fcsi.hill.1"
    # , "fcsi.hill.2"
  ))

# check NA
plot_data %>% 
  filter(is.na(fcsi))
# no NA

get_box_stats <- function(y, upper_limit = stats_height * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y), "\n",
      "FCSi =", round(median(y), 2), "\n"
    )
  ))
}

stats_height <- max(plot_data$fcsi, na.rm = T)

### by site type ####
plot_data %>% 
  ggplot(aes(x = site.type, y = fcsi, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  # labs(x = "order q", y = "Effective number of target species")+
  facet_grid(cols=vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-site_type.jpg")

#### and region ####
plot_data %>% 
  # filter(system == "richness_ffh_calth_btt") %>% 
  ggplot(aes(x = site.type, y = fcsi, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-sitetype_region.jpg")

### by restoration method ####
plot_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "mga", "res", "dih")) %>% 
  ggplot(aes(x = rest.meth,
             y = fcsi, fill = rest.meth)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-rest_meth.jpg")

#### and region ####
plot_data %>% 
  # filter(system == "fcsi_ffh_calth_btt") %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = rest.meth, y = fcsi, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-restmeth_region.jpg")


#### restored (by method) and reference sites ####
plot_data %>% 
  filter(rest.meth.type != "dih&seed") %>%
  ggplot(aes(x = fct_relevel(rest.meth.type, "negative", "cus", "mga", "res", "dih", "positive"),
             y = fcsi, fill = rest.meth.type)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system)) +
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-rest_meth+site_type.jpg"
       # , dpi = 300, width = 30, height = 12, units = "cm"
)


### by region ####
plot_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = region, y = fcsi, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-region.jpg")


### by age of restoration ####
plot_data %>% 
  filter(rest.age != "NA") %>% 
  ggplot(aes(x = rest.age, y = fcsi, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Restoration Age", y = "Index of Favourable Conservation Status (FCSi)") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-age.jpg")

#### and region ####
plot_data %>% 
  filter(rest.age != "NA") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = rest.age, y = fcsi, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-age_region.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm"
)

#### and method ####
plot_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = rest.age, y = fcsi, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-age_restmeth.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm"
)


### by land use history ####
plot_data %>% 
  filter(land.use.hist != "NA") %>% 
  ggplot(aes(x = land.use.hist, y = fcsi, fill = land.use.hist)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-land_use_history.jpg")

#### and method ####
plot_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(land.use.hist)) %>% 
  ggplot(aes(x = land.use.hist, y = fcsi, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-land_use_hist_restmeth.jpg")

### by hydrology ####
#### factorial ####
plot_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = fcsi, fill = hydrology)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-hydrology.jpg")

##### and region ####
plot_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = fcsi, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-hydrology_region.jpg")

##### and method ####
plot_data %>% 
  filter(!is.na(rest.meth)) %>%
  ggplot(aes(x = hydrology, y = fcsi, fill = rest.meth)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-hydrology_method.jpg")

#### CWM Ellenberg F abundance ####
plot_data %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = fcsi, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-CWM_abu_oek_f.jpg")

##### and region ####
plot_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = fcsi, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-CWM_abu_oek_f_region.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm"
)

##### and method ####
plot_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = fcsi, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/fsci-CWM_abu_oek_f_restmeth.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm
)

# #### CWM Ellenberg F presence ###
# plot_data %>% 
#   ggplot(aes(x = site_cwm_pres_oek_f, y = fcsi, colour = system)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   # labs(x = "Restoration Age", y = "Target Species Richness") +
#   scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
#                          name = "system",
#   )
# ggsave("outputs/figures/plants_species_diversity/fcsi-CWM_pres_oek_f.jpg")

### by LUI ####
plot_data %>% 
  ggplot(aes(x = lui, y = fcsi, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-LUI.jpg")

#### and region ####
plot_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = fcsi, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-LUI_region.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm"
)

#### and method ####
plot_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = fcsi, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/fcsi-LUI_restmeth.jpg"
       # , dpi = 300, width = 25, height = 12, units = "cm"
)


### by management ####
plot_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = fcsi, fill = mngm.type)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-mngm.jpg")

#### and region ####
plot_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = fcsi, fill = region)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-mngm_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = fcsi, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/fcsi-mngm_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


## 3 - Hill numbers -----------------------------------------------------------


plot_hill_data <- sites %>% 
  pivot_longer(
    cols = starts_with("tot"),
    names_to = "system",
    values_to = "richness"
  ) 



get_box_stats <- function(y, upper_limit = stats_height * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y), "\n"
    )
  ))
}

stats_height <- max(plot_hill_data$richness, na.rm = T)

# set colour palette
palette("Tableau 10")


### by site type ####
plot_hill_data %>% 
  ggplot(aes(x = site.type, y = richness, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  # labs(x = "order q", y = "Effective number of target species")+
  facet_grid(cols=vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette())

# scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-sitetype.jpg")


#### and region ####
plot_hill_data %>% 
  # filter(system == "richness_ffh_calth_btt") %>% 
  ggplot(aes(x = site.type, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette())
# scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-sitetype_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

### by restoration method ####
plot_hill_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(
    # x = fct_reorder(rest.meth, richness, median),
    x = rest.meth,
    y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system),) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette()[c(1:3,5)])

# scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-rest_meth.jpg")

#### and region ####
plot_hill_data %>% 
  # filter(system == "richness_ffh_calth_btt") %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = rest.meth, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-restmeth_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


#### restored (by method) and reference sites ####
plot_hill_data %>% 
  filter(rest.meth.type != "dih&seed") %>%
  ggplot(aes(x = fct_relevel(rest.meth.type, "negative", "cus", "mga", "res", "dih", "positive"),
             y = richness, fill = rest.meth.type)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols = vars(system)) +
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9, size = 2
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_manual(values = palette()[c(1:3,5,10,10)])
# scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-rest_meth+site_type.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")



### by region ####
plot_hill_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = region, y = richness, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-region.jpg")


### by age of restoration ####
plot_hill_data %>% 
  filter(rest.age != "NA") %>% 
  ggplot(aes(x = rest.age, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-age.jpg")

#### and region ####
plot_hill_data %>% 
  filter(rest.age != "NA") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = rest.age, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_manual(values = palette())

# scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
# name = "Region",
# )
ggsave("outputs/figures/plants_species_diversity/hill012-age_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_hill_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = rest.age, y = richness, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-age_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

### by land use history ####
plot_hill_data %>% 
  filter(land.use.hist != "NA") %>% 
  ggplot(aes(x = land.use.hist, y = richness, fill = land.use.hist)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-land_use_history.jpg")

#### and method ####
plot_hill_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(land.use.hist)) %>% 
  ggplot(aes(x = land.use.hist, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-land_use_hist_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


### by hydrology ####
#### factorial ####
plot_hill_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = richness, fill = hydrology)) +
  geom_boxplot(show.legend = FALSE) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-hydrology.jpg")

##### and region ####
plot_hill_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = hydrology, y = richness, fill = region)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-hydrology_region.jpg")

##### and method ####
plot_hill_data %>% 
  filter(!is.na(rest.meth)) %>%
  ggplot(aes(x = hydrology, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = T) +
  facet_grid(cols=vars(system)) + 
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-hydrology_method.jpg")

#### CWM Ellenberg F abundance ####
plot_hill_data %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-CWM_abu_oek_f.jpg")

##### and region ####
plot_hill_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = hydrology, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-CWM_abu_oek_f_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

##### and method ####
plot_hill_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = site_cwm_abu_oek_f, y = richness, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  # labs(x = "Restoration Age", y = "Target Species Richness") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-CWM_abu_oek_f_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

# #### CWM Ellenberg F presence ###
# plot_hill_data %>% 
#   ggplot(aes(x = site_cwm_pres_oek_f, y = richness, colour = system)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   # labs(x = "Restoration Age", y = "Target Species Richness") +
#   scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
#                          name = "system",
#   )
# ggsave("outputs/figures/plants_species_diversity/hill012-CWM_pres_oek_f.jpg")


### by LUI ####
plot_hill_data %>% 
  ggplot(aes(x = lui, y = richness, colour = system)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "system",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-LUI.jpg")

#### and region ####
plot_hill_data %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = richness, colour = region)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-LUI_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_hill_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  # filter(system == "richness_ffh") %>% 
  ggplot(aes(x = lui, y = richness, colour = rest.meth)) +
  geom_point() +
  facet_grid(cols=vars(system)) + 
  geom_smooth(method = "lm", se = F) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration Method",
  )
ggsave("outputs/figures/plants_species_diversity/hill012-LUI_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")


### by management ####
plot_hill_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = mngm.type)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-mngm.jpg")

#### and region ####
plot_hill_data %>%
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = region)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-mngm_region.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")

#### and method ####
plot_hill_data %>%
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(!is.na(mngm.type)) %>% 
  ggplot(aes(x = mngm.type, y = richness, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  facet_grid(cols=vars(system)) + 
  stat_summary(
    fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/plants_species_diversity/hill012-mngm_restmeth.jpg",
       dpi = 300, width = 25, height = 12, units = "cm")
