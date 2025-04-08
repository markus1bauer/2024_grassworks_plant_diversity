#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation diversity analysis
# Descriptive figures
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(hilldiv)
library(hillR)
library(vegan)
library(nlme)
library(lme4)
library(lmerTest)
library(DHARMa)
library(gridExtra)


### Start ###
rm(list = ls())

## load data -------------------------------------------------------------------
### site environment data ####

sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) %>%
  # calculate mean of CWM Ellenberg F per site
  group_by(id.site) %>% 
  mutate(site_cwm_abu_oek_f = mean(cwm_abu_oek_f),
         site_cwm_pres_oek_f = mean(cwm_pres_oek_f)) %>%
  ungroup() %>%
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"),
         region = fct_relevel(region, "north", "centre", "south")) %>% 
  # add a variable with restoration method and site type
  mutate(
    rest.meth.type = if_else(
      site.type == "restored", rest.meth, site.type),
    # add a variable with management type
    mngm.type = if_else(
      mngm.mow == 1 & mngm.graz == 1, "both", if_else(
        mngm.mow == 1, "mowing", if_else(
          mngm.graz == 1, "grazing", if_else(
            mngm.no == 1, "none", NA
          )
        )
      ))) %>% 
  select(
    id.site, site.type, rest.meth, rest.meth.type, region, start.rest,
    land.use.hist, hydrology, site_cwm_abu_oek_f, site_cwm_pres_oek_f, lui,
    mngm.type
  ) %>%
  distinct() %>% 
  # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
  filter(!is.na(id.site))



### species abundance data ####
# to do: change input data, use not-preliminary data
abundances_site <- read_csv(
  here(
    "data", "processed", "data_processed_plants_nms_site_abundances.csv"
  ),
  col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
) %>% 
  # select(id.site, subtransect, name.tnrs, cover.median) %>% 
  rename(name.plant = name.tnrs)


# # to do: change input data, use processed data
# abundances_plot <- read_csv(
#   here(
#     "data", "raw", "species_vegetation",
#     "data_raw_plants_nms_abundances_20240122.csv"
#   ),
#   col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
# ) %>%
#   # remove duplicated rows
#   # to do: check in not-preliminary data if there are still duplicates
#   distinct() %>%
#   # remove mistakes: double entries with different values
#   # to do: check input before analysis with not-preliminary data
#   distinct(id.site, subtransect, name.plant, .keep_all = TRUE) %>% 
#   mutate(id.plot = str_c(id.site, subtransect, sep = "_")) %>%
#   arrange(name.plant) %>%
#   ### Check that each species occurs at least one time ###
#   group_by(name.plant) %>%
#   mutate(
#     total = sum(cover.median, na.rm = T),
#     presence = if_else(total > 0, 1, 0)
#   ) %>%
#   # filter only species which occur at least one time
#   # but keep all species in subtransect "T"
#   filter(!(presence == 0 & subtransect != "T")) %>% 
#   ungroup() %>%
#   select(id.site, subtransect, name.plant, cover.median)





# transform input data ----------------------------------------------------




# ## transform into site_species data frame for hill_div package:
# # converting long into wide data (columns = sites, rows = species)
# veg_spp_site_mat <-
#   pivot_wider(
#     cover.mean,
#     id_cols = plant.name,
#     names_from = site.ID,
#     values_from = cover.median.site,
#     values_fill = 0
#   )
# veg_spp_site_mat <- as.data.frame(veg_spp_site_mat)
# rownames(veg_spp_site_mat) <- veg_spp_site_mat[, 1]  ## set rownames by using first column
# veg_spp_site_mat <- veg_spp_site_mat[, -1]           ## remove the first variable


# Descriptive analysis ----------------------------------------------------


species_abundances_plot <- abundances_plot %>% 
  # exclude plants identified on genus- or family-level
  filter(rank %in% c("species", "subspecies", "variety")) %>%
  # exclude cf, keep NA (TO DO: check in not-preliminary data if NAs are still there and what to do)
  filter(cf == "FALSCH" | is.na(cf))

## species list & total species number
species_list <- species_abundances_plot %>% 
  distinct(name.plant)
# 680 species in total found in all sites

## number of species in subtransects A1-A4 (no random walks)
species_abundances_plot %>% 
  filter(subtransect != "T") %>% 
  distinct(name.plant)
# 565 species found in vegetation plots (subtransects A1-A4)

## number of species in regions
species_abundances_plot %>%
  distinct(name.plant, region) %>% 
  group_by(region) %>% 
  count()
# region     n
# <chr>  <int>
# 1 M        489
# 2 N        291
# 3 S        407

richness <- species_abundances_plot %>% 
  distinct(name.plant, id.site) %>% 
  group_by(id.site) %>% 
  count()

  write_csv(
    richness,
    here(
      "outputs", "table_plants_richness_A1-A4_T.csv"
    )
  )




# Hill numbers -------------------------------------------------------------



## Plots ####
# transform into long format
plot_div_data <- sites %>% 
  pivot_longer(
    cols = starts_with("hill"),
    names_to = "q.order",
    values_to = "eff.no"
  )


### by site type ####
plot_div_data %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(site.type)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
# ggsave("outputs/figures/hill_012-site_type.jpg")

plot_div_data %>% 
  ggplot(aes(x = site.type, y = eff.no, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "site type", y = "Effective number of species")+
  # scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(q.order)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-site_type_2.jpg")


# by site type and region
plot_div_data %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = region)) +
  geom_boxplot(show.legend = T) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(site.type)) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-site_type+region.jpg")


## hill 0
type_aov <- aov(hill.0 ~ site.type,
                data = sites
)
summary(type_aov)
# F(2,183) = 26.77; p < 0.001
shapiro.test(log1p(type_aov$residuals))
bartlett.test(hill.0 ~ site.type, data = sites)
TukeyHSD(type_aov)
# $site.type
# diff      lwr      upr   p adj
# restored-negative 17.482955 11.30382 23.66209 0.00000
# positive-negative 21.331439 13.61931 29.04357 0.00000
# positive-restored  3.848485 -2.25616  9.95313 0.29825

## hill 1
type_aov <- aov(hill.1 ~ site.type,
                data = sites
)
summary(type_aov)
# F(2,183) = 22; p < 0.001
shapiro.test(log1p(type_aov$residuals))
bartlett.test(hill.0 ~ site.type, data = sites)
TukeyHSD(type_aov)
# $site.type
# diff        lwr      upr     p adj
# restored-negative  8.012384  4.6651488 11.35962 0.0000002
# positive-negative 11.075100  6.8974424 15.25276 0.0000000
# positive-restored  3.062717 -0.2441665  6.36960 0.0758369

## hill 2
type_aov <- aov(hill.2 ~ site.type,
                data = sites
)
summary(type_aov)
# F(2,183) = 16.97; p < 0.001
shapiro.test(log1p(type_aov$residuals))
bartlett.test(hill.0 ~ site.type, data = sites)
TukeyHSD(type_aov)
# $site.type
# diff        lwr       upr     p adj
# restored-negative 5.256095  2.7330715  7.779118 0.0000057
# positive-negative 7.368573  4.2196074 10.517539 0.0000003
# positive-restored 2.112478 -0.3801294  4.605086 0.1144536

### by restoration method ####
plot_div_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = rest.meth)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols = vars(
    fct_reorder(rest.meth, eff.no, median)),
    # labeller = labeller(rest.meth = rest.meth.labels)
    ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-rest_meth.jpg")

# Hill numbers in facet_grid
rest.meth.labels <- c(cus = "Cultivar Seed Mixture", mga = "Management Adaptation",
                      res = "Regional Seed Mixture", dih = "Direct Harvesting")
q.order.labels <- c(hill.0 = "Hill's q=0", hill.1 = "Hill's q=1", hill.2 = "Hill's q=2")
plot_div_data %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = fct_reorder(rest.meth, eff.no, median), y = eff.no, fill = rest.meth)) +
  geom_boxplot(show.legend = T) +
  labs(x = "Restoration method", y = "Effective number of species")+
  # scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols = vars(
    q.order),
    labeller = labeller(q.order = q.order.labels)
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       labels = rest.meth.labels,
                       name = "Restoration Method",
                       breaks = c("cus", "mga", "res", "dih")
                       ) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = c(0, 80))
ggsave("outputs/figures/hill_012-rest_meth_2.jpg",
       dpi = 300, width = 16.5, height = 11, units = "cm")

# Anova
rest_sites_noNA <- sites %>% 
  filter(rest.meth != "NA",
         rest.meth != "dih&seed")
## hill 0
rest_aov <- aov(hill.0 ~ rest.meth,
                data = rest_sites_noNA
)
summary(rest_aov)
# F(3,114) = 15.45; p < 0.001
shapiro.test(rest_aov$residuals)
bartlett.test(richness_all ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = hill.0 ~ rest.meth, data = rest_sites_noNA)
# 
# $rest.meth
# diff         lwr       upr     p adj
# dih-cus  20.7725753  12.4785658 29.066585 0.0000000
# mga-cus   8.4229249  -0.9850446 17.830894 0.0963213
# res-cus   9.0191816   0.5019403 17.536423 0.0335059
# mga-dih -12.3496503 -20.7613968 -3.937904 0.0011975
# res-dih -11.7533937 -19.1554872 -4.351300 0.0003849
# res-mga   0.5962567  -8.0356774  9.228191 0.9979141

## hill 1
rest_aov <- aov(hill.1 ~ rest.meth,
                data = rest_sites_noNA
)
summary(rest_aov)
# F(3,114) = 9.391; p < 0.001
shapiro.test(rest_aov$residuals)
bartlett.test(richness_all ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = hill.1 ~ rest.meth, data = rest_sites_noNA)
# 
# $rest.meth
# diff        lwr        upr     p adj
# dih-cus  9.394608  4.7027252 14.0864914 0.0000048
# mga-cus  4.816884 -0.5051614 10.1389290 0.0908782
# res-cus  4.635352 -0.1828125  9.4535161 0.0639431
# mga-dih -4.577724 -9.3362109  0.1807619 0.0639614
# res-dih -4.759256 -8.9465868 -0.5719261 0.0191354
# res-mga -0.181532 -5.0645775  4.7015135 0.9996719

## hill 2
rest_aov <- aov(hill.2 ~ rest.meth,
                data = rest_sites_noNA
)
summary(rest_aov)
# F(3,114) = 8.352; p < 0.001
shapiro.test(rest_aov$residuals)
bartlett.test(richness_all ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = hill.2 ~ rest.meth, data = rest_sites_noNA)
# 
# $rest.meth
# diff        lwr        upr     p adj
# dih-cus  6.6670415  3.1082900 10.2257931 0.0000201
# mga-cus  2.9079018 -1.1288220  6.9446257 0.2431754
# res-cus  3.3837955 -0.2707393  7.0383302 0.0801173
# mga-dih -3.7591397 -7.3684092 -0.1498701 0.0377976
# res-dih -3.2832460 -6.4592988 -0.1071932 0.0398630
# res-mga  0.4758936 -3.2278529  4.1796402 0.9869780

## only North sites
rest_sites_noNA <- sites %>% 
  filter(rest.meth != "NA",
         rest.meth != "dih&seed") %>% 
  filter(region == "north")

## hill 0
rest_aov <- aov(hill.0 ~ rest.meth,
                data = rest_sites_noNA
)
summary(rest_aov)
# F(3,36) = 3.545; p = 0.0239
shapiro.test(rest_aov$residuals)
bartlett.test(hill.0 ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# dih-cus p = 0.0199265

## only Centre sites
rest_sites_noNA <- sites %>% 
  filter(rest.meth != "NA",
         rest.meth != "dih&seed") %>% 
  filter(region == "centre")

## hill 0
rest_aov <- aov(hill.0 ~ rest.meth,
                data = rest_sites_noNA
)
summary(rest_aov)
# F(3,35) = 5.848; p = 0.00239
shapiro.test(rest_aov$residuals)
bartlett.test(hill.0 ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# dih-cus p = 0.0086
# mga-cus p = 0.0037

#### by region ####


plot_div_data %>% 
  filter(region == "north") %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = rest.meth)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols = vars(
    fct_reorder(rest.meth, eff.no, median)),
    # labeller = labeller(rest.meth = rest.meth.labels)
  ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-rest_meth.jpg")


### by region ####
plot_div_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(region)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-region.jpg")

plot_div_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = region, y = eff.no, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "region", y = "Effective number of species")+
  # scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(q.order)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-region_2.jpg")

### by age of restoration ####
#### by restoration method ####
plot_div_data %>% 
  filter(start.rest != "NA",
         rest.meth != "dih&seed",
         rest.meth != "NA") %>% 
  filter(q.order == "hill.0") %>% 
  ggplot(aes(x = start.rest, y = eff.no, colour = rest.meth)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Start of restoration", y = "Species richness q_0") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration method",
                         )
ggsave("outputs/figures/hill_0-start_method.jpg")

#### by region ####
plot_div_data %>% 
  filter(start.rest != "NA") %>% 
  filter(q.order == "hill.0") %>% 
  ggplot(aes(x = start.rest, y = eff.no, colour = region)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Start of restoration", y = "Species richness q_0") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Region",
                         breaks = c("South", "Centre", "North") # changes order of legend
  )
ggsave("outputs/figures/hill_0-start_region.jpg")

### by land use history ####
plot_div_data %>% 
  filter(land.use.hist != "NA") %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = land.use.hist)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(land.use.hist)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-land_use_history.jpg")

### by hydrology ####
#### CWM Ellenberg F abundance ####
plot_div_data %>% 
  # filter(q.order == "hill.0") %>%
  ggplot(aes(x = site_cwm_abu_oek_f, y = eff.no, colour = q.order)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "site CWM Ellenberg F", y = "Effective number of species") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "q order",
  )
ggsave("outputs/figures/hill_012-CWM_abu_oek_f.jpg")

##### and by restoration method ####
plot_div_data %>% 
  filter(rest.meth != "dih&seed",
         rest.meth != "NA") %>% 
    filter(q.order == "hill.0") %>%
  ggplot(aes(x = site_cwm_abu_oek_f, y = eff.no, colour = rest.meth)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "site CWM Ellenberg F", y = "Effective number of species") +
  coord_cartesian(ylim = c(0, 80)) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration method",
  )
ggsave("outputs/figures/hill_012-CWM_abu_oek_f_method.jpg")

#### factorial ####
plot_div_data %>% 
  # filter(hydrology != "NA") %>% 
  ggplot(aes(x = q.order, y = eff.no, fill = hydrology)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "order q", y = "Effective number of species")+
  scale_x_discrete(labels = c("0", "1", "2"))+
  facet_grid(cols=vars(hydrology)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-hydrology.jpg")


### by LUI ####
plot_div_data %>% 
  filter(site.type == "restored") %>% 
  # filter(q.order == "hill.0") %>%
  ggplot(aes(x = lui, y = eff.no, colour = q.order)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "LUI (Land Use Intensity Index)", y = "Effective number of species") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "q order",
  )
ggsave("outputs/figures/hill_012-LUI.jpg")

#### and by restoration method ####
plot_div_data %>%
  filter(site.type == "restored") %>% 
  filter(rest.meth != "dih&seed",
         rest.meth != "NA") %>% 
  filter(q.order == "hill.0") %>%
  ggplot(aes(x = lui, y = eff.no, colour = rest.meth)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "LUI (Land Use Intensity Index)", y = "Effective number of species") +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Restoration method",
  )
ggsave("outputs/figures/hill_012-LUI_method.jpg")

#### and by management type ####
plot_div_data %>% 
  filter(site.type == "restored") %>% 
  filter(rest.meth != "dih&seed",
         rest.meth != "NA") %>% 
  filter(q.order == "hill.0") %>%
  ggplot(aes(x = lui, y = eff.no, colour = mngm.type)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  labs(x = "LUI (Land Use Intensity Index)", y = "Effective number of species") +
  coord_cartesian(ylim = c(0, 80)) +
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                         name = "Management",
  )
ggsave("outputs/figures/hill_012-LUI_management.jpg")

### by management type ####
plot_div_data %>% 
  filter(site.type == "restored") %>% 
  filter(mngm.type != "check",
         mngm.type != "none",
         !is.na(mngm.type)) %>%
  ggplot(aes(x = mngm.type, y = eff.no, fill = mngm.type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "management", y = "Effective number of species")+
  facet_grid(cols=vars(q.order)) + 
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/hill_012-management.jpg")

## Hill number profiles ####

### Calculation ####
# calculate sequence of order q (0 - 5 in 0.1 steps)
hilldiv_profile <- div_profile(veg_spp_site_mat)
# transform and as tibble
hilldiv_profile_tbl <- as_tibble(t(hilldiv_profile),
                                 rownames = "site.ID")

# div_profile_data <- left_join(env_data[c("site.ID","rest.meth", "site.type")],
#                               hilldiv_profile_tbl)
env_data <- left_join(env_data, hilldiv_profile_tbl)


### profile by site type####
# calculate mean effective numbers in every q-order per site type
div_mean_site_type <- env_data %>% 
  # TO DO: resolve NAs and non-finite values (M_STN, M_STP) in Hill numbers
  filter("0" != "NA", site.ID != "M_STN", site.ID != "M_STP") %>% 
  group_by(site.type) %>% 
  summarise(across("0":"5", \(x) mean(x, na.rm = TRUE))) %>%
  # summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  pivot_longer(cols= "0":"5",
               names_to= "order.q",
               values_to = "spec.eff",
               names_transform = list(order.q = as.numeric))
# plot profile
ggplot(div_mean_site_type, aes(x = order.q, y = spec.eff, colour = site.type))+
  geom_line(linewidth = 1.5)+
  xlab("Diversity order")+
  ylab("Effective number of species")+
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.95,
                         name = "Site Type",
                         breaks = c("positive", "restored", "negative")) # changes order of legend
ggsave("outputs/figures/hill_profile-site_type.jpg")



### profile by restoration method ####
# calculate mean effective numbers in every q-order per restoration method
div_mean_rest_meth <- env_data %>% 
  # TO DO: resolve NAs and non-finite values (M_STN, M_STP) in Hill numbers
  filter("0" != "NA", site.ID != "M_STN", site.ID != "M_STP") %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  group_by(rest.meth) %>% 
  summarise(across("0":"5", \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols= "0":"5",
               names_to= "order.q",
               values_to = "spec.eff",
               names_transform = list(order.q = as.numeric))
# plot profile
ggplot(div_mean_rest_meth, aes(x = order.q, y = spec.eff, colour = rest.meth))+
  geom_line(linewidth = 1.5)+
  xlab("Diversity order")+
  ylab("Effective number of species")+
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.95,
                         name = "Restoration method",
                         breaks = c("dih", "res", "mga", "cus") # changes order of legend
                         ) 
ggsave("outputs/figures/hill_profile-rest_meth.jpg")

### profile by region ####
div_mean_region <- env_data %>% 
  # TO DO: resolve NAs and non-finite values (M_STN, M_STP) in Hill numbers
  filter("0" != "NA", site.ID != "M_STN", site.ID != "M_STP") %>% 
  group_by(region) %>% 
  summarise(across("0":"5", \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(cols= "0":"5",
               names_to= "order.q",
               values_to = "spec.eff",
               names_transform = list(order.q = as.numeric))
# plot profile
ggplot(div_mean_region, aes(x = order.q, y = spec.eff, colour = region))+
  geom_line(linewidth = 1.5)+
  xlab("Diversity order")+
  ylab("Effective number of species")+
  scale_colour_viridis_d(option="viridis", begin = 0.25, end = 0.95,
                         name = "Region",
                         breaks = c("South", "Centre", "North") # changes order of legend
  ) 
ggsave("outputs/figures/hill_profile-region.jpg")



# Species richness --------------------------------------------------------

## Calculation ####
# total species richness
n_all <- abundances %>% 
  group_by(id.site) %>% 
  count() %>% 
  rename(richness_all = n)


# add species richness to sites data
sites <- sites %>% 
  left_join(n_all, by = "id.site")

# #calculate species richness by counting cover entries
# # add column with species richness to env data
# env_data <- env_data %>% 
#   left_join(
#     cover_mean %>% 
#       group_by(site.ID) %>% 
#       count()
#   ) %>% 
#   rename(richness.pl = n)

## Plots ####

rm(list = setdiff(ls(), c("sites")))

get_box_stats <- function(y, upper_limit = stats_height * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste(
      "n =", length(y), "\n",
      "Median =", round(median(y), 2), "\n"
      
    )
  ))
}

stats_height <- max(sites$richness_all)


### plant species richness by restoration method ####
sites %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = fct_reorder(rest.meth, richness_all, median),
             y = richness_all, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Restoration method", y = "Plant Species Richness")+
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       name = "Restoration method",
                       breaks = c("dih", "res", "mga", "cus") # changes order of legend
  )
ggsave("outputs/figures/richness_pl-rest_meth.jpg")

# deutsch ###
sites %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = fct_reorder(rest.meth, richness_all, median),
             y = richness_all, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Renaturierungsmethode", y = "Anzahl Pflanzenarten")+
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       name = "Renaturierungsmethode",
                       breaks = c("cus", "res", "mga", "dih"),
                       labels = c ("Regel-Saatgut-Mischung", "Regio-Saatgut-Mischung",
                                   "Management Anpassung", "Direkternte Methoden") 
  ) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = c(0, 80))
ggsave("outputs/figures/richness_pl-rest_meth-deutsch.jpg",
       dpi = 300, width = 16.5, height = 11, units = "cm")


#### and by management type ####
rest.meth.labels <- c(cus = "Cultivar Seed Mixture", mga = "Management Adaptation",
                      res = "Regional Seed Mixture", dih = "Direct Harvesting")
sites %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  filter(mngm.type != "none") %>% 
  ggplot(aes(x = mngm.type,
             y = richness_all, fill = mngm.type)) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Restoration method", y = "Plant Species Richness") +
  facet_grid(cols = vars(rest.meth),
             labeller = labeller(rest.meth = rest.meth.labels)) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       name = "Management",
                       # breaks = c("dih", "res", "mga", "cus") # changes order of legend
  )
ggsave("outputs/figures/richness_pl-rest_meth_mngm.jpg",
       dpi = 300, width = 20, height = 11, units = "cm")




  

# Anova

rest_sites_noNA <- sites %>% 
  filter(rest.meth != "NA",
         rest.meth != "dih&seed")

rest_aov <- aov(richness_all ~ rest.meth,
               data = rest_sites_noNA
)
summary(rest_aov)
# F(3,114) = 16.85; p < 0.001
shapiro.test(rest_aov$residuals)
bartlett.test(richness_all ~ rest.meth, data = rest_sites_noNA)
TukeyHSD(rest_aov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = richness_all ~ rest.meth, data = rest_sites_noNA)
# 
# $rest.meth
# diff        lwr        upr     p adj
# dih-cus  45.67001  27.829061  63.510961 0.0000000
# mga-cus  13.19565  -7.041498  33.432802 0.3283976
# res-cus  27.07801   8.756869  45.399142 0.0010937
# mga-dih -32.47436 -50.568569 -14.380149 0.0000468
# res-dih -18.59201 -34.514388  -2.669625 0.0151389
# res-mga  13.88235  -4.685495  32.450201 0.2135558


### by region ####
sites %>% 
  ggplot(aes(x = region, y = richness_all, fill = region)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Region", y = "Plant Species Richness")+
  # facet_grid(cols=vars(system)) +
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  # theme(strip.text.x = element_text(size = 10, face = "bold")) +
  coord_cartesian(ylim = c(0, 80)) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/richness-region.jpg",
       dpi = 300, width = 12, height = 11, units = "cm")

#### by region and method ####
env_data %>% 
  # TO DO: resolve NAs and non-finite values (M_STN, M_STP) in Hill numbers
  filter("0" != "NA", site.ID != "M_STN", site.ID != "M_STP") %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = region, y = richness.pl, fill = fct_reorder(rest.meth, richness.pl, median))) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Region", y = "Plant Species Richness")+
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       name = "Restoration method",
                       breaks = c("dih", "res", "mga", "cus") # changes order of legend
  )#+
  #geom_jitter()
ggsave("outputs/figures/richness_pl-region+method.jpg")

### by site type ####
sites %>% 
  ggplot(aes(x = site.type, y = richness_all, fill = site.type)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Site Type", y = "Plant Species Richness")+
  # facet_grid(cols=vars(system)) +
  # stat_summary(
  #   fun.data = get_box_stats, geom = "text", hjust = 0.5, vjust = 0.9
  # ) +
  theme(strip.text.x = element_text(size = 10, face = "bold")) +
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97)
ggsave("outputs/figures/richness-site_type.jpg")

#### by site type and hydrology ####
env_data %>% 
  # TO DO: resolve NAs and non-finite values (M_STN, M_STP) in Hill numbers
  filter("0" != "NA", site.ID != "M_STN", site.ID != "M_STP") %>% 
  ggplot(aes(x = site.type, y = richness.pl, fill = hydrology)) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Site type", y = "Plant Species Richness")+
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       #name = "Restoration method",
                       #breaks = c("dih", "res", "mga", "cus") # changes order of legend
  )#+
  #geom_jitter()
ggsave("outputs/figures/richness_pl-site_type+hydrology.jpg")


### plant species richness by observation year ####
env_data %>% 
  ggplot(aes(x = obs.year, y = richness.pl)) +
  geom_boxplot() +
  geom_jitter()

### plant species richness by hydrology ####
env_data %>% 
  ggplot(aes(x = hydrology, y = richness.pl)) +
  geom_boxplot() +
  geom_jitter()
# heavy bias to dry sites
bartlett.test(richness.pl ~ hydrology, data = env_data)


### plant species richness by age of restoration ####
env_data %>% 
  filter(start.rest != "NA") %>% 
  ggplot(aes(x = start.rest, y = richness.pl, colour = rest.meth)) +
  #geom_boxplot() +
  geom_point()
#geom_jitter()


