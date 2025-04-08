#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Vegetation analysis
# Analysis forb grass ratio
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

# packages ----------------------------------------------------------------
library(tidyverse)
library(here)
library(vegan)
library(nlme)
library(lme4)
library(lmerTest)
library(DHARMa)
library(gridExtra)
rm(list = ls())


### Start ###
rm(list = ls())

# load data -----------------------------------------------------------------


## vegetation survey head ####
veghead_data <- readxl::read_excel(
  here(
    "data", "raw",
    "species_vegetation",
    "data_raw_plants_nms_head_20240222.xlsx"
  ), na =c("NA", "")) %>%
  rename(id.site = site.ID) %>% 
  select(id.site, subtransect, starts_with("cover"))



## site specific environment data ####
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
  mutate(rest.meth.type = if_else(
    site.type == "restored", rest.meth, site.type
  )) %>% 
  select(
    id.site, site.type, rest.meth, rest.meth.type, region, start.rest, 
    land.use.hist, hydrology, site_cwm_abu_oek_f, site_cwm_pres_oek_f
  ) %>% 
  distinct() %>% 
  # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
  filter(!is.na(id.site))

# forb-grass ratio --------------------------------------------------------
## calculation -------------------------------------------------------------


# calculate the forb grass ratio and forb+legum grass ratio at each subtransect
veghead_data <- veghead_data %>% 
  mutate(forb.grass.ratio = cover.forbs / cover.grass) %>% 
  mutate(forblegu.grass.ratio = (cover.forbs + cover.legumes) / cover.grass)

# calculate the mean forb grass ratio at each site
forb_gr_ratio <- veghead_data %>% 
  group_by(id.site) %>% 
  summarize(forb.grass.ratio.site = mean(forb.grass.ratio, na.rm = TRUE)) 


# calculate the mean forb+legumes grass ratio at each site
forblegu_gr_ratio <- veghead_data %>% 
  group_by(id.site) %>% 
  summarize(forblegu.grass.ratio.site = mean(forblegu.grass.ratio, na.rm = TRUE)) 


# add column to env data
sites <- sites %>% 
  left_join(forb_gr_ratio, by = "id.site") %>% 
  left_join(forblegu_gr_ratio, by = "id.site")




## plots -------------------------------------------------------------------
### by restoration method ####
sites %>% 
  filter(rest.meth != "NA", rest.meth != "dih&seed") %>% 
  ggplot(aes(x = fct_reorder(rest.meth, forblegu.grass.ratio.site, median),
             y = forblegu.grass.ratio.site, fill = rest.meth)) +
  geom_boxplot(show.legend = TRUE) +
  labs(x = "Restoration method", y = "Forb/grass ratio")+
  scale_fill_viridis_d(option="viridis", begin = 0.25, end = 0.97,
                       name = "Restoration method",
                       breaks = c("dih", "res", "mga", "cus") # changes order of legend
  )
ggsave("outputs/figures/forb_grass_ratio-rest_meth.jpg")






# models ------------------------------------------------------------------



## model forb-grass ratio --------------------------------------------------


## test normal distribution
hist(env_data$forblegu.grass.ratio.site)
# -> looks bad
hist(log(env_data$forblegu.grass.ratio.site))
# --> looks better
hist(sqrt(env_data$forblegu.grass.ratio.site))
shapiro.test(log(env_data$forblegu.grass.ratio.site))
shapiro.test(sqrt(env_data$forblegu.grass.ratio.site))
# sg. -> not ok


## test normal distribution
hist(env_data$forb.grass.ratio.site)
# -> looks bad
hist(log(env_data$forb.grass.ratio.site))
# --> looks better
hist(sqrt(env_data$forb.grass.ratio.site))
shapiro.test(log(env_data$forb.grass.ratio.site))
shapiro.test(sqrt(env_data$forb.grass.ratio.site))
# sg. -> not ok

# standardise explanatory variable (only numerical variables)
# so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”)
# It ensures that the estimated coefficients are all on the same scale, making it easier to compare effect sizes.



lm_forbgrass1 <- lm(forblegu.grass.ratio.site ~ 
               rest.cus + 
               rest.res + 
               rest.dih + 
               rest.mga +
               rest.cus*rest.dih + 
               rest.res*rest.dih, 
             data = env_data,
             na.action = na.omit
)
summary(lm_forbgrass1)


## model check ----------------------------------------------------------------


# check assumptions
plot(lm_forbgrass1, which = 1) #plot residuals
plot(lm_forbgrass1, which = 2) #plot qqplot
# --> doesn't look that good

# independence of variables?
env_data %>% 
  filter(region != "NA") %>% 
  ggplot(aes(x = region, y = forblegu.grass.ratio.site)) +
  geom_boxplot() +
  coord_cartesian(
    ylim = c(0,5)
  ) +
  geom_jitter()
# not that much of a difference between the regions

simulationOutput <- simulateResiduals(lm_forbgrass1, plot = TRUE)
par(mfrow = c(2, 2))
data <- env_data %>% 
  filter(rest.meth != "NA")
plotResiduals(main = "region", simulationOutput$scaledResiduals, data$region)
plotResiduals(main = "rest.cus", simulationOutput$scaledResiduals, data$rest.cus)
plotResiduals(main = "rest.res", simulationOutput$scaledResiduals, data$rest.res)
plotResiduals(main = "rest.dih", simulationOutput$scaledResiduals, data$rest.dih)
plotResiduals(main = "rest.mga", simulationOutput$scaledResiduals, data$rest.mga)


## plot -----------------------------------------------------------------------


env_data %>% 
  ggplot(aes(x = hydrology, y = forblegu.grass.ratio.site)) +
  geom_boxplot() +
  geom_jitter()
# bias to dry sites --> random factor



