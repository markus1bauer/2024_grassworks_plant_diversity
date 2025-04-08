#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# PCA soil variables
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

### Packages ###

library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())

## load data -------------------------------------------------------------------
### site environment data ####
sites <- read_csv(
  here::here("data", "processed", "sites_processed_environment_nms.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) %>%
  # calculate mean of CWM Ellenberg F per site
  group_by(id.site) %>% 
  mutate(site_cwm_abu_oek_f = mean(cwm_abu_oek_f),
         site_cwm_pres_oek_f = mean(cwm_pres_oek_f),
         c.n.B10 = c.perc.B10 / n.perc.B10,
         c.n.B30 = c.perc.B30 / n.perc.B30,
         c.n = mean(c(c.perc.B10, c.perc.B30)),
         ph.value = mean(c(ph.value.B10, ph.value.B30)),
         c.perc = mean(c(c.perc.B10, c.perc.B30)),
         tic.perc = mean(c(tic.perc.B10, tic.perc.B30)),
         toc.perc = mean(c(toc.perc.B10, toc.perc.B30)),
         n.perc = mean(c(n.perc.B10, n.perc.B30)),
         hydrology_cont = if_else(
           hydrology == "dry", 0.33, if_else(
             hydrology == "fresh", 0.66, if_else(
               hydrology == "moist", 1, -99
             )
           ))) %>%
  ungroup() %>%
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"),
         region = fct_relevel(region, "north", "centre", "south")) %>% 
  # add a variable with restoration method and site type
  mutate(rest.meth.type = if_else(
    site.type == "restored", rest.meth, site.type
  )) %>% 
  select(
    id.site, lui, ph.value.B10, ph.value.B30, ph.value, c.perc.B10, c.perc.B30,
    c.perc, tic.perc.B10, tic.perc.B30, tic.perc, toc.perc.B10, toc.perc.B30,
    toc.perc, n.perc.B10, n.perc.B30, n.perc, soil.texture.B10, soil.texture.B30,
    site_cwm_abu_oek_f, site_cwm_pres_oek_f,
    c.n.B10, c.n.B30, c.n, hydrology, hydrology_cont
  ) %>%
  distinct() %>%
  # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
  filter(!is.na(id.site))

# calculate mean proportions of clay and silt from soil texture
## source: Ad-hoc-AG Boden (2005): Bodenkundliche Kartieranleitung, 
## 5. Auflage, 438 S., 41 Abb., 103 Tab. ., 31 Listen. BGR, Hannover.


clay.valu <- c("Lehm Lts"= mean(c(25,45)),
               "lehmiger Sand" = mean(c(5,17)),
               "lehmiger Sand Sl" = mean(c(5,17)),
               "lehmiger Ton Tl" = mean(c(45,65)), 
               "Reinsand"= mean(c(0,5)), 
               "sandiger Lehm Ls3" = mean(c(17,25)), 
               "sandiger Schluff" = mean(c(0,8)),
               "sandiger Schluff Us" = mean(c(0,8)),
               "sandiger Ton" = mean(c(25,65)), 
               "sandiger Ton Ts" = mean(c(25,65)), 
               "Schluff" = mean(c(0,8)),
               "Schluff Uu"= mean(c(0,8)),
               "schluffiger Lehm Lu"= mean(c(17,30)),
               "schluffiger Sand" = mean(c(0,8)),
               "schluffiger Sand Su" = mean(c(0,8)),
               "schluffiger Ton Tu"= mean(c(25,65)), 
               "toniger Lehm" = mean(c(25,45)),
               "toniger Lehm Lt" = mean(c(25,45)), 
               "toniger Sand" =mean(c(5,25)),
               "toniger Sand St" = mean(c(5,25)), 
               "toniger Schluff" = mean(c(8,25)),
               "toniger Schluff Ut" = mean(c(8,25)),
               "sandiger Lehm" = mean(c(17,35)))

silt.valu <- c("Lehm Lts" = mean(c(15,30)),
               "lehmiger Sand" = mean(c(10,40)),
               "lehmiger Sand Sl" = mean(c(10,40)),
               "lehmiger Ton Tl" = mean(c(15,30)), 
               "Reinsand" = mean(c(0,10)),
               "sandiger Lehm Ls3" = mean(c(30,40)), 
               "sandiger Schluff" = mean(c(50,80)),
               "sandiger Schluff Us" = mean(c(50,80)),
               "sandiger Ton" = mean(c(0,15)),
               "sandiger Ton Ts"= mean(c(0,15)),
               "Schluff" = mean(c(80,100)),
               "Schluff Uu" = mean(c(80,100)),
               "schluffiger Lehm Lu" = mean(c(50,65)),
               "schluffiger Sand" = mean(c(10,50)),
               "schluffiger Sand Su" = mean(c(10,50)),
               "schluffiger Ton Tu" = mean(c(30,75)),
               "toniger Lehm" = mean(c(30,50)),
               "toniger Lehm Lt" = mean(c(30,50)),
               "toniger Sand" = mean(c(0,15)),
               "toniger Sand St" = mean(c(0,15)),
               "toniger Schluff" = mean(c(65,92)),
               "toniger Schluff Ut" = mean(c(65,92)),
               "sandiger Lehm" = mean(c(15,50)))

sand.valu <- c("Lehm Lts" = mean(c(25,60)),
               "lehmiger Sand" = mean(c(43,85)),
               "lehmiger Sand Sl" = mean(c(43,85)),
               "lehmiger Ton Tl" = mean(c(5,40)),
               "Reinsand" = mean(c(85,100)),
               "sandiger Lehm Ls3" = mean(c(35,53)),
               "sandiger Schluff" = mean(c(12,50)),
               "sandiger Schluff Us" = mean(c(12,50)),
               "sandiger Ton" = mean(c(20,75)),
               "sandiger Ton Ts" = mean(c(20,75)),
               "Schluff" = mean(c(0,20)),
               "Schluff Uu" = mean(c(0,20)), 
               "schluffiger Lehm Lu" = mean(c(5,33)),
               "schluffiger Sand" = mean(c(42,90)),
               "schluffiger Sand Su" = mean(c(42,90)),
               "schluffiger Ton Tu" = mean(c(0,25)), 
               "toniger Lehm" = mean(c(5,45)),
               "toniger Lehm Lt" = mean(c(5,45)),
               "toniger Sand" = mean(c(60,95)), 
               "toniger Sand St" = mean(c(60,95)), 
               "toniger Schluff" = mean(c(0,27)),
               "toniger Schluff Ut" = mean(c(0,27)),
               "sandiger Lehm" = mean(c(15,68)))

library(plyr)
data <- sites %>% 
  mutate(clay.perc.B10 = revalue(soil.texture.B10, clay.valu),
         clay.perc.B30 = revalue(soil.texture.B30, clay.valu),
         silt.perc.B10 = revalue(soil.texture.B10, silt.valu),
         silt.perc.B30 = revalue(soil.texture.B30, silt.valu),
         sand.perc.B10 = revalue(soil.texture.B10, sand.valu),
         sand.perc.B30 = revalue(soil.texture.B30, sand.valu)) %>% 
  mutate_at(c("clay.perc.B10", "clay.perc.B30",
              "silt.perc.B10", "silt.perc.B30",
              "sand.perc.B10", "sand.perc.B30"), as.numeric) %>% 
  mutate(clay.perc = rowMeans(select(., c(clay.perc.B10, clay.perc.B30))),
         silt.perc = rowMeans(select(., c(silt.perc.B10, silt.perc.B30))),
         sand.perc = rowMeans(select(., c(sand.perc.B10, sand.perc.B30))))
 
sites <- data




data <- sites %>%
  select(
    id.site, 
    # ph.value.B10, 
    # tic.perc.B10, 
    # toc.perc.B10, 
    # n.perc.B10, 
    # c.n.B10, 
    # clay.perc.B10,
    # silt.perc.B10,
    # sand.perc.B10,
    # ph.value.B30, 
    # tic.perc.B30, 
    # toc.perc.B30, 
    # n.perc.B30, 
    # c.n.B30, 
    # clay.perc.B30,
    # silt.perc.B30,
    # sand.perc.B30,
    ph.value, 
    tic.perc, 
    toc.perc, 
    n.perc, 
    # c.n, 
    clay.perc,
    # silt.perc,
    sand.perc,
    site_cwm_abu_oek_f,
    # site_cwm_pres_oek_f,
    # hydrology,
    hydrology_cont
  ) %>%
  filter(!is.na(site_cwm_abu_oek_f)) %>%
  column_to_rownames(var = "id.site")

#### * PCA ####
pca <- rda(data, scale = TRUE)

summary(pca, axes = 0)
# cumulative proportion explained of two first axes: 0.575 = 57,5 %
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))

# Plots using biplot.rda 
par(mfrow = c(1, 2)) 
biplot(pca, scaling = 1, type = c("text", "points"), main = "PCA - scaling 1") 
biplot(pca, type = c("text", "points"), main = "PCA - scaling 2") # Default scaling 2





min(data$clay.perc.B30)
levels(as.factor(sites$soil.texture.B30))
mean(c(0, 5))
> levels(as.factor(sites$soil.texture.B10))
[1] "Lehm Lts"            "lehmiger Sand"       "lehmiger Sand Sl"    "lehmiger Ton Tl"     "Reinsand"           
[6] "sandiger Lehm Ls3"   "sandiger Schluff"    "sandiger Schluff Us" "sandiger Ton"        "sandiger Ton Ts"    
[11] "Schluff"             "Schluff Uu"          "schluffiger Lehm Lu" "schluffiger Sand"    "schluffiger Sand Su"
[16] "schluffiger Ton Tu"  "toniger Lehm"        "toniger Lehm Lt"     "toniger Sand"        "toniger Sand St"    
[21] "toniger Schluff"     "toniger Schluff Ut" 