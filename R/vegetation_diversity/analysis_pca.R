#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# PCA all variables
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/115-famd-factor-analysis-of-mixed-data-in-r-essentials/

### Packages ###

library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())

## load data -------------------------------------------------------------------
### site environment data ####
sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20240813.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))

# todo: resolve infinite numbers in c.n.B30 and c.n in sites N_RET and N_HOR
# quick and dirty:
sites[sites$id.site == "N_RET", "c.n"] <- NA
sites[sites$id.site == "N_HOR", "c.n"] <- NA
sites[sites$id.site == "N_RET", "c.n.B30"] <- NA
sites[sites$id.site == "N_HOR", "c.n.B30"] <- NA

# choose one restoration method for "dih&seed"
# --> dih because it is most likely that dih has the greatest influence
# --> Line said for M_SKF: ReS
sites %>% 
  filter(rest.meth == "dih&seed")
sites[sites$id.site == "S_FHZ", "rest.meth"] <- "dih"
sites[sites$id.site == "S_GIG", "rest.meth"] <- "dih"
sites[sites$id.site == "M_SKF", "rest.meth"] <- "res"
sites[sites$id.site == "S_WTZ", "rest.meth"] <- "dih"

### diversity data ####

diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20240814.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))


data <- sites %>% 
  left_join(diversity, by = "id.site") %>% 
  select(
    id.site,
    tot.hill.0,
    tot.hill.1,
    tot.hill.2,
    target.hill.0,
    target.hill.1,
    target.hill.2,
    fcsi.hill.0,
    # hydrology,
    # region,
    # rest.meth,
    rest.age,
    # land.use.hist,
    # site.type,
    # site.cwm.abu.oek.f,
    lui,
    # mngm.type,
    c.n,
    c.perc,
    n.perc,
    ph.value,
    tic.perc,
    toc.perc,
    clay.perc,
    silt.perc,
    sand.perc
  ) %>% 
  distinct() %>% 
  column_to_rownames(var = "id.site")

data <- na.omit(data) 

  
#### * PCA ####
pca <- rda(data, scale = TRUE)

summary(pca, axes = 0)
# cumulative proportion explained of two first axes: 0.5596 = 56.0 %
screeplot(pca, bstick = TRUE, npcs = length(pca$CA$eig))

# Plots using biplot.rda 
par(mfrow = c(1, 2)) 
biplot(pca, scaling = 1, type = c("text", "points"), main = "PCA - scaling 1") 
biplot(pca, type = c("text", "points"), main = "PCA - scaling 2") # Default scaling 2





# FAMD --------------------------------------------------------------------

library("FactoMineR")
library("factoextra")

data <- sites %>% 
  left_join(diversity, by = "id.site") %>% 
  select(
    id.site,
    # tot.hill.0,
    # tot.hill.1,
    # tot.hill.2,
    # target.hill.0,
    # target.hill.1,
    # target.hill.2,
    # fcsi.hill.0,
    hydrology,
    region,
    rest.meth,
    rest.age,
    land.use.hist,
    # site.type,
    site.cwm.abu.oek.f,
    lui,
    mngm.type,
    c.n,
    c.perc,
    n.perc,
    ph.value,
    tic.perc,
    toc.perc,
    clay.perc,
    silt.perc,
    sand.perc
  ) %>% 
  distinct() %>% 
  column_to_rownames(var = "id.site")
data <- na.omit(data) 



# compute FAMD
res.famd <- FAMD(data, graph = FALSE)

# proportion of variances by dimensions
eig.val <- get_eigenvalue(res.famd)
head(eig.val)

# scree plot (the percentages of inertia explained by each FAMD dimension)
fviz_screeplot(res.famd)

## all variables ###

# Plot of variables
fviz_famd_var(res.famd, repel = TRUE)
# Contribution to the first dimension
fviz_contrib(res.famd, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.famd, "var", axes = 2)

## quantitative variables ###

# Plot of variables
fviz_famd_var(res.famd, "quanti.var", repel = TRUE,
              col.var = "black")
# highlight contribution of variables
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

## qualitative variables ###

# Plot of variables
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

fviz_famd_ind(res.famd, col.ind = "cos2", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

fviz_mfa_ind(res.famd, 
             habillage = "rest.meth", # color by groups 
             # palette = c("#00AFBB", "#E7B800", "#FC4E07", "lightgreen"),
             addEllipses = TRUE, ellipse.type = "confidence", 
             repel = TRUE # Avoid text overlapping
) 
fviz_ellipses(res.famd, c("rest.meth", "hydrology"), repel = TRUE)
