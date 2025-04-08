#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Plot figure: Species abundance accumulation curve
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(tidyverse)
library(here)
library(vegan)
library(BiodiversityR)


### Start ###
rm(list = ls())


## load data -------------------------------------------------------------------


### species abundance data ####
abundances <- read_csv(
  here(
    "data", "processed", "data_processed_plants_nms_site_abundances.csv"
  ),
  col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
) %>% 
  # todo: check name input
  rename(name.plant = name.tnrs)


abundances_plot <- read_csv(
  here(
    "data", "processed", "data_processed_plants_nms_plot_abundances.csv"
  ),
  col_names = TRUE, na = c("", "NA", "na"), col_types = cols(.default = "?")
) %>% 
  # todo: check name input
  rename(name.plant = name.original)

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
  mutate(rest.meth.type = if_else(
    site.type == "restored", rest.meth, site.type
  )) %>% 
  select(
    id.site, id.plot, site.type, rest.meth, rest.meth.type, region, start.rest, 
    land.use.hist, hydrology, site_cwm_abu_oek_f, site_cwm_pres_oek_f
  ) %>% 
  # distinct() %>% 
  # remove row with no values (only NAs) --> should be resolved when M_WDG issue is gone
  filter(!is.na(id.site))


## transform data -------------------------------------------------------------------
abundances_wide <- abundances %>% 
  pivot_wider(
    id_cols = id.site,
    names_from = name.plant,
    values_from = cover.mean,
    values_fill = 0
  ) #%>% 
  # column_to_rownames("id.site")

abundances_wide <- abundances_plot %>% 
  filter(subtransect != "T") %>% 
  # quick and dirty remove duplicates 
  # todo: solve later if issue is not gone with German SL name resolving
  # distinct(id.plot, name.plant) %>% 
    pivot_wider(
    id_cols = id.plot,
    names_from = name.plant,
    values_from = presence,
    values_fill = 0
  )


  


## combine site information to abundance data
sites_subset <- sites %>% 
  select(id.plot, region)
abundances_wide <- abundances_wide %>% 
  left_join(sites_subset, by = "id.plot") %>% 
  select(id.plot, region, everything())

# subsets
abundances_wide_N <- abundances_wide %>% 
  filter(region == "north")
abundances_wide_C <- abundances_wide %>% 
  filter(region == "centre")
abundances_wide_S <- abundances_wide %>% 
  filter(region == "south")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Calculation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


curve_all <- specaccum(abundances_wide[, 3:680])
curve_N <- specaccum(abundances_wide_N[, 3:680])
curve_C <- specaccum(abundances_wide_C[, 3:680])
curve_S <- specaccum(abundances_wide_S[, 3:680])

curve_all <- specaccum(abundances_wide[, 2:680])
curve_all <- accumresult(abundances_wide[, 2:680])


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Plots #####################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot(curve_all)

plot(curve_all,
     col="grey70", lwd=2,
     ci.lty=0,
     # xlim=c(1,744), ylim=c(1,700),
     main="Species accumulation",
     ylab="Number of species",
)
plot(curve_N, 
     add = T,
     col="blue", lwd=2,
     ci.lty=0,
)
plot(curve_C, 
     add = T,
     col="red", lwd=2,
     ci.lty=0,
    )
plot(curve_S, 
     add = T,
     col="green", lwd=2,
     ci.lty=0,
     )


curve_all_pool <- poolaccum(abundances_wide)

jpeg(file = "outputs/figures/species_accumulation_chao.jpg")
plot(curve_all_pool, display = "chao")
dev.off()


ggsave(
  here("outputs", "figures", "figure_species_accumulation.jpg"),
  dpi = 300, width = 16.5, height = 12, units = "cm"
)


abundances


veg.list <- list(
  Beaver = t(beaver[,4:ncol(beaver)]), 
  Control = t(control[,4:ncol(control)]))


# Loading data
data(BCI)

# Using vegan package
sp1 <- specaccum(BCI[1:25,], method = "exact")
plot(sp1)
sp2 <- specaccum(BCI, method = "random")

# Using BiodiversityR package
sp3 <- accumresult(BCI, method = "exact")
sp4 <- accumresult(BCI, method = "random")
