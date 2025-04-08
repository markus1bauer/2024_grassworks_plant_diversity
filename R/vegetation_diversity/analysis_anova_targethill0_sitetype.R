#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# Target species richness ~ site type
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(ggpubr)
library(rstatix)



### Start ###
rm(list = ls())

## load data -------------------------------------------------------------------

# data = data
# groups = rest.meth.type
# variable = target.hill.0

## sites data
sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20240813.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, rest.meth,
    # region, rest.age, 
    # land.use.hist, hydrology, hydrology.cont, site.cwm.abu.oek.f, lui,
    # mngm.type, obs.year,
    # c.n, c.perc, n.perc, ph.value, tic.perc, toc.perc,
    # clay.perc, silt.perc, sand.perc,
  ) %>% 
  distinct()

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



## Diversity data
diversity <- read_csv(
  here("data", "processed", "data_processed_plants_site_diversity_20240814.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))

# add to environment data
data <- sites %>% 
  left_join(diversity, by = "id.site")


### Visualization ###
data %>% 
  ggplot(aes(x = site.type,
             y = target.hill.0)) +
  geom_boxplot()

# relevel according to type and median of rest.meth
data <- data %>% 
  mutate(site.type = fct_relevel(site.type, "negative", "restored", "positive"))

### summary statistics ###
data %>%
  group_by(site.type) %>%
  get_summary_stats(target.hill.0, type = "full")
# site.type variable          n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <fct>     <fct>         <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 negative  target.hill.0    32     8    31     14  10.8    22  11.2  5.93  16.2  6.57  1.16  2.37
# 2 restored  target.hill.0   122     6    65     31  23      38  15   11.9   31.3 12.3   1.11  2.20
# 3 positive  target.hill.0    33    19    62     34  28      41  13   10.4   36.4 12.1   2.11  4.29

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, site.type, target.hill.0) %>% 
  group_by(site.type) %>%
  identify_outliers(target.hill.0)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(target.hill.0 ~ site.type, data = data, na.action = na.omit)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, as all the points fall approximately along the reference line,
# we can assume normality. This conclusion is supported by the Shapiro-Wilk test.
# The p-value is not significant (p = 0.581), so we can assume normality.

# Check normality assumption by groups
data %>%
  filter(!is.na(site.type)) %>% 
  group_by(site.type) %>%
  shapiro_test(target.hill.0)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "target.hill.0", facet.by = "site.type")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(target.hill.0 ~ site.type)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(target.hill.0 ~ site.type)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd      F        p p<.05   ges
# 1 site.type   2 184 28.883 1.23e-11     * 0.239

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(2, 184) = 28.883, p < 0.001, eta2[g] = 0.239



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(target.hill.0 ~ site.type)
pwc
# term      group1   group2   null.value estimate conf.low conf.high    p.adj p.adj.signif
# * <chr>     <chr>    <chr>         <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#   1 site.type negative restored          0    15.1     9.73       20.5 1.07e- 9 ****        
#   2 site.type negative positive          0    20.1    13.4        26.9 9.06e-11 ****        
#   3 site.type restored positive          0     5.03   -0.293      10.3 6.84e- 2 ns    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "site.type")
ggboxplot(data, x = "site.type", y = "target.hill.0") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# they can be wrong positioned if the order is changed manually



# and beautiful

# analysis of variance
anova <- aov(target.hill.0 ~ site.type, data = data)
# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey, reversed = T)

# table with factors and median
dt <- data %>% 
  group_by(site.type) %>% 
  get_summary_stats(target.hill.0, type = "median") %>% 
  # arrange(median)
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$site.type)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "site.type", y = "target.hill.0",
          fill = "site.type",
          legend = "none",
) +
  labs(
    x = "Site type",
    y = "Characteristic Species Richness (16 m²)",
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc),
  ) +
  xlab(NULL) +
  # theme(
  # axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("Negative Reference",
                              "Restored",
                              "Positive Reference")) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(40, 40, 40, 40)) +
  rotate_x_text(45) +
  geom_text(
    data = dt,
    aes(x = site.type, y = 80, label = cld), vjust = -0.5) +
  scale_fill_manual(
    values = c("darkgrey", "darkgrey", "darkgrey", "darkgrey"),
    name = "Site type",
    labels = c("Negative Reference",
               "Cultivar Seed Mixture",
               "Regional Seed Mixture",
               "Management Adaptation",
               "Direct Harvesting",
               "Positive Reference"))

ggsave("outputs/figures/plants_species_diversity/anova_targethill0-sitetype.jpg",
       dpi = 300, width = 16, height = 20, units = "cm")




