#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# Red List Germany Species Richness ~ restmeth-sitetype
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
# variable = rlg.hill.0

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


sites <- sites %>% 
  # add a variable with restoration method and site type
  mutate(rest.meth.type = if_else(
    site.type == "restored", rest.meth, site.type)) %>% 
  filter(!is.na(rest.meth.type))


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
  ggplot(aes(x = rest.meth.type,
             y = rlg.hill.0)) +
  geom_boxplot()

# relevel according to type and median of rest.meth
data <- data %>% 
  mutate(rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "mga", "dih", "positive"))

### summary statistics ###
data %>%
  group_by(rest.meth.type) %>%
  get_summary_stats(rlg.hill.0, type = "full")
# rest.meth.type variable       n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <fct>          <fct>      <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 negative       rlg.hill.0    32     0     4      0     0  1     1     0    0.594  1.10 0.195 0.398
# 2 cus            rlg.hill.0    21     0     2      0     0  1     1     0    0.476  0.68 0.148 0.309
# 3 res            rlg.hill.0    38     0     8      1     0  2.75  2.75  1.48 1.74   2.16 0.351 0.711
# 4 mga            rlg.hill.0    22     0    29      3     1 10.8   9.75  4.45 6.54   7.71 1.64  3.42 
# 5 dih            rlg.hill.0    40     0    14      3     1  5     4     2.96 3.85   3.77 0.595 1.20 
# 6 positive       rlg.hill.0    33     0    29      3     1  6     5     2.96 5      6.60 1.15  2.34 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, rest.meth.type, rlg.hill.0) %>% 
  group_by(rest.meth.type) %>%
  identify_outliers(rlg.hill.0)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(rlg.hill.0 ~ rest.meth.type, data = data, na.action = na.omit)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, as all the points fall approximately along the reference line,
# we can assume normality. This conclusion is supported by the Shapiro-Wilk test.
# The p-value is not significant (p = 0.581), so we can assume normality.

# Check normality assumption by groups
data %>%
  filter(!is.na(rest.meth.type)) %>% 
  group_by(rest.meth.type) %>%
  shapiro_test(rlg.hill.0)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "rlg.hill.0", facet.by = "rest.meth.type")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(rlg.hill.0 ~ rest.meth.type)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(rlg.hill.0 ~ rest.meth.type)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd     F        p p<.05   ges
# 1 rest.meth.type   5 180 8.614 2.45e-07     * 0.193

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(5, 189) = 22.745, p < 0.001, eta2[g] = 0.387.



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(rlg.hill.0 ~ rest.meth.type)
pwc
# term           group1   group2   null.value estimate conf.low conf.high     p.adj p.adj.signif
# * <chr>          <chr>    <chr>         <dbl>    <dbl>    <dbl>     <dbl>     <dbl> <chr>       
# 1 rest.meth.type negative cus               0   -0.118 -3.64        3.41  1         ns          
# 2 rest.meth.type negative res               0    1.14  -1.87        4.15  0.883     ns          
# 3 rest.meth.type negative mga               0    5.95   2.48        9.43  0.0000268 ****        
# 4 rest.meth.type negative dih               0    3.26   0.280       6.23  0.023     *           
# 5 rest.meth.type negative positive          0    4.41   1.29        7.52  0.000953  ***         
# 6 rest.meth.type cus      res               0    1.26  -2.15        4.67  0.895     ns          
# 7 rest.meth.type cus      mga               0    6.07   2.24        9.90  0.000132  ***         
# 8 rest.meth.type cus      dih               0    3.37  -0.00738     6.76  0.0509    ns          
# 9 rest.meth.type cus      positive          0    4.52   1.02        8.03  0.00356   **          
# 10 rest.meth.type res      mga               0    4.81   1.45        8.17  0.000806  ***         
# 11 rest.meth.type res      dih               0    2.11  -0.729       4.96  0.271     ns          
# 12 rest.meth.type res      positive          0    3.26   0.278       6.25  0.0232    *           
# 13 rest.meth.type mga      dih               0   -2.70  -6.03        0.635 0.187     ns          
# 14 rest.meth.type mga      positive          0   -1.55  -5.00        1.91  0.791     ns          
# 15 rest.meth.type dih      positive          0    1.15  -1.80        4.10  0.871     ns  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "rest.meth.type")
ggboxplot(data, x = "rest.meth.type", y = "rlg.hill.0") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# they can be wrong positioned if the order is changed manually



# and beautiful

# analysis of variance
anova <- aov(rlg.hill.0 ~ rest.meth.type, data = data)
# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey, reversed = T)

# table with factors and median
dt <- data %>% 
  group_by(rest.meth.type) %>% 
  get_summary_stats(rlg.hill.0, type = "median") %>% 
  # arrange(median)
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$rest.meth.type)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "rest.meth.type", y = "rlg.hill.0",
          fill = "rest.meth.type",
          legend = "none",
) +
  labs(
    x = "Site type",
    y = "Red List Germany Species Richness (16 m²)",
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc),
  ) +
  xlab(NULL) +
  # theme(
  # axis.text.x = element_blank()) +
  scale_x_discrete(labels = c("Negative Reference",
                              "Cultivar Seed Mixture",
                              "Regional Seed Mixture",
                              "Management Adaptation",
                              "Direct Harvesting",
                              "Positive Reference")) +
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin = margin(40, 40, 40, 40)) +
  rotate_x_text(45) +
  geom_text(
    data = dt,
    aes(x = rest.meth.type, y = 30, label = cld), vjust = -0.2) +
  scale_fill_manual(
    values = c("lightgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "lightgrey"),
    name = "Site type",
    labels = c("Negative Reference",
               "Cultivar Seed Mixture",
               "Regional Seed Mixture",
               "Management Adaptation",
               "Direct Harvesting",
               "Positive Reference")) +
  coord_cartesian(ylim = c(0,30))

ggsave("outputs/figures/plants_species_diversity/anova_rlghill0-meth&type.jpg",
       dpi = 300, width = 16, height = 20, units = "cm")




