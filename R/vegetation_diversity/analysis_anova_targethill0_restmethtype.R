#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# Characteristic Species Richness ~ restmeth-sitetype
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
  here("data", "processed", "sites_processed_environment_nms_20240923.csv"),
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
  here("data", "processed", "data_processed_plants_site_diversity_20241111_agg.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  ))

# add to environment data
data <- sites %>% 
  left_join(diversity, by = "id.site")


### Visualization ###
data %>% 
  ggplot(aes(x = rest.meth.type,
             y = target.hill.0)) +
  geom_boxplot()

# relevel according to type and median of rest.meth
data <- data %>% 
  mutate(rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "mga", "dih", "positive"))

### summary statistics ###
data %>%
  group_by(rest.meth.type) %>%
  get_summary_stats(target.hill.0, type = "full")
# rest.meth.type variable          n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <fct>          <fct>         <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 negative       target.hill.0    32     8    31   14    10.8  22    11.2  5.93  16.2  6.57  1.16  2.37
# 2 cus            target.hill.0    21     8    41   19    15    25    10    7.41  21.6  8.61  1.88  3.92
# 3 res            target.hill.0    38     6    60   29.5  23    35.8  12.8  9.64  29.8  9.66  1.57  3.17
# 4 mga            target.hill.0    22    10    51   29    22.2  35.2  13   10.4   28.9  9.53  2.03  4.23
# 5 dih            target.hill.0    40    16    65   37.5  31    47.2  16.2 13.3   39.7 12.4   1.96  3.97
# 6 positive       target.hill.0    33    19    62   34    28    41    13   10.4   36.4 12.1   2.11  4.29


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, rest.meth.type, target.hill.0) %>% 
  group_by(rest.meth.type) %>%
  identify_outliers(target.hill.0)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(target.hill.0 ~ rest.meth.type, data = data, na.action = na.omit)
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
  shapiro_test(target.hill.0)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "target.hill.0", facet.by = "rest.meth.type")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(target.hill.0 ~ rest.meth.type)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(target.hill.0 ~ rest.meth.type)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd      F        p p<.05   ges
# 1 rest.meth.type   5 180 24.199 1.37e-18     * 0.402

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(5, 189) = 22.745, p < 0.001, eta2[g] = 0.387.



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(target.hill.0 ~ rest.meth.type)
pwc
# term           group1   group2   null.value estimate conf.low conf.high    p.adj p.adj.signif
# * <chr>          <chr>    <chr>         <dbl>    <dbl>    <dbl>     <dbl>    <dbl> <chr>       
#   1 rest.meth.type negative cus               0    5.40    -2.88      13.7  4.19e- 1 ns          
#   2 rest.meth.type negative res               0   13.5      6.47      20.6  1.79e- 6 ****        
#   3 rest.meth.type negative mga               0   12.7      4.52      20.9  1.94e- 4 ***         
#   4 rest.meth.type negative dih               0   23.5     16.5       30.5  0        ****        
#   5 rest.meth.type negative positive          0   20.1     12.8       27.5  3.32e-12 ****        
#   6 rest.meth.type cus      res               0    8.14     0.124     16.2  4.43e- 2 *           
#   7 rest.meth.type cus      mga               0    7.29    -1.71      16.3  1.86e- 1 ns          
#   8 rest.meth.type cus      dih               0   18.1     10.2       26.1  8.2 e- 9 ****        
#   9 rest.meth.type cus      positive          0   14.7      6.51      23.0  9.62e- 6 ****        
#   10 rest.meth.type res      mga               0   -0.854   -8.76       7.05 1   e+ 0 ns          
#   11 rest.meth.type res      dih               0    9.96     3.28      16.6  4.06e- 4 ***         
#   12 rest.meth.type res      positive          0    6.60    -0.418     13.6  7.84e- 2 ns          
#   13 rest.meth.type mga      dih               0   10.8      2.99      18.6  1.38e- 3 **          
#   14 rest.meth.type mga      positive          0    7.45    -0.664     15.6  9.19e- 2 ns          
#   15 rest.meth.type dih      positive          0   -3.36   -10.3        3.58 7.29e- 1 ns      

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "rest.meth.type")
ggboxplot(data, x = "rest.meth.type", y = "target.hill.0") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# they can be wrong positioned if the order is changed manually



# and beautiful

# analysis of variance
anova <- aov(target.hill.0 ~ rest.meth.type, data = data)
# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey, reversed = T)

# table with factors and median
dt <- data %>% 
  group_by(rest.meth.type) %>% 
  get_summary_stats(target.hill.0, type = "median") %>% 
  # arrange(median)
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$rest.meth.type)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "rest.meth.type", y = "target.hill.0",
          fill = "rest.meth.type",
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
    aes(x = rest.meth.type, y = 80, label = cld), vjust = -0.2) +
  scale_fill_manual(
    values = c("lightgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "lightgrey"),
    name = "Site type",
    labels = c("Negative Reference",
               "Cultivar Seed Mixture",
               "Regional Seed Mixture",
               "Management Adaptation",
               "Direct Harvesting",
               "Positive Reference"))

ggsave("outputs/figures/plants_species_diversity/anova_targethill0-meth&type.jpg",
       dpi = 300, width = 16, height = 20, units = "cm")




