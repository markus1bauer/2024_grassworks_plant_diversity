#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# Species Richness ~ restmeth-sitetype
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
# variable = tot.hill.0

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
             y = tot.hill.0)) +
  geom_boxplot()

# relevel according to type and median of rest.meth
data <- data %>% 
  mutate(rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "mga", "dih", "positive"))

### summary statistics ###
data %>%
  group_by(rest.meth.type) %>%
  get_summary_stats(tot.hill.0, type = "full")
# rest.meth.type variable       n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <fct>          <fct>      <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 negative       tot.hill.0    32    11    47   20    15.8  26.5  10.8  8.15  21.9  8.55  1.51  3.08
# 2 cus            tot.hill.0    21    13    45   29    24    33     9    5.93  28.2  8.84  1.93  4.02
# 3 res            tot.hill.0    38     9    64   38    28.2  45.2  17   13.3   37.1 11.4   1.85  3.75
# 4 mga            tot.hill.0    22    14    64   41    29    42    13   14.8   37.0 12.1   2.58  5.36
# 5 dih            tot.hill.0    40    24    81   48.5  39.5  56.5  17   13.3   49.0 14.1   2.23  4.52
# 6 positive       tot.hill.0    33    20    73   43    35    49    14   11.9   43.3 13.4   2.33  4.75


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, rest.meth.type, tot.hill.0) %>% 
  group_by(rest.meth.type) %>%
  identify_outliers(tot.hill.0)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(tot.hill.0 ~ rest.meth.type, data = data, na.action = na.omit)
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
  shapiro_test(tot.hill.0)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "tot.hill.0", facet.by = "rest.meth.type")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(tot.hill.0 ~ rest.meth.type)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(tot.hill.0 ~ rest.meth.type)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd      F        p p<.05   ges
# 1 rest.meth.type   5 180 22.745 1.17e-17     * 0.387

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(5, 189) = 22.745, p < 0.001, eta2[g] = 0.387.



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(tot.hill.0 ~ rest.meth.type)
pwc
# term           group1   group2   null.value estimate conf.low conf.high         p.adj p.adj.signif
# * <chr>          <chr>    <chr>         <dbl>    <dbl>    <dbl>     <dbl>         <dbl> <chr>       
#   1 rest.meth.type negative cus               0    6.36    -3.22      15.9  0.398         ns          
#   2 rest.meth.type negative res               0   15.3      7.07      23.4  0.00000359    ****        
#   3 rest.meth.type negative mga               0   15.1      5.63      24.5  0.000116      ***         
#   4 rest.meth.type negative dih               0   27.1     19.0       35.2  0             ****        
#   5 rest.meth.type negative positive          0   21.4     13.0       29.9  0.00000000014 ****        
#   6 rest.meth.type cus      res               0    8.89    -0.384     18.2  0.0686        ns          
#   7 rest.meth.type cus      mga               0    8.72    -1.69      19.1  0.158         ns          
#   8 rest.meth.type cus      dih               0   20.7     11.5       29.9  0.0000000121  ****        
#   9 rest.meth.type cus      positive          0   15.1      5.54      24.6  0.000138      ***         
#   10 rest.meth.type res      mga               0   -0.177   -9.32       8.96 1             ns          
#   11 rest.meth.type res      dih               0   11.8      4.09      19.5  0.00026       ***         
#   12 rest.meth.type res      positive          0    6.17    -1.95      14.3  0.248         ns          
#   13 rest.meth.type mga      dih               0   12.0      2.94      21.1  0.00254       **          
#   14 rest.meth.type mga      positive          0    6.35    -3.04      15.7  0.377         ns          
#   15 rest.meth.type dih      positive          0   -5.65   -13.7        2.38 0.331         ns      

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "rest.meth.type")
ggboxplot(data, x = "rest.meth.type", y = "tot.hill.0") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# they can be wrong positioned if the order is changed manually



# and beautiful

# analysis of variance
anova <- aov(tot.hill.0 ~ rest.meth.type, data = data)
# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey, reversed = T)

# table with factors and median
dt <- data %>% 
  group_by(rest.meth.type) %>% 
  get_summary_stats(tot.hill.0, type = "median") %>% 
  # arrange(median)
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$rest.meth.type)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "rest.meth.type", y = "tot.hill.0",
          fill = "rest.meth.type",
          legend = "none",
) +
  labs(
    x = "Site type",
    y = "Species Richness (16 m²)",
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
    aes(x = rest.meth.type, y = 80, label = cld), vjust = -0.7) +
  scale_fill_manual(
    values = c("lightgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "lightgrey"),
                    name = "Site type",
                    labels = c("Negative Reference",
                               "Cultivar Seed Mixture",
                               "Regional Seed Mixture",
                               "Management Adaptation",
                               "Direct Harvesting",
                               "Positive Reference"))

ggsave("outputs/figures/plants_species_diversity/anova_hill0-meth&type.jpg",
       dpi = 300, width = 16, height = 20, units = "cm")




