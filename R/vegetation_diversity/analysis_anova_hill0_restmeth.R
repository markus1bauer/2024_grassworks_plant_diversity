#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# Species Richness ~ restmeth
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
# groups = rest.meth
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
  ggplot(aes(x = rest.meth,
             y = tot.hill.0)) +
  geom_boxplot()

# remove NA
data <- data %>% 
  filter(!is.na(rest.meth)) %>% 
  mutate(rest.meth = fct_relevel(rest.meth, "cus", "res", "mga", "dih"))

### summary statistics ###
data %>%
  group_by(rest.meth) %>%
  get_summary_stats(tot.hill.0, type = "full")
# rest.meth variable       n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <fct>     <fct>      <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#   1 cus       tot.hill.0    21    13    45   29    24    33       9  5.93  28.2  8.84  1.93  4.02
# 2 dih       tot.hill.0    40    24    81   48.5  39.5  56.5    17 13.3   49.0 14.1   2.23  4.52
# 3 mga       tot.hill.0    22    14    64   41    29    42      13 14.8   37.0 12.1   2.58  5.36
# 4 res       tot.hill.0    38     9    64   38    28.2  45.2    17 13.3   37.1 11.4   1.85  3.75


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, rest.meth, tot.hill.0) %>% 
  group_by(rest.meth) %>%
  identify_outliers(tot.hill.0)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(tot.hill.0 ~ rest.meth, data = data, na.action = na.omit)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, as all the points fall approximately along the reference line,
# we can assume normality. This conclusion is supported by the Shapiro-Wilk test.
# The p-value is not significant (p = 0.581), so we can assume normality.

# Check normality assumption by groups
data %>%
  filter(!is.na(rest.meth)) %>% 
  group_by(rest.meth) %>%
  shapiro_test(tot.hill.0)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "tot.hill.0", facet.by = "rest.meth")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(tot.hill.0 ~ rest.meth)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(tot.hill.0 ~ rest.meth)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd      F        p p<.05   ges
# 1 rest.meth   3 117 14.918 2.75e-08     * 0.277

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(3, 117) = 14.918, p < 0.001, eta2[g] = 0.277.



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(tot.hill.0 ~ rest.meth)
pwc
#     term      group1 group2 null.value estimate conf.low conf.high        p.adj p.adj.signif
# * <chr>     <chr>  <chr>       <dbl>    <dbl>    <dbl>     <dbl>        <dbl> <chr>       
#   1 rest.meth cus    dih             0   20.7     12.2       29.2  0.0000000276 ****        
#   2 rest.meth cus    mga             0    8.72    -0.933     18.4  0.0919       ns          
#   3 rest.meth cus    res             0    8.89     0.293     17.5  0.0398       *           
#   4 rest.meth dih    mga             0  -12.0    -20.4       -3.60 0.00171      **          
#   5 rest.meth dih    res             0  -11.8    -19.0       -4.65 0.000207     ***         
#   6 rest.meth mga    res             0    0.177   -8.30       8.65 1            ns    

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "rest.meth")
ggboxplot(data, x = "rest.meth", y = "tot.hill.0") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# I don't know why they are sometimes wrong..



# and beautiful

# analysis of variance
anova <- aov(tot.hill.0 ~ rest.meth, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

# table with factors and median
dt <- data %>% 
  group_by(rest.meth) %>% 
  get_summary_stats(tot.hill.0, type = "median") %>% 
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$rest.meth)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "rest.meth", y = "tot.hill.0",
          fill = "rest.meth",
          legend = "right",
) +
  geom_text(
    data = dt,
    aes(x = rest.meth, y = 80, label = cld), vjust = -0.5) +
  labs(
    x = "Restoration Method",
    y = "Species Richness",
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc),
    ) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = palette(),
                    name = "Restoration Method",
                    labels = c("Cultivar Seed Mixture",
                               "Regional Seed Mixture",
                               "Management Adaptation",
                               "Direct Harvesting"
                               ))







