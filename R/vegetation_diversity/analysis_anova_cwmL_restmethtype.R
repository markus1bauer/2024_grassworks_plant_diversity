#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# ANOVA Species Diversity
# CWM Ellenberg L-value ~ restmeth-sitetype
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
# variable = site.cwm.abu.oek.l

## sites data
sites <- read_csv(
  here("data", "processed", "sites_processed_environment_nms_20240923.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )) %>%
  dplyr::select(
    id.site, site.type, rest.meth, site.cwm.abu.oek.l
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

data <- sites

### Visualization ###
data %>% 
  ggplot(aes(x = rest.meth.type,
             y = site.cwm.abu.oek.l)) +
  geom_boxplot()


### summary statistics ###
data %>%
  group_by(rest.meth.type) %>%
  get_summary_stats(site.cwm.abu.oek.l, type = "full") %>% 
  arrange(median)# rest.meth.type variable               n   min   max median    q1    q3   iqr   mad  mean    sd    se    ci
# <chr>          <fct>              <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
# 1 negative       site.cwm.abu.oek.l    32  6.09  7.90   6.98  6.74  7.40 0.654 0.45   7.05 0.493 0.087 0.178
# 2 dih            site.cwm.abu.oek.l    40  6.59  7.76   7.00  6.85  7.18 0.334 0.243  7.03 0.258 0.041 0.082
# 3 positive       site.cwm.abu.oek.l    33  6.61  7.91   7.08  6.88  7.30 0.422 0.315  7.10 0.292 0.051 0.103
# 4 cus            site.cwm.abu.oek.l    21  6.42  7.84   7.20  6.98  7.32 0.343 0.182  7.19 0.32  0.07  0.146
# 5 res            site.cwm.abu.oek.l    38  6.46  7.99   7.21  7.03  7.35 0.324 0.237  7.2  0.279 0.045 0.092
# 6 mga            site.cwm.abu.oek.l    22  6.54  8.07   7.34  7.08  7.71 0.635 0.517  7.36 0.441 0.094 0.196

# relevel according to type and median of rest.meth
data <- data %>% 
  mutate(rest.meth.type = fct_relevel(rest.meth.type, "negative", "cus", "res", "mga", "dih", "positive"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Assumptions ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Outliers ###
data %>% 
  select(id.site, rest.meth.type, site.cwm.abu.oek.l) %>% 
  group_by(rest.meth.type) %>%
  identify_outliers(site.cwm.abu.oek.l)
# no extreme outliers

### Check normality assumption ###
# Build the linear model
model  <- lm(site.cwm.abu.oek.l ~ rest.meth.type, data = data, na.action = na.omit)
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
  shapiro_test(site.cwm.abu.oek.l)
# Note that, if your sample size is greater than 50, the normal QQ plot is preferred
# because at larger sample sizes the Shapiro-Wilk test becomes very sensitive
# even to a minor deviation from normality:
ggqqplot(data, "site.cwm.abu.oek.l", facet.by = "rest.meth.type")

### Check homogeneity of variance assumption ###
# residuals versus fits plot
plot(model, 1)
# In the plot above, there is no evident relationships between residuals and 
# fitted values (the mean of each groups), which is good. So, we can assume the 
# homogeneity of variances.

# Levene's test #
data %>% levene_test(site.cwm.abu.oek.l ~ rest.meth.type)
# From the output above, we can see that the p-value is > 0.05, which is not 
# significant. This means that, there is not significant difference between 
# variances across groups. Therefore, we can assume the homogeneity of variances
# in the different treatment groups.


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Computation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res.aov <- data %>% 
  anova_test(site.cwm.abu.oek.l ~ rest.meth.type)
res.aov
# ANOVA Table (type II tests)
# 
# Effect DFn DFd     F     p p<.05   ges
# 1 rest.meth.type   5 180 3.316 0.007     * 0.084

# From the above ANOVA table, it can be seen that there are significant 
# differences between groups (p < 0.001), which are highlighted with “*“, 
# F(5, 189) = 22.745, p < 0.001, eta2[g] = 0.387.



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# D Post-hoc tests ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Pairwise comparisons
pwc <- data %>% 
  tukey_hsd(site.cwm.abu.oek.l ~ rest.meth.type)
pwc
# term           group1   group2   null.value estimate conf.low conf.high   p.adj p.adj.signif
# * <chr>          <chr>    <chr>         <dbl>    <dbl>    <dbl>     <dbl>   <dbl> <chr>       
# 1 rest.meth.type negative cus               0  0.142    -0.140     0.424  0.694   ns          
# 2 rest.meth.type negative res               0  0.152    -0.0887    0.393  0.455   ns          
# 3 rest.meth.type negative mga               0  0.307     0.0295    0.585  0.0207  *           
# 4 rest.meth.type negative dih               0 -0.0166   -0.255     0.221  1       ns          
# 5 rest.meth.type negative positive          0  0.0490   -0.200     0.298  0.993   ns          
# 6 rest.meth.type cus      res               0  0.00992  -0.263     0.283  1       ns          
# 7 rest.meth.type cus      mga               0  0.165    -0.141     0.471  0.63    ns          
# 8 rest.meth.type cus      dih               0 -0.159    -0.429     0.112  0.539   ns          
# 9 rest.meth.type cus      positive          0 -0.0933   -0.373     0.187  0.93    ns          
# 10 rest.meth.type res      mga               0  0.155    -0.114     0.424  0.558   ns          
# 11 rest.meth.type res      dih               0 -0.169    -0.396     0.0586 0.273   ns          
# 12 rest.meth.type res      positive          0 -0.103    -0.342     0.136  0.814   ns          
# 13 rest.meth.type mga      dih               0 -0.324    -0.590    -0.0576 0.00753 **          
# 14 rest.meth.type mga      positive          0 -0.258    -0.535     0.0178 0.0811  ns          
# 15 rest.meth.type dih      positive          0  0.0655   -0.170     0.302  0.967   ns     

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# E Visualization ############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
palette("Tableau 10")
library(multcompView)

# Visualization: box plots with p-values
# quick
pwc <- pwc %>% add_xy_position(x = "rest.meth.type")
ggboxplot(data, x = "rest.meth.type", y = "site.cwm.abu.oek.l") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) 

# !!! check if p-values are positioned correctly at the groups!!!
# they can be wrong positioned if the order is changed manually



# and beautiful

# analysis of variance
anova <- aov(site.cwm.abu.oek.l ~ rest.meth.type, data = data)
# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey, reversed = T)

# table with factors and median
dt <- data %>% 
  group_by(rest.meth.type) %>% 
  get_summary_stats(site.cwm.abu.oek.l, type = "median") %>% 
  # arrange(median)
  arrange(desc(median))

# extracting the compact letter display (cld) and adding to the dt table
cld <- as.data.frame.list(cld$rest.meth.type)
dt$cld <- cld$Letters

print(dt)

ggboxplot(data, x = "rest.meth.type", y = "site.cwm.abu.oek.l",
          fill = "rest.meth.type",
          legend = "none",
) +
  labs(
    x = "Site type",
    y = "CWM Ellenberg L value (EIV-L)",
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
    aes(x = rest.meth.type, y = 8, label = cld), vjust = -1.2) +
  scale_fill_manual(
    values = c("lightgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "lightgrey"),
    name = "Site type",
    labels = c("Negative Reference",
               "Cultivar Seed Mixture",
               "Regional Seed Mixture",
               "Management Adaptation",
               "Direct Harvesting",
               "Positive Reference"))


ggsave("outputs/figures/plants_species_diversity/anova_cwmL-meth&type.jpg",
       dpi = 300, width = 16, height = 18, units = "cm")




