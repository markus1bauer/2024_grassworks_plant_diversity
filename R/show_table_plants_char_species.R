#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GRASSWORKS Project
# Output table characteristic grassland species
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# author: Christin Juno Laschke




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(gt)


### Start ###
rm(list = ls())


### Load data ###
# used for diversity calculation of vegetation plots (A1-A4)
# in skript _prepare_data_3b_site_environment.R (aggregated species level)
abundances <- read_csv(
  here("data", "processed", "data_processed_abundances_A1A4_20250306.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types = cols(
    .default = "?"
  )
) 

species <- abundances %>% 
  distinct(name.plant.agg, .keep_all = T) %>% 
  select(name.plant.agg, target.species.agg, r.all.diagnostic, ffh.lrt, calth.btt) %>% 
  arrange(name.plant.agg) %>% 
  mutate(name.plant.agg = str_replace(name.plant.agg, "Gruppe", "group"),
         # add "spec." to genus taxa
         name.plant.agg = if_else(
           !str_detect(name.plant.agg, " ") & !str_detect(name.plant.agg, "group"),
           paste0(name.plant.agg, " spec."),
           name.plant.agg
         ))




species %>% 
  filter(target.species.agg == 1)
# 335 characteristic species



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(table <- species %>%
   gt() %>%
   
   ### 1 General settings ####
 opt_table_lines(
   extent = "none"
 ) %>%
   tab_options(
     table.font.style = "Arial",
     table.font.size = px(12),
     table.font.color = "black",
     data_row.padding = px(4),
     table.align = "left",
     column_labels.border.top.style = "solid",
     table_body.border.bottom.style = "solid",
     table_body.border.top.style = "solid",
     column_labels.border.top.color = "black",
     table_body.border.bottom.color = "black",
     table_body.border.top.color = "black",
     column_labels.border.top.width = px(2),
     table_body.border.bottom.width = px(2),
     table_body.border.top.width = px(1)
   )
 %>%
   
   ### 2 Bold characteristic species ####
   tab_style(
     style = list(cell_text(weight = "bold")),  # Apply bold text
     locations = cells_body(
       columns = name.plant.agg,  # Column to style
       rows = target.species.agg == 1  # Condition for bolding
     )
   ) %>% 
   cols_hide(columns = target.species.agg)
 %>% 
   
   ### 3 Species names in italic ####
 tab_style(
   style = cell_text(style = "italic"),
   locations = cells_body(
     columns = name.plant.agg,  # Column to style
     )
 ) %>% 
   cols_hide(columns = target.species.agg)
 %>% 
   
   
   ### 4 Column labels ####
 cols_label(
   name.plant.agg = md("**Plant**"),
   # target.species.agg = md("**Characteristic Grassland Species**"),
   r.all.diagnostic = md("**EUNIS-ESy**"),
   ffh.lrt = md("**Habitats Directive**"),
   calth.btt = md("**Calthion**")
 )
 )
    


### Save ###
gtsave(table, here("outputs", "tables", "table_char_spec_list_20250407.docx"))
