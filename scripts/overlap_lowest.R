library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(eulerr)

#read and the lowest files
dz_lowest <- read.csv('/Users/Ewann/splicing_comparison/majiq_dz_curves/Control-TDP43KD002_annotated_junctions.csv')
sy_lowest <- read.csv('/Users/Ewann/splicing_comparison/majiq_sy5y_curves/noDox-dox00125_annotated_junctions.csv')


dz_cryptic = dz_lowest %>% 
  select(gene_name,paste_into_igv_junction,junc_cat,gene_id,control_mean_psi,tdp43kd002_mean_psi,strand) %>% 
  filter(tdp43kd002_mean_psi > 0.1 & control_mean_psi < 0.05)  %>% 
  select(paste_into_igv_junction) %>% 
  mutate(dz = TRUE)

  
  
sy_cryptic = sy_lowest %>% 
  select(gene_name,paste_into_igv_junction,junc_cat,gene_id,no_dox_mean_psi,dox00125_mean_psi,strand) %>% 
  filter(dox00125_mean_psi > 0.1 & no_dox_mean_psi < 0.05)  %>%
  select(paste_into_igv_junction) %>% 
  mutate(sy5y = TRUE)
  
  
overlap_lowest = full_join(dz_cryptic, sy_cryptic)

overlap_df = overlap_lowest %>% 
  replace(is.na(.), FALSE) %>% 
  unique() %>% tibble::remove_rownames() %>% 
  tibble::column_to_rownames('paste_into_igv_junction')

overlap_df %>% euler() %>% plot(quantities = TRUE)
