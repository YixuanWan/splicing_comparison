library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(eulerr)

#read and the lowest files
dz <- read.csv('/Users/Ewann/splicing_comparison/data/majiq/curves_dz/Control-TDP43KD1_annotated_junctions.csv')
sy <- read.csv('/Users/Ewann/splicing_comparison/data/majiq/curves_sh/noDox-dox0075_annotated_junctions.csv')


dz_cryptic = dz %>% 
  select(gene_name,paste_into_igv_junction,junc_cat,gene_id,control_mean_psi,tdp43kd1_mean_psi,strand) %>% 
  filter(tdp43kd1_mean_psi > 0.1 & control_mean_psi < 0.05)  %>% 
  select(paste_into_igv_junction, strand) %>% 
  mutate(dz = TRUE) |> 
  unique()

  
  
sy_cryptic = sy %>% 
  select(gene_name,paste_into_igv_junction,junc_cat,gene_id,no_dox_mean_psi,dox0075_mean_psi,strand) %>% 
  filter(dox0075_mean_psi > 0.1 & no_dox_mean_psi < 0.05)  %>%
  select(paste_into_igv_junction, strand) %>% 
  mutate(sy5y = TRUE) |> 
  unique()
  
  
overlap = full_join(dz_cryptic, sy_cryptic)

overlap_df = overlap %>% 
  replace(is.na(.), FALSE) %>% 
  unique() %>% tibble::remove_rownames() %>% 
  tibble::column_to_rownames('paste_into_igv_junction')

overlap_df %>% eulerr::euler() %>% plot(quantities = TRUE)
