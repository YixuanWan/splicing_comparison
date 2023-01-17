library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(eulerr)

  
#read and the highest files
dz_highest <- read.csv('/Users/Ewann/splicing_comparison/majiq_dz_curves/Control-TDP43KD1_annotated_junctions.csv')
sy_highest <- read.csv('/Users/Ewann/splicing_comparison/majiq_sy5y_curves/noDox-dox0075_annotated_junctions.csv')


dz_cryptic_high = dz_highest %>% 
  filter(tdp43kd1_mean_psi > 0.1 & control_mean_psi < 0.05)  %>% 
  mutate(cryptic_events = paste(paste_into_igv_junction, gene_name)) %>% 
  select(cryptic_events) %>% 
  mutate(dz_high = TRUE)


sy_cryptic_high = sy_highest %>% 
  filter(dox0075_mean_psi > 0.1 & no_dox_mean_psi < 0.05)  %>%
  mutate(cryptic_events = paste(paste_into_igv_junction, gene_name)) %>% 
  select(cryptic_events) %>%
  mutate(sy5y_high = TRUE)


overlap_highest = full_join(dz_cryptic_high, sy_cryptic_high)

#Plot the Venn diagram
overlap_df_high = overlap_highest %>% 
  replace(is.na(.), FALSE) %>%
  unique() %>% 
  tibble::remove_rownames()  %>%  
  tibble::column_to_rownames('cryptic_events')

overlap_df_high %>% euler() %>% plot(quantities = TRUE, main = "Highest TDP KD", labels = FALSE)

#Identify the cryptic events in highest TDP KD
cryptic_high = overlap_highest %>% 
  filter(dz_high == TRUE & sy5y_high == TRUE)


#read and the lowest files
dz_lowest <- read.csv('/Users/Ewann/splicing_comparison/majiq_dz_curves/Control-TDP43KD002_annotated_junctions.csv')
sy_lowest <- read.csv('/Users/Ewann/splicing_comparison/majiq_sy5y_curves/noDox-dox00125_annotated_junctions.csv')


dz_cryptic_low = dz_lowest %>% 
  filter(tdp43kd002_mean_psi > 0.1 & control_mean_psi < 0.05)  %>% 
  mutate(cryptic_events = paste(paste_into_igv_junction, gene_name)) %>% 
  select(cryptic_events) %>% 
  mutate(dz_low = TRUE)



sy_cryptic_low = sy_lowest %>% 
  filter(dox00125_mean_psi > 0.1 & no_dox_mean_psi < 0.05)  %>%
  mutate(cryptic_events = paste(paste_into_igv_junction, gene_name)) %>% 
  select(cryptic_events) %>%
  mutate(sy5y_low = TRUE)


overlap_lowest = full_join(dz_cryptic_low, sy_cryptic_low)

overlap_df_low = overlap_lowest %>% 
  replace(is.na(.), FALSE) %>% 
  unique() %>% tibble::remove_rownames() %>% 
  tibble::column_to_rownames('cryptic_events')



#plot venn diagram
overlap_df_low %>% euler() %>% plot(quantities = TRUE, main = "Lowest TDP KD", legend = list(labels = c("DZ", "SY5Y")))

#Identify shared cryptic events in lowest TDP KD
cryptic_low = overlap_lowest %>% 
  filter(dz_low == TRUE & sy5y_low == TRUE) 
  
  

#identify shared common cryptic events in both conditions
comparison = full_join(overlap_highest, overlap_lowest) %>%
  replace(is.na(.), FALSE) %>% 
  unique() %>% tibble::remove_rownames() %>% 
  #filter((dz_high == TRUE & sy5y_high == TRUE) &(dz_low == FALSE | sy5y_low == FALSE))  
  mutate(cryptic_shared = (dz_high == TRUE & sy5y_high == TRUE) &(dz_low == TRUE & sy5y_low == TRUE)) %>% 
  select(cryptic_events, cryptic_shared)

overlap_df = full_join(overlap_highest, comparison) %>% 
  replace(is.na(.), FALSE) %>% 
  unique() %>% tibble::remove_rownames() %>% 
  tibble::column_to_rownames('cryptic_events') 

overlap_df %>% euler() %>% plot(quantities = TRUE, main = "Shared Common Cryptic Events", legend = list(labels = c("DZ", "SY5Y", "Shared")))

#Identiy shared cryptic events across TDP KD levels
cryptic_shared = comparison %>% 
  filter(cryptic_shared == TRUE)  
  



