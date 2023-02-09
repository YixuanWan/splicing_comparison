library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)


cluster_mount_point = "/Users/Ewann/cluster"



get_sample_name <- function(metadatapath, control){
  metadata_filepath = file.path(cluster_mount_point, metadatapath)
  metadata_orig = read.csv(metadata_filepath,header=TRUE)
  
  # pull out the control samples
  name = metadata_orig |> 
    filter(group == control) |> 
    pull(sample_name)|> tolower() 
  name = chartr("-", "_", name)
  
  return(name)
}

my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

# import files
filepath = list.files("/Users/Ewann/splicing_comparison/results/gene_counts",
                      pattern = ".gene_counts.csv",
                      full.names = TRUE)

# get control names
humphreyi3_control = get_sample_name("first_weeks/splicing_comparisons/symbams/humphrey_i3cortical/majiq_sample_sheet.csv", "control")
dz_control = get_sample_name("first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/sample_sheet_dz_curves.csv", "DZ_curves_0")
sy5y_control = get_sample_name("first_weeks/TDP_CHX_CLONES_GLIA/sample_sheet.csv", "no_dox")

controls_id = c(humphreyi3_control, dz_control, sy5y_control) 
  


# select controls
controls = purrr::map(filepath, my_clean_reader) |>  
  purrr::reduce(full_join, by = c("gene", "gene_name")) 
  
genecounts = controls |> select(all_of(controls_id))
genenames = controls |> select(gene_name)

# Named vector of size factors for each sample in matrix
sfs <- DESeq2::estimateSizeFactorsForMatrix(genecounts)


# Divide each count column (sample) by its corresponding size factor
norm_counts <- sweep(x = genecounts,
                     MARGIN = 2, # operate on columns
                     STATS = sfs,
                     FUN = '/')


# add a column with celltypes
norm_counts_long = norm_counts |> 
  cbind(genenames) |> 
  as.data.table() |> 
  melt() 
  
norm_counts_celltype = norm_counts_long |> 
  mutate(celltype = case_when(variable %in% humphreyi3_control ~ 'Humphrey-i3',
                            variable %in% dz_control ~ 'SK-N-DZ',
                            variable%in% sy5y_control ~ 'SH-SY-5Y',
                              TRUE ~ NA_character_))  |> 
  filter(!is.na(celltype))

# plot the results
norm_counts_celltype |>
  filter(gene_name == "TARDBP") |> 
  ggplot(aes(x = celltype, y = value)) +
  geom_boxplot(colour = "darkgrey") +
  geom_jitter(aes(colour = celltype)) +
  theme_minimal() +
  labs(
    x = "Celltype",
    y = "TDP-43/Normalised gene counts",
    colour = "Celltype",
  ) 


