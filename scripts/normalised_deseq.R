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
  name = chartr(".", "_", name)
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
sedigghii3_control = get_sample_name("first_weeks/WARD_BAMS_NEW/ward_i3_newer_longer_bams_sample_sheet.csv", "CTL")
browni3_control = get_sample_name("alb_projects/data/ward_bams/ward_bams.csv", "control")

controls_id = c(humphreyi3_control, dz_control, sy5y_control, sedigghii3_control, browni3_control) 
  
# get tdpkd names
humphreyi3_tdp = get_sample_name("first_weeks/splicing_comparisons/symbams/humphrey_i3cortical/majiq_sample_sheet.csv", "tdp43KD")
dz_tdp = get_sample_name("first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/sample_sheet_dz_curves.csv", "DZ_curves_1")
sy5y_tdp = get_sample_name("first_weeks/TDP_CHX_CLONES_GLIA/sample_sheet.csv", "dox_075")
sedigghii3_tdp = get_sample_name("first_weeks/WARD_BAMS_NEW/ward_i3_newer_longer_bams_sample_sheet.csv", "TDP43_KD")
browni3_tdp = get_sample_name("alb_projects/data/ward_bams/ward_bams.csv", "tdpKD")

tdpkd_id = c(humphreyi3_tdp, dz_tdp, sy5y_tdp, sedigghii3_tdp, browni3_tdp)



metatable = purrr::map(filepath, my_clean_reader) |>   
  purrr::reduce(full_join, by = c("gene", "gene_name")) 

# select controls  
genecounts = metatable |> select(all_of(controls_id))
genenames = metatable |> select(gene_name)

# select tdpkds
genecounts_tdp = metatable |> select(all_of(tdpkd_id))


# Named vector of size factors for each sample in matrix
sfs <- DESeq2::estimateSizeFactorsForMatrix(genecounts)

# same thing for tdp
sfs_tdp <-DESeq2::estimateSizeFactorsForMatrix(genecounts_tdp)


# Divide each count column (sample) by its corresponding size factor
norm_counts <- sweep(x = genecounts,
                     MARGIN = 2, # operate on columns
                     STATS = sfs,
                     FUN = '/')

# same for tdp kd group
norm_counts_tdp <- sweep(x = genecounts_tdp,
                     MARGIN = 2, # operate on columns
                     STATS = sfs_tdp,
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
                            variable%in% sedigghii3_control ~ 'Sedigghi-i3',
                            variable%in% browni3_control ~ 'Brown-i3',
                              TRUE ~ NA_character_))  |> 
  filter(!is.na(celltype))

# same for tdp_kd
norm_counts_long_tdp = norm_counts_tdp |> 
  cbind(genenames) |> 
  as.data.table() |> 
  melt() 

norm_counts_celltype_tdp = norm_counts_long_tdp |> 
  mutate(celltype = case_when(variable %in% humphreyi3_tdp ~ 'Humphrey-i3',
                              variable %in% dz_tdp ~ 'SK-N-DZ',
                              variable%in% sy5y_tdp ~ 'SH-SY-5Y',
                              variable%in% sedigghii3_tdp ~ 'Sedigghi-i3',
                              variable%in% browni3_tdp ~ 'Brown-i3',
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
    title = "Baseline"
  ) 

# same for tdp kd
norm_counts_celltype_tdp |>
  filter(gene_name == "TARDBP") |> 
  ggplot(aes(x = celltype, y = value)) +
  geom_boxplot(colour = "darkgrey") +
  geom_jitter(aes(colour = celltype)) +
  theme_minimal() +
  labs(
    x = "Celltype",
    y = "TDP-43/Normalised gene counts",
    colour = "Celltype",
    title = "TDP KD"
  ) 
