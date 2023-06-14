library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(ggrepel)

# import files
human_rbp = data.table::fread("/Users/Ewann/splicing_comparison/data/human_rbps.csv",header=TRUE)




# RBP_TDP correlation
corr = as.data.table(readxl::read_excel('/Users/Ewann/splicing_comparison/data/rbp_deseq/correlations_with_tdp_level.xlsx')) #genes RNA level correlated across SH-curve with TDP level 
output_path = "/Users/Ewann/splicing_comparison/data/rbp_deseq"
 
tdp_corr = function (filepath, experiment){
  deseq = fread(filepath)
  rbp_corr_tdp = corr |> left_join(human_rbp,by =c ('variable' = 'ensmbleID')) |> 
    filter(!is.na(gene_name)) |> 
    filter(p_value < 0.01) |> 
    left_join(deseq, by = c("variable" = "ensgene")) |> 
    filter(estimate * log2FoldChange > 0) #correlation with TDP-43 level and effect in Humphrey datsets goes in the same direction (neg * neg >0;pos * pos > 0, neg * pos < 0) 
  
  fwrite(rbp_corr_tdp, paste0(output_path, "/", experiment, ".rbp_corr_tdp.csv"))
}

tdp_corr("/Users/Ewann/splicing_comparison/data/deseq2/seddighi.ipscCortical_neuron.DESEQ2_results.csv", "seddighi.ipscCortical_neuron")
tdp_corr("/Users/Ewann/splicing_comparison/data/deseq2/klim.ipscMotor_neuron.DESEQ2_results.csv", "klim.ipscMotor_neuron")
tdp_corr("/Users/Ewann/splicing_comparison/data/deseq2/brown.ipscCortical_neuron.DESEQ2_results.csv", "brown.ipscCortical_neuron")





# Correlating ELAL3 expression levels and # cyptic events
my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

log2fc = function(file){
  data = file |> 
    filter(gene_name == "ELAVL3") 
  return(data$log2fold_change)
}

padj = function(file){
  data = file |> 
    filter(gene_name == "ELAVL3") 
  return(data$padj)
}

name = function(file){
  data = file |> 
    filter(gene_name == "ELAVL3") 
  return(data$experiment)
}

experiment_table <- data.table::fread("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated.csv", header = T)

# get a list of filepaths
files = list.files("/Users/Ewann/splicing_comparison/data/deseq2",
                   pattern = "DESEQ2_results.csv",
                   full.names = TRUE)

# DZ_curve file modification
dz_curve = fread("/Users/Ewann/splicing_comparison/data/deseq2/DZcurve_1|DZcurve_control.DESEQ2_results.csv")
dz = dz_curve |> 
  mutate(gene_name = symbol, .keep = 'unused') |> 
  mutate(experiment = "DZcurve_1|DZcurve_control")

fwrite(dz, "/Users/Ewann/splicing_comparison/data/deseq2/DZcurve_1|DZcurve_control.DESEQ2_results.csv")

# clean the data
deseq2 = purrr::map(files, my_clean_reader)



# extract ELAVL3 expression and experiment name from deseq results
fc_elavl3 = purrr::map(deseq2, log2fc) |> purrr::simplify()
padj_elavl3 = purrr::map(deseq2, padj)|> purrr::simplify()
experiment = purrr::map(deseq2, name)|> purrr::simplify()

elavl3 = tibble(fc_elavl3, padj_elavl3, experiment)

# merge into the tdp experiment table
metatable = left_join(experiment_table, elavl3, by = c('name' = 'experiment'))

# write the results to a csv file
fwrite(metatable, "/Users/Ewann/splicing_comparison/data/rbp_deseq/elavl3_tdp_experiment_table.csv")

metatable = fread("/Users/Ewann/splicing_comparison/samplesheet/elavl3_tdp_experiment_table.csv", header = T)


# plot the relationship
# elavl3 & cyptic junctions
metatable |> 
  filter(!is.na(fc_elavl3)) |>  
  ggplot(aes(x = fc_elavl3, y = n_cryptic_junctions, label = name)) +
  geom_point() +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 3) +
  theme_minimal() 
  
# tdp43 & cryptic junctions
metatable |> 
  ggplot(aes(x = log2fold_change_tdp, y = n_cryptic_junctions, label = name)) +
  geom_point(aes(colour = cell.type)) +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 3) +
  theme_minimal() +
  labs(
    x = "log2FC/TARDBP",
    y = "N cryptic junctions",
    colour = "Cell Type"
  )

# elavl3 & tdp43
tst = cor.test(metatable$fc_elavl3, metatable$log2fold_change_tdp, method = "pearson") |> broom::tidy()
estim = tst |> pull(estimate)
pval = tst |> pull(p.value)

metatable |> 
  filter(!is.na(fc_elavl3)) |> 
  ggplot(aes(x = fc_elavl3, y = log2fold_change_tdp, label = name)) +
  geom_point() +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 3) +
  theme_minimal() +
  labs(
    subtitle = glue::glue("Pearson, r = {estim}, p = {pval}")
  )

# library size & cryptic junctions
metatable |> 
  ggplot(aes(x = average.library.size, y = n_cryptic_junctions, label = name)) +
  geom_point(aes(colour = cell.type)) +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 3) +
  theme_minimal()  +
  labs(
    x = "Library Sizes",
    y = "N cryptic junctions",
    colour = "Cell Type"
  )

# read length $ cryptic junctions
metatable |> 
  ggplot(aes(x = read.length, y = n_cryptic_junctions, label = name)) +
  geom_point(aes(colour = cell.type)) +
  geom_text_repel(box.padding = 0.2, max.overlaps = Inf, size = 3) +
  theme_minimal()   +
  labs(
    x = "Read Length",
    y = "N cryptic junctions",
    colour = "Cell Type"
  )

# cell type & cryptic junctions
metatable |> 
  ggplot(aes(x = cell.type, y = n_cryptic_junctions)) +
  geom_boxplot(colour = "darkgrey") +
  geom_jitter(aes(colour = cell.type)) +
  theme_minimal()   +
  labs(
    x = "Cell type",
    y = "N cryptic junctions",
    colour = "Cell Type"
  )



