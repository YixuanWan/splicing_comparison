library(dplyr)
library(tidyr)

my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

tdp_deseq <- function(file){
    if("symbol" %in% colnames(file)){
      df = file |> 
        filter(symbol == "TARDBP") |> 
        dplyr::select(log2fold_change, padj, deseq2_table_name)
    }
    else{
      df = file |> 
        filter(gene_name == "TARDBP") |> 
        dplyr::select(log2fold_change, padj, deseq2_table_name) 
    }
    return(df)
}

# loading the tables
splicingTable = read.csv("/Users/Ewann/splicing_comparison/data/majiq/splicing_full_delta_psi_tables.csv")
metatable = read.csv("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated-tdp_experiments_updated.csv")

# work out the number of cryptic events
crypticTable = splicingTable |> 
  filter(comparison %in% metatable$comparison_majiq) |> 
  select(gene_name,paste_into_igv_junction,junc_cat,baseline_PSI,contrast_PSI,comparison,strand) |>  
  mutate(is_cryptic = contrast_PSI > 0.1 & baseline_PSI < 0.05)  |> 
  group_by(comparison,is_cryptic) |>   
  summarize(n_jun = n_distinct(paste_into_igv_junction)) |> 
  ungroup() |>  
  filter(is_cryptic == TRUE) |>  
  right_join(metatable) |> 
  mutate(n_cryptic_junctions = n_jun, .keep = "unused") |> 
  select(-is_cryptic) |>  
  arrange(experiment)

metadata = crypticTable

# add deseq names to the table
filepath = list.files("/Users/Ewann/splicing_comparison/data/deseq2", full.names = TRUE)
suffix = ".csv"
filenames = list.files("/Users/Ewann/splicing_comparison/data/deseq2")
experiment_names = gsub(suffix,"",filenames)
estimate_files = purrr::map(filepath,my_clean_reader)
fc_tdp = purrr::map2(estimate_files, experiment_names, ~cbind(.x, deseq2_table_name = .y))

# add tdp43 deseq2 results to the metadata
fc_tdp = purrr::map(fc_tdp, tdp_deseq)  
fc_tdp = data.table::rbindlist(fc_tdp)
metatable = metatable |> 
  select(-padj, -log2fold_change_tdp) |> 
  mutate(deseq2_table_name = fc_tdp$deseq2_table_name) |> 
  left_join(fc_tdp)


# write out the table
write.csv(metadata, "/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated - tdp_experiments_updated.csv", row.names = FALSE)

# plot the number of cryptic events for selective experiments |> 
metatable |> 
  filter(comparison %in% c("controlbrowncorticalneuron-tdp43kdbrowncorticalneuron", "nodox-dox0075", "control-tdp43kd1")) |> 
  ggplot(aes(x = comparison, y = n_cryptic_junctions, fill = comparison)) +
  geom_col(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    x = "Experiment",
    y = "Number of Cryptic Junctions"
  ) +
  scale_x_discrete(labels = c("SK-N-DZ", "i3Neuron", "SH-SY5Y")) +
  scale_fill_manual(values = c("#E76BF3", "salmon", "#619CFF"))
