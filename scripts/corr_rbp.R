library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gridExtra)

my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}


fc <- function(file){
  df = file |> 
    left_join(human_rbp, by = "gene_name") |> 
    filter(is.na(ensmbleID) == FALSE) 
  return(df$log2fold_change)
}

padj <- function(file){
  df = file |> 
    left_join(human_rbp, by = "gene_name") |> 
    filter(is.na(ensmbleID) == FALSE) 
  return(df$padj)
}

gene_name <- function(file){
  df = file |> 
    left_join(human_rbp, by = "gene_name") |> 
    filter(is.na(ensmbleID) == FALSE) 
  return(df$gene_name)
}

experiment <- function(file){
  df = file |> 
    left_join(human_rbp, by = "gene_name") |> 
    filter(is.na(ensmbleID) == FALSE) 
  return(df$experiment)
}

correlate_function_with_nesting = function(df){
  df |>  
    mutate(cor.test(fc , n_cryptic_junctions) |>  
             broom::tidy()) |> janitor::clean_names()
}
safe_func = purrr::possibly(correlate_function_with_nesting,otherwise = NA)

# import files
human_rbp = data.table::fread("/Users/Ewann/splicing_comparison/data/rbp_deseq/human_rbps.csv",header=TRUE)
experiment_table <- data.table::fread("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated.csv", header = T)

# get a list of filepaths
files = list.files("/Users/Ewann/splicing_comparison/data/deseq2",
                   pattern = "DESEQ2_results.csv",
                   full.names = TRUE)

# clean the files
deseq2 = purrr::map(files, my_clean_reader)


# get log2fc, padj, gene_name, experiment from the files
fc = purrr::map(deseq2, fc) |> purrr::simplify()
padj = purrr::map(deseq2, padj) |> purrr::simplify()
gene_name = purrr::map(deseq2, gene_name) |> purrr::simplify()
experiment = purrr::map(deseq2, experiment) |> purrr::simplify()

# merge into a table
rbp_table_long = tibble(experiment, gene_name, fc, padj) |> 
  left_join(experiment_table, by = c('experiment' = 'name'))

# perform correlation by rbp
rbp_table = rbp_table_long |> 
  # filter(padj < 0.1) |>
  group_by(gene_name) |> 
  nest() |> 
  mutate(correlation_result = purrr::map(data,safe_func)) |>  
  unnest(cols = correlation_result)

#write out the results
rbp_corr_significant = rbp_table |> 
  select(gene_name, p_value, estimate) |>
  filter(!is.na(p_value)) |> 
  # filter(p_value < 0.05 | gene_name == "TARDBP") |> 
  unique() |> 
  arrange(p_value)

fwrite(rbp_corr_significant, "/Users/Ewann/splicing_comparison/data/rbp_deseq/corr_rbp_cryptic.csv")

# make plots for HNRNPF, SRSF3, SYNCRIP
corr_plot <- function(RBP) {
  RBP_table = rbp_table_long |> 
    filter(gene_name == RBP) |> 
    left_join(rbp_table, by = c("gene_name", "experiment") ) |> 
    select(gene_name, experiment, fc.x, n_cryptic_junctions.x, p_value) |> 
    unique() |> 
    mutate(significance = ifelse((p_value<0.05 & !is.na(p_value)), 'Yes', 'No'), .keep = 'unused')
    
  RBP_table |>   
    ggplot(aes(x = fc.x, y = n_cryptic_junctions.x, colour = significance)) +
    geom_smooth(data = RBP_table[RBP_table$significance == "Yes",], method = 'lm', colour = "orange", linewidth = 0.7, alpha = 0.5) +
    geom_point() +
    scale_colour_manual(values = c("grey", "darkorange")) +
    labs(colour = "padj<0.1?") +
    xlab(paste0("log2FoldChange/", RBP)) +
    ylab("N cryptic junctions") +
    theme_minimal()
}

corr_plot("HNRNPF")
corr_plot("SRSF3")
corr_plot("SYNCRIP")

corr_full <- function (RBP) {
rbp_table |> 
    filter(gene_name == RBP) |> 
    select(gene_name, experiment, fc, n_cryptic_junctions, p_value) |>  
    filter(!is.na(fc)) |> 
    ggplot(aes(x = fc, y = n_cryptic_junctions)) +
    geom_smooth(method = 'lm') +
    geom_point() +
    xlab(paste0("log2FoldChange/", RBP)) +
    ylab("N cryptic junctions") +
    theme_minimal()
}


# Plotting for all RBPs
RBP = human_rbp$gene_name
pdf("/Users/Ewann/splicing_comparison/data/rbp_deseq/corr_rbp_cryptic.pdf")
for (r in RBP){
  plot = corr_full(r)
  print(plot)
}
dev.off()



# check the baseMean of HNRNPCL2
