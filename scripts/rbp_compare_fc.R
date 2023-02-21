library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

make_plots <- function(i){
  deseq2_rbp |> 
    filter(experiment == i) |> 
    ggplot(aes(x = log2fold_change, colour = RBP, fill = RBP)) +
    geom_area(stat = "bin", alpha = 0.5) +
    # geom_density(alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.001, base = 2)) +
    labs(
      x = paste0("log2FoldChange/", toupper(i)),
      y = "Counts"
    ) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 10)
    )
}

# import deseq files, human_rbp files, and experiment table
files = list.files("/Users/Ewann/splicing_comparison/data/deseq2",
                   pattern = "DESEQ2_results.csv",
                   full.names = TRUE)

human_rbp = data.table::fread("/Users/Ewann/splicing_comparison/data/human_rbps.csv",header=TRUE)

metatable = data.table::fread("/Users/Ewann/splicing_comparison/samplesheet/tdp_experiments_updated.csv", header=TRUE)

# identify RBPs
deseq2 = purrr::map(files, my_clean_reader)  |> data.table::rbindlist(fill = TRUE) 

deseq2_rbp = deseq2 |> 
  filter(!is.na(log2fold_change)) |> 
  select(-geneid) |> 
  mutate(RBP = ifelse((gene_name %in% human_rbp$gene_name), "Yes", "No"))  

# identify cell types
metadata = metatable |> 
  select(name, cell.type)

deseq2_rbp = deseq2_rbp |> 
  left_join(metadata, by = c("experiment" = "name"))

# make boxplots to plot RBPs and other genes
celltype = "iPSC MN"
deseq2_rbp|> 
  filter(cell.type == celltype) |> 
  ggplot(aes(x = experiment, y = log2fold_change, fill = RBP)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  ggpubr::stat_compare_means(method = "wilcox.test") +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  labs(
    x = "Experiment",
    y = "Log2FoldChange",
    subtitle = celltype
  )

# make frequency plots with RBPs and other genes
EXPERIMENT = deseq2_rbp$experiment |> unique()
pdf(width = 7, height = 5, "/Users/Ewann/splicing_comparison/results/rbp_deseq/rbp_frequency_plots.pdf")
for (i in EXPERIMENT){
  plot = make_plots(i)
  print(plot)
}
dev.off()

  
         