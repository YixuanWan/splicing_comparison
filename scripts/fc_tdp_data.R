library(dplyr)
library(tidyr)
my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}

experiment_table <- read.csv("TDP-43 KD experiments - Sheet1.csv")

#REPLACE WITH YOUR DOWNLOAD FILE
file_path = file.path(here::here(),'deseq2_tables_firstpanel')
suffix = "DESEQ2_results.csv"
estimate_files = list.files(file_path,
                            pattern = suffix,
                            full.names = TRUE)

fc_tdp = purrr::map(estimate_files,my_clean_reader)
fc_tdp = data.table::rbindlist(fc_tdp)


# tdp = fc_tdp %>%
#     filter(gene_name == "TARDBP") %>%
#     select(log2fold_change, gene_name, experiment) %>% 
#     rename(name = experiment)

tdp = fc_tdp %>%
    filter(gene_name == "TARDBP") %>%
    select(log2fold_change, padj, experiment) 


updated_table = experiment_table %>% 
  left_join(tdp, by = c("name" = "experiment")) %>% 
  mutate(TARDBP = NULL, log2FoldChange.TARDBP = NULL, log2fold_change_tdp = log2fold_change, .keep = "unused") %>% 
  as.data.table() %>%
  arrange(log2fold_change_tdp)


write.csv(updated_table, 'tdp_experiments_updated.csv')


