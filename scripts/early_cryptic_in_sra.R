source("scripts/snapcount_sra_cryptics.R")


early_cryptics = fread("/Users/Ewann/splicing_comparison/data/early_cryptics_to_query.csv")
early_cryptics = early_cryptics |> select(-n_tdp_studies, -mean_coverage_tdp, -n_non_tdp_studies, -mean_coverage_non_tdp)

early_cryptics = early_cryptics |> 
  mutate(count = 0)

early_cryptics_query = list(data.table(matrix(ncol = 8)))

for (val in rownames(early_cryptics)){
  id = as.numeric(val)
  early_cryptics_query[[id]] = query_specific_junction_snapcount(data_source = "srav3h", 
                                                              coords = early_cryptics[id]$snapcount_input, 
                                                              strand_code = early_cryptics[id]$strand,
                                                              gene_name = early_cryptics[id]$gene_id) 
}
early_cryptics_query_long = rbindlist(early_cryptics_query)

# summarise the number of studies containing TDP43
tdp_studies = early_cryptics_query_long |> 
  filter(coverage > 2) |> 
  filter(grepl("TDP", study_title, ignore.case = TRUE)) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_tdp_studies = n(), mean_coverage_tdp = mean(coverage))

tdp_samples = early_cryptics_query_long |> 
  filter(coverage > 2) |> 
  filter(!grepl("TDP", study_title, ignore.case = TRUE)) |> 
  filter(grepl("TDP", sample_title, ignore.case = TRUE)) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_tdp_studies = n(), mean_coverage_tdp = mean(coverage))

tardbp_samples = early_cryptics_query_long |> 
  filter(coverage > 2) |> 
  filter(!grepl("TDP", study_title, ignore.case = TRUE)) |> 
  filter(grepl("TARDBP", sample_title, ignore.case = TRUE)) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_tdp_studies = n(), mean_coverage_tdp = mean(coverage))

tdp_studies = rbind(tdp_studies, tdp_samples, tardbp_samples)

early_cryptics = left_join(early_cryptics, tdp_studies, by = c("snapcount_input" = "coords"))

# summarise the number of studies not containing TDP43
non_tdp_studies = early_cryptics_query_long |> 
  filter(coverage > 2) |> 
  filter(!grepl("TDP", study_title, ignore.case = TRUE)) |> 
  filter(!grepl("TARDBP", sample_title, ignore.case = TRUE)) |> 
  filter(!grepl("TDP", sample_title, ignore.case = TRUE)) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_non_tdp_studies = n(), mean_coverage_non_tdp = mean(coverage))
early_cryptics = left_join(early_cryptics, non_tdp_studies, by = c("snapcount_input" = "coords"))

# Examine co-occurrence of cryptic events
tmp = early_cryptics_query_long |> 
  filter(coverage > 2) |> 
  filter(grepl("TDP", study_title, ignore.case = TRUE)) |>
  group_by(study_title, coords, gene) |> 
  summarise(mean_coverage = mean(coverage))|> 
  pivot_wider(id_cols = study_title, names_from = c(gene, coords), values_from = mean_coverage) |> 
  column_to_rownames(var = "study_title") |>
  t() |> as.matrix()

tmp = tmp |> replace(is.na(tmp), 0)

pheatmap(tmp, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, 
         fontsize = 8, treeheight_row = 60, colorRampPalette(c("azure", "orange"))(10))

fwrite(early_cryptics, "/Users/Ewann/splicing_comparison/data/early_cryptics_to_query.csv")


