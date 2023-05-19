source("scripts/snapcount_sra_cryptics.R")



# load in selective cryptics
early_cryptics = fread("/Users/Ewann/splicing_comparison/data/early_cryptics_to_query.csv")
# early_cryptics = early_cryptics |> select(-n_tdp_studies, -mean_coverage_tdp, -n_non_tdp_studies, -mean_coverage_non_tdp) |> distinct(snapcount_input, .keep_all = TRUE)
cryptics_inKD = fread("/Users/Ewann/splicing_comparison/data/nygc/frequently_observed_events_in_kds.csv")
cryptics_nygc = fread("/Users/Ewann/splicing_comparison/results/nygc/most_corr_genes.csv")

early_cryptics_snapcount = early_cryptics |> 
  select(gene_id, snapcount_input, strand) 

nygc_snapcount = data.table(str_split_fixed(cryptics_nygc$junc, "[_:-]",4)) |> 
  mutate(start = as.numeric(V3) + 1, end = as.numeric(V4) -1) |> 
  mutate(snapcount_input = paste0(V2, ":", start, "-", end)) |> 
  select(snapcount_input) |> 
  cbind(cryptics_nygc) |> 
  mutate(gene_id = gene, strand = "-", .keep = "unused")|> 
  select(gene_id, snapcount_input, strand)

inKD_snapcount = data.table(str_split_fixed(cryptics_inKD$paste_into_igv_junction, "[:-]",3)) |> 
  mutate(start = as.numeric(V2) + 1, end = as.numeric(V3) -1) |> 
  mutate(snapcount_input = paste0(V1, ":", start, "-", end)) |> 
  select(snapcount_input) |> 
  cbind(cryptics_inKD) |> 
  mutate(gene_id = gene, .keep = "unused") |> 
  select(gene_id, snapcount_input, strand) 

snapcount_metalist = rbind(early_cryptics_snapcount, nygc_snapcount, inKD_snapcount) |> unique()


# query for specific junctions
  tmp = list(data.table(matrix(ncol = 8)))
  for (val in rownames(nygc_revert)){
    id = as.numeric(val)
    tmp[[id]] = query_specific_junction_snapcount(data_source = "srav3h", 
                                                  coords = nygc_revert[id]$snapcount_input, 
                                                  strand_code = nygc_revert[id]$strand,
                                                  gene_name = nygc_revert[id]$gene_id)}
  
  

# identify junctions not found in database - probably due to wrong strand_code provided 
is_null = purrr::map(tmp, function(df){is.null(dim(df))}) 
query_metatable = tmp[which(is_null == FALSE)] |> rbindlist()

# rerun the query with strand '-'
nygc_revert = nygc_snapcount[which(is_null == TRUE) - 24] |> mutate(strand = "-")

is_null = purrr::map(tmp, function(df){is.null(dim(df))}) 
query_nygc = tmp[which(is_null == FALSE)] |> rbindlist()

query_metatable = rbind(query_metatable, query_nygc)

# summarise the number of studies containing TDP43
tdp_studies = query_metatable |> 
  filter(coverage > 2) |> 
  filter((grepl("TDP", study_title, ignore.case = TRUE)) | 
         (grepl("TARDBP", study_title, ignore.case = TRUE)) | 
         (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
         (grepl("TDP", study_title, ignore.case = TRUE))) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_tdp_studies = n(), mean_coverage_tdp = mean(coverage))

snapcount_metalist = left_join(snapcount_metalist, tdp_studies, by = c("snapcount_input" = "coords"))

# summarise the number of studies not containing TDP43
non_tdp_studies = query_metatable |> 
  filter(coverage > 2) |> 
  filter((!grepl("TDP", study_title, ignore.case = TRUE)) | 
         (!grepl("TARDBP", sample_title, ignore.case = TRUE)) |
         (!grepl("TARDBP", study_title, ignore.case = TRUE)) |
         (!grepl("TDP", sample_title, ignore.case = TRUE))) |> 

  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(n_non_tdp_studies = n(), mean_coverage_non_tdp = mean(coverage))

snapcount_metalist = left_join(snapcount_metalist, non_tdp_studies, by = c("snapcount_input" = "coords"))


# summarise RPM for TDP and non-TDP studies
tdp_rpm = query_metatable |> 
  filter(coverage > 2) |> 
  filter((grepl("TDP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (grepl("TDP", study_title, ignore.case = TRUE))) |> 
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> 
  summarise(mean_rpm_tdp = mean(coverage/star.all_mapped_reads)*(10^6)) |> 
  mutate(z_score_tdp = (mean_rpm_tdp - median(mean_rpm_tdp))/sd(mean_rpm_tdp))
snapcount_metalist = left_join(snapcount_metalist, tdp_rpm, by = c("snapcount_input" = "coords"))

non_tdp_rpm = query_metatable |> 
  filter(coverage > 2) |> 
  filter((!grepl("TDP", study_title, ignore.case = TRUE)) | 
           (!grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (!grepl("TARDBP", study_title, ignore.case = TRUE)) |
           (!grepl("TDP", sample_title, ignore.case = TRUE))) |>
  group_by(coords) |> 
  distinct(study_title, .keep_all = TRUE) |> summarise(mean_rpm_non_tdp = mean(coverage)/mean(star.all_mapped_reads)*(10^6)) |>
  mutate(z_score_non_tdp = (mean_rpm_non_tdp - median(mean_rpm_non_tdp))/sd(mean_rpm_non_tdp))
snapcount_metalist = left_join(snapcount_metalist, non_tdp_rpm, by = c("snapcount_input" = "coords"))

# snapcount_metalist = snapcount_metalist |> select(-tdp_mean_rpm, -non_tdp_mean_rpm, -tdp_z_score, -non_tdp_z_score)


# sample-wise z-score calculation
std_metatable = query_metatable |> 
  filter(coverage > 2) |> 
  filter(star.all_mapped_reads > 1e+8) |> 
  mutate(rpm = coverage/star.all_mapped_reads*(1e+6)) |> 
  mutate(z_scores = scale(rpm, center = TRUE, scale = TRUE)) |> 
  mutate(tdp_studies = ifelse(is_tdp, TRUE, FALSE))

# create a function to determine tdp and non-tdp studies - DOESN'T WORK YET :(
# is_tdp <- function(x){
#   
#   if(grepl("TDP", df$study_title, ignore.case = TRUE)) 
#      {return(TRUE)} 
#   if(grepl("TARDBP", df$study_title, ignore.case = TRUE)) 
#     {return(TRUE)} 
#   if(grepl("TARDBP", df$sample_title, ignore.case = TRUE)) 
#     {return(TRUE)}
#   if(grepl("TDP", df$study_title, ignore.case = TRUE))
#     {return(TRUE)}
#   
#   else
#     {return(FALSE)}
# }


# comparing tdp_z_score vs non_tdp_z_score
snapcount_metalist |> 
  pivot_longer(cols = c(z_score_tdp, z_score_non_tdp), names_to = "TDP", values_to = "z_score") |> 
  ggplot(aes(x = TDP, y = z_score)) +
  geom_boxplot() +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.001)) +
  theme_minimal()

snapcount_metalist |> 
  filter(gene_id != "ZNF383") |> 
  ggplot(aes(x = tdp_z_score, y = non_tdp_z_score)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey") +
  # geom_smooth(method = "lm") +
  geom_point(na.rm = TRUE) +
  ggrepel::geom_text_repel(aes(label = gene_id), na.rm = TRUE, size = 3) + 
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.001)) +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.001)) +
  scale_color_gradient(low = "lightblue", high = "purple") +
  ggpubr::stat_cor() +
  theme_minimal()

# Examine co-occurrence of cryptic events
cryptic_map = query_metatable |> 
  filter(coverage > 2) |> 
  filter((grepl("TDP", study_title, ignore.case = TRUE)) | (grepl("TDP", sample_title, ignore.case = TRUE)) | (grepl("TARDBP", sample_title, ignore.case = TRUE))) 
  group_by(study_title, coords, gene) |> 
  summarise(mean_coverage = mean(coverage))|> 
  pivot_wider(id_cols = study_title, names_from = c(gene, coords), values_from = mean_coverage) |> 
  column_to_rownames(var = "study_title") |>
  t() |> as.matrix()

cryptic_map = cryptic_map |> replace(is.na(tmp), 0)

pheatmap(cryptic_map, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, 
         fontsize = 8, treeheight_row = 60, colorRampPalette(c("azure", "orange"))(10))

fwrite(snapcount_metalist, "/Users/Ewann/splicing_comparison/data/selective_cryptics_sra.csv")


