source("scripts/snapcount_sra_cryptics.R")



# load in selective cryptics
early_cryptics = fread("/Users/Ewann/splicing_comparison/data/early_cryptics_to_query.csv")
# early_cryptics = early_cryptics |> select(-n_tdp_studies, -mean_coverage_tdp, -n_non_tdp_studies, -mean_coverage_non_tdp) |> distinct(snapcount_input, .keep_all = TRUE)
cryptics_inKD = fread("/Users/Ewann/splicing_comparison/data/nygc/frequently_observed_events_in_kds.csv")
cryptics_nygc = fread("/Users/Ewann/splicing_comparison/results/nygc/most_corr_genes.csv")
all_nygc_junctions = fread("/Users/Ewann/splicing_comparison/data/nygc/all_nygc_junctions.csv")

early_cryptics_snapcount = early_cryptics |> 
  dplyr::select(gene_id, snapcount_input, strand) 

nygc_snapcount = cryptics_nygc |> 
  separate(junc, c("gen", "seq", "start", "end")) |> 
  mutate(junc = glue::glue("{seq}:{start}-{end}")) |> 
  select(-gene, -gen) |> 
  left_join(all_nygc_junctions, by = c('junc' = "paste_into_igv_junction")) |> 
  mutate(start = as.numeric(start) + 1, end = as.numeric(end) -1) |> 
  mutate(snapcount_input = paste0(seq, ":", start, "-", end)) |> 
  mutate(gene_id = gene)|> 
  dplyr::select(gene_id, snapcount_input, strand)

inKD_snapcount = data.table(str_split_fixed(cryptics_inKD$paste_into_igv_junction, "[:-]",3)) |> 
  mutate(start = as.numeric(V2) + 1, end = as.numeric(V3) -1) |> 
  mutate(snapcount_input = paste0(V1, ":", start, "-", end)) |> 
  dplyr::select(snapcount_input) |> 
  cbind(cryptics_inKD) |> 
  mutate(gene_id = gene, .keep = "unused") |> 
  dplyr::select(gene_id, snapcount_input, strand) 

snapcount_metalist = rbind(early_cryptics_snapcount, nygc_snapcount, inKD_snapcount) |> distinct(snapcount_input, .keep_all = TRUE)

fwrite(snapcount_metalist, "~/splicing_comparison/data/snapcount_metalist.csv")
snapcount_metalist = fread("~/splicing_comparison/data/cryptic_queries/snapcount_metalist.csv")

# generate table for querying functions
ints = GenomicFeatures::intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene) |> makeGRangesFromDataFrame() |> unique() 
gr_metalist = snapcount_metalist |> 
  filter(type %in% c("novel_donor", 'novel_acceptor')) |> 
  separate(snapcount_input, c("seqnames", "start", "end")) |> 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

two_junc = subsetByOverlaps(ints, gr_metalist, type = "any") |> 
  annotatr::annotate_regions(annotation = gr_metalist) |> 
  as.data.table() |> 
  mutate(excl = glue::glue("{seqnames}:{start}-{end}"), 
         incl = glue::glue("{annot.seqnames}:{annot.start}-{annot.end}")) |> 
  filter(incl != excl) |> 
  dplyr::select(excl, strand, incl, annot.strand, annot.gene_id)

saveRDS(two_junc, "data/cryptic_queries/two_junc.rds")
two_junc = readRDS("data/cryptic_queries/two_junc.rds")

excl = two_junc |> dplyr::select(excl, strand) |> separate(excl, c("seqname", "start", "end"))
incl = two_junc |> dplyr::select(incl, annot.strand) |> separate(incl, c("seqname", "start", "end"))

rtracklayer::export(incl, "/Users/Ewann/splicing_comparison/data/cryptic_queries/cryptic_junction_inclusions.bed")
rtracklayer::export(excl, "/Users/Ewann/splicing_comparison/data/cryptic_queries/cryptic_junction_exclusions.bed")

fwrite(two_junc, "data/cryptic_queries/two_junc.csv")

# ----------------------------------querying function-----------------------------------------------------------
tmp = list(data.table(matrix()))
for (val in rownames(two_junc)){
  id = as.numeric(val)
  tmp[[id]] = combine_two_junctions(inc1 = two_junc[id]$incl,
                                    excl = two_junc[id]$excl,
                                    strand_code = two_junc[id]$annot.strand,
                                    data_source = 'encode1159')} 

is_null = purrr::map(tmp, function(df){is.null(dim(df))}) 
encode_metatable = tmp[which(is_null == FALSE)] |> rbindlist() |> filter(!is.na(psi)) |> unique()


encode_psi = encode_metatable |> 
  # filter(inclusion_count1 > 2) |>
  # filter(junction_coverage > 1e+6) |>
  # filter(psi > 0.05) |>
  filter(!is.na(inclusion_count1)) |> 
  mutate(excl = glue::glue("{chromosome}:{start}-{end}")) |> 
  filter(excl %in% two_junc$excl) |> 
  filter(!(excl %in% two_junc[incl == excl]$excl)) |>
  right_join((two_junc |> dplyr::select(excl, annot.gene_id))) |> 
  unique()
# --------------------------------------------------------------------------------------------------------------

saveRDS(encode_psi,"~/splicing_comparison/data/cryptic_queries/encode_psi.rds")
saveRDS(gtex_psi, "~/splicing_comparison/data/cryptic_queries/gtex_psi.rds")
sra_annotable = readRDS("~/splicing_comparison/data/cryptic_queries/srav3h_psi_annotable.rds")
encode_psi = readRDS("~/splicing_comparison/data/cryptic_queries/encode_psi.rds")
gtex_psi = readRDS("~/splicing_comparison/data/cryptic_queries/gtex_psi.rds")

# gene count by experiment
encode_psi |>
  filter(incl_1 == "chr3:124701599-124702037") |> 
  filter(!is.na(inclusion_count1)) |> 
  mutate(control = ifelse((Experiment.target %in% c("Non-specific target control-human", "")), "Yes", "No")) |> 
  mutate(condition = ifelse((Experiment.target %in% c("Non-specific target control-human", "")), "Control", gsub("-human", "", Experiment.target))) |> 
  group_by(Experiment.accession) |> 
  summarise(average_count = mean(inclusion_count1), control, condition) |> 
  ungroup() |> 
  unique() |> 
  mutate(Experiment.accession = forcats::fct_reorder(Experiment.accession, average_count)) |> 
  slice_max(average_count,n = 75) |>
  ggplot(aes(x = Experiment.accession, y = average_count, fill = control)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = condition), size = 2, nudge_y = -1.8) +
  geom_hline(yintercept = 20, colour = "turquoise", linetype = "dashed") +
  geom_hline(yintercept = 3, color = "grey80", linetype = "dashed") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(
    title = "ENCODE KALRN_chr3:124701599-124702037",
    x = "Experiment",
    y = "Count",
    fill = "Control") 
  
# SNRNP70
encode_psi |> 
  filter(grepl("TARDBP", Experiment.target)) |> 
  group_by(incl_1) |> 
  summarise(incl_count = mean(inclusion_count1), annot.gene_id) |> 
  unique() |> 
  ungroup() |> 
  mutate(junc = glue::glue("{annot.gene_id}_{incl_1}")) |> 
  mutate(junc = forcats::fct_reorder(junc, incl_count)) |>
  ggplot(aes(x = junc, y = incl_count)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "TARDBP-human",
    subtitle = "ENCODE KD experiment",
    x = "Junction",
    y = "Counts"
  )

# average count by KD condition
rbp_count_by_gene = encode_psi |> 
  mutate(condition = ifelse((Experiment.target %in% c("Non-specific target control-human", "")), "Control", gsub("-human", "", Experiment.target))) |> 
  group_by(condition, incl_1) |> 
  summarise(average_count = mean(inclusion_count1), annot.gene_id) |> 
  ungroup() |> 
  unique()

fwrite(rbp_count_by_gene, "~/splicing_comparison/results/rbpKD_count_by_gene.csv")
rbp_count_by_gene = fread("~/splicing_comparison/results/rbpKD_count_by_gene.csv")

# hclustering - RBP
rbp_count_by_gene |> 
  mutate(incl_1 = glue::glue("{annot.gene_id}_{incl_1}")) |> 
  select(incl_1,condition,average_count) |> unique() |>  
  pivot_wider(names_from = "condition",values_from = 'average_count',values_fill = 0) |> 
  tibble::column_to_rownames('incl_1') |> t() |> 
  dist() |> hclust() |> as.dendrogram() |> plot(cex = 0.4, xlim = c(254, 256))

# hclustering - cryptic junction
rbp_count_by_gene |> 
  mutate(incl_1 = glue::glue("{annot.gene_id}_{incl_1}")) |> 
  select(incl_1,condition,average_count) |> unique() |>  
  pivot_wider(names_from = "condition",values_from = 'average_count',values_fill = 0) |> 
  tibble::column_to_rownames('incl_1') |> 
  dist() |> hclust() |> plot(ylim = c(60,100))

# plotting expression level by experimental condition - KALRN
rbp_count_by_gene |> 
  filter(incl_1 == "chr3:124701256-124702037") |> 
  slice_max(average_count,n = 30) |>
  mutate(condition = forcats::fct_reorder(condition, average_count)) |> 
  ggplot(aes(x = condition, y = average_count)) + 
  geom_col() + 
  coord_flip() +
  geom_hline(yintercept = 20, color = "turquoise", linetype = "dashed") +
  geom_hline(yintercept = 3, color = "grey80", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = 'ENCODE KALRN_chr3:124701256-124702037',
    x = "Experiment condition",
    y = "Average count") 



# KALRN junctions in gtex database
# chr3:124700034-124701092
# chr3:124700034-124700976
# chr3:124701599-124702037
# chr3:124701256-124702037

gtex_psi |> 
  group_by(incl_1, Body_Site) |> 
  summarise(average_count = mean(inclusion_count1), annot.gene_id) |> 
  ungroup() |> 
  unique() |> 
  filter(incl_1 == "chr3:124701256-124702037") |> 
  mutate(Body_Site = fct_reorder(Body_Site, average_count)) |> 
  ggplot(aes(x = average_count, y = Body_Site)) + 
  geom_col() +
  labs(
    title = "GTEx expression by body site", 
    subtitle = "KALRN_chr3:124701256-124702037"
  ) +
  theme_minimal()


# query for new junctions found in TDP(longer acceptor) and SNRP70KD
new_snrnp70 = combine_two_junctions(inc1 = "chr3:124700034-124700848",
                                   excl = "chr3:124700034-124702037",
                                   strand_code = "+",
                                   data_source = 'encode1159') 

# query for new splicing events of KALRN in ENCODE
count_query(snapcount_coords = "chr3:124701599-124702037", strand_code = "+", data_source = 'encode1159')
# all junctions found in ENCODE:
# chr3:124700034-124700976 
# chr3:124700034-124702037
# chr3:124701256-124702037
# chr3:124701599-124702037

kalrn_junc |> 
  # filter(start == "124700034", end == "124700976") |> 
  filter(Experiment.accession == "ENCSR635BOO") |> 
  view()

  




# =====================================================================================================================
# calculating sum_psi for tdp and non_tdp studies
sum_psi_encode = encode_psi |> 
  group_by(Experiment.accession) |> 
  summarise(sum_psi = sum(psi), Experiment.target, Biosample.term.name) |> 
  unique() |> 
  ungroup() |> 
  mutate(tdp_studies = (grepl("TDP", Experiment.target, ignore.case = TRUE)) | 
           (grepl("TARDBP", Experiment.target, ignore.case = TRUE)))

encode_psi |> 
  filter(Experiment.accession == "ENCSR927SLP") |> 
  ggplot(aes(x = annot.gene_id, y = psi)) +
  geom_point() +
  theme_minimal() +
  theme(
    axis.text.x = element_blank()
  ) + 
  ggrepel::geom_text_repel(aes(label = annot.gene_id), max.overlaps = 5, size = 3) +
  labs(
    title = "LIN28B",
    subtitle = "Experiment: ENCSR927SLP"
  )

sum_psi_encode |> 
  # mutate(tdp_studies = (grepl("TDP", study_title, ignore.case = TRUE)) | 
  #          (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
  #          (grepl("TARDBP", study_title, ignore.case = TRUE)) |
  #          (grepl("TDP", sample_title, ignore.case = TRUE))) |>
  ggplot(aes(x = sum_psi,
             fill = tdp_studies)) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  labs(
    title = "Compilation: ENCODE"
  ) +
  theme_minimal()

# comparing sum_z with sum_psi in non_tdp studies
# pull experiment names of non-TDP studies with high sum_z
high_z = sample_sum |> 
  mutate(tdp_studies = (grepl("TDP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (grepl("TARDBP", study_title, ignore.case = TRUE)) |
           (grepl("TDP", sample_title, ignore.case = TRUE))) |> 
  filter(sum_z > 25 && tdp_studies == FALSE)  

psi_z = sum_psi |> 
  select(sum_psi, experiment_acc) |> 
  right_join(high_z, sum_psi, by = c("experiment_acc")) |> 
  replace(is.na(psi_z), 0) 

psi_z|>   
  ggplot(aes(x = sum_psi, y = sum_z)) + 
  geom_point() + 
  ggrepel::geom_text_repel(aes(label = study_title), max.overlaps = 3, size = 3) +
  theme_minimal()

# clustering psi
filtered_sample = sum_psi |> 
  slice_max(order_by = sum_psi, n = 10000) |> 
  pull(experiment_acc)

psi_cluster = psi_annotable |> 
  filter(experiment_acc %in% filtered_sample) |> 
  mutate(coords = glue::glue("{annot.gene_id}_{junc}.psi")) |> 
  pivot_wider(names_from = 'coords',
              values_from = 'psi',
              values_fill = 0) |> 
  distinct(experiment_acc, .keep_all = TRUE) |> 
  dplyr::select(study_title,sample_title,experiment_acc,contains('psi')) |> 
  mutate(tdp_studies = (grepl("TDP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (grepl("TARDBP", study_title, ignore.case = TRUE)) |
           (grepl("TDP", sample_title, ignore.case = TRUE)))

# UMAP
psi_umap = umap::umap(psi_cluster|> dplyr::select(contains("psi")) )
# saveRDS(std_umap,"~/splicing_comparison/data/std_umap.rds")
# std_umap = readRDS("~/splicing_comparison/data/std_umap.rds")

psi_umap$layout |> 
  as.data.frame() |> 
  cbind(psi_cluster) |>
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "none")
# scale_color_gradient(low = "turquoise", high = "salmon")



# ========================================================================================
# query for specific junctions
tmp = list(data.table(matrix(ncol = 8)))
for (val in rownames(snapcount_metalist)){
  id = as.numeric(val)
  tmp[[id]] = query_specific_junction_snapcount(data_source = "srav3h", 
                                                coords = snapcount_metalist[id]$snapcount_input, 
                                                strand_code = snapcount_metalist[id]$strand,
                                                gene_name = snapcount_metalist[id]$gene_id)}  


# sample-wise z-score calculation
std_metatable = query_metatable |> 
  filter(coverage > 2) |> 
  filter(star.all_mapped_reads > 1e+6) |> 
  mutate(rpm = coverage/star.all_mapped_reads*(1e+6)) |> 
  # filter(rpm > 1) |>
  mutate(plot_col = glue::glue("{gene}_{coords}")) |> 
  pivot_wider(names_from = 'plot_col',
              values_from = 'rpm',
              values_fill = 0) |> 
  unique() |> 
  mutate(across(contains('chr'), ~ scale(.), .names = "{.col}_zscore")) |>
  dplyr::select(study_title,sample_title,experiment_acc,contains('zscore')) |> 
  mutate(tdp_studies = (grepl("TDP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (grepl("TARDBP", study_title, ignore.case = TRUE)) |
           (grepl("TDP", sample_title, ignore.case = TRUE)))

# Identify non-tdp samples with high zscores
sample_sum = std_metatable |> 
  melt(id.vars = c("study_title","sample_title","experiment_acc")) |>  
  group_by(experiment_acc) |> 
  summarize(sum_z = sum(value),sample_title,study_title) |> 
  unique()

sample_sum |> 
  mutate(tdp_studies = (grepl("TDP", study_title, ignore.case = TRUE)) | 
           (grepl("TARDBP", sample_title, ignore.case = TRUE)) |
           (grepl("TARDBP", study_title, ignore.case = TRUE)) |
           (grepl("TDP", sample_title, ignore.case = TRUE))) |> 
  ggplot(aes(x = sum_z,
             fill = tdp_studies)) + 
  geom_density(alpha = 0.3) + 
  scale_x_continuous(trans = scales::pseudo_log_trans()) +
  theme_minimal()

sample = "SRX1966109"
tmp = std_metatable |> 
  filter(experiment_acc == sample) |> 
  pivot_longer(cols = contains("zscore"), names_to = "gene", values_to = "zscore", names_pattern = "(\\w+)_.*") 

tmp |> 
  ggplot(aes(x = gene, y = zscore)) +
  ggrepel::geom_text_repel(aes(label = gene), max.overlaps = 5, size = 3) +
  geom_point() +
  theme_minimal() +
  labs(
    title = unique(tmp$sample_title),
    subtitle = ""
  ) +
  theme(
    axis.text.x = element_blank()) 



# write out RDS objects
saveRDS(query_metatable,"~/splicing_comparison/data/cryptic_queries/query_metatable.rds")
saveRDS(std_metatable, "~/splicing_comparison/data/cryptic_queries/std_metatable.rds")
saveRDS(sample_sum, "~/splicing_comparison/data/cryptic_queries/sample_sum.rds")

query_metatable = readRDS("~/splicing_comparison/data/cryptic_queries/query_metatable.rds")
std_metatable = readRDS("~/splicing_comparison/data/cryptic_queries/std_metatable.rds")
sample_sum = readRDS("~/splicing_comparison/data/cryptic_queries/sample_sum.rds")


# CLUSTERING
std_cluster = std_metatable |> dplyr::select(contains("zscore")) 

# UMAP
std_umap = umap::umap(std_cluster)
# saveRDS(std_umap,"~/splicing_comparison/data/std_umap.rds")
# std_umap = readRDS("~/splicing_comparison/data/std_umap.rds")

std_umap$layout |> 
  as.data.frame() |> 
  cbind(std_metatable) |>
  left_join(sample_sum, by = ("experiment_acc")) |>
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = "none")
# scale_color_gradient(low = "turquoise", high = "salmon")

# kmeans clustering - decide the number of clusters

k <- 2:10

# silhouette scores
avg_sil = purrr::map(k, function(k){
  km = stats::kmeans(std_cluster, centers = k, nstart = 25)
  ss = cluster::silhouette(km$cluster, dist(std_cluster))
  mean(ss[, 3])
})

plot(k, type='b', avg_sil, 
     xlab='Number of clusters', 
     ylab='Average Silhouette Scores', 
     frame=FALSE) 


# clustering - actual clustering
kmeans = kmeans(std_cluster, 2)

kmeans_std = data.frame(kmean = rowMeans(std_cluster), cluster = kmeans$cluster)


# plot the results with ggplot2
pca = prcomp(std_cluster, center = TRUE)
pc_std = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])

kmeans_std_plot = cbind(kmeans_std, pc_std, std_cluster)
kmeans_std_plot$cluster <- factor(kmeans_std_plot$cluster)
kmeans_std_plot |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(
    # fill = mean,
    # shape = cluster,
    colour = cluster), 
    # colour = "grey21",
    size = 0.3) +
  labs(title = "K-Means clustering of z-scores",
       x = "PC1",
       y = 'PC2',
       colour =
         #   "average PSI",
         # shape = 
         "Cluster") +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  # scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 2)) +
  # scale_fill_gradient(
  #   high = "turquoise",
  #   low = "salmon",
  #   space = "Lab",
  #   guide = "colourbar",
  #   aesthetics = "fill") +
  # scale_shape_manual(values = c(21:25)) +
  theme_minimal()

# clustering with DBSCAN
dbscan_std = dbscan::dbscan(std_cluster, eps = 1, minPts = 50)
dbscan_std_plot = data.frame(rowMeans(std_cluster), cluster = dbscan_std$cluster)
dbscan_std_plot$cluster = factor (dbscan_std_plot$cluster)

dbscan_std_plot |> 
  cbind(pc_std) |> 
  cbind(std_metatable) |> 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = tdp_studies)) +
  labs(
    x = "PC1",
    y = 'PC2',
    colour = "TDP studies?") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.01, base = 2)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.1, base = 2)) +
  scale_color_manual(values = c("grey80", "orange")) +
  theme_minimal()

# DBSCAN + UMAP
dbscan_std_plot |> 
  cbind(std_metatable) |> 
  cbind(as.data.table(std_umap$layout)) |> 
  ggplot(aes(x = V1, y = V2, colour = cluster)) +
  geom_point() +
  theme_minimal()

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
           (grepl("TDP", sample_title, ignore.case = TRUE))) |> 
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

# ============ USELESS CODES (FOR NOW)==================

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

fwrite(snapcount_metalist, "/Users/Ewann/splicing_comparison/data/cryptic_queries/selective_cryptics_sra.csv")