library(GenomicFeatures)
library(annotatr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(GenomicRanges)
library(ggplot2)

source("scripts/snapcount_sra_cryptics.R")

# summarise type of splicing
tdpkd = fread('~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/majiq/delta_psi_voila_tsv/Control-TDP43KD_annotated_junctions.csv')
srsf3kd = fread("~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/majiq/delta_psi_voila_tsv/Control-SRSF3KD_annotated_junctions.csv")
tdpsrsf3kd = fread("~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/majiq/delta_psi_voila_tsv/Control-TDP43KDSRSF3KD_annotated_junctions.csv")

tdp_cj = tdpkd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,tdp43kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & tdp43kd_mean_psi > 0.10) |> 
  arrange(-tdp43kd_mean_psi) |> 
  filter(junc_cat != 'annotated') |> 
  pull(paste_into_igv_junction) |> unique()

shared_cj = srsf3kd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,srsf3kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & srsf3kd_mean_psi > 0.10) |> 
  arrange(-srsf3kd_mean_psi) |> 
  filter(paste_into_igv_junction %in% tdp_cj) |> 
  distinct(paste_into_igv_junction, .keep_all = TRUE) |> 
  filter(junc_cat != "none")
# |> filter(junc_cat %in% c("novel_donor","novel_acceptor")) 

# Type of cryptic junctions
shared_cj |> 
  ggplot(aes(x = "", y = srsf3kd_mean_psi, fill = junc_cat)) +
  geom_bar(stat = "identity") +
  coord_polar("y", start = 0) + 
  theme_void() + 
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Shared cryptic junctions in TDP43KD and SRSF3KD",
    fill = "Junction type")

# PSI comparison
psi_cj = tdpkd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,tdp43kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & tdp43kd_mean_psi > 0.10) |> 
  arrange(-tdp43kd_mean_psi) |> 
  filter(junc_cat != 'annotated') |> 
  distinct(paste_into_igv_junction, junc_cat, .keep_all = TRUE) |> 
  select(-probability_changing, -control_mean_psi, -strand) |> 
  right_join(shared_cj, by = c("paste_into_igv_junction", "junc_cat", "gene_name")) |> 
  filter(!is.na(tdp43kd_mean_psi)) |> 
  filter(junc_cat != "none") |> 
  select(-probability_changing, -control_mean_psi, -strand) 

psi_cj |> 
  pivot_longer(cols = ends_with("mean_psi"), names_to = "condition", values_to = "mean_psi") |> 
  ggplot(aes(x = condition, y = mean_psi, fill = condition)) + 
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "PSI levels - Shared cryptic junctions in TDP43KD and SRSF3KD",
    fill = "Junction type",
    x = "Condition",
    y = "PSI"
  ) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("SRSF3KD", "TDP43KD")) +
  scale_fill_brewer(palette = "Set3")

psi_cj |> 
  mutate(psi_com = srsf3kd_mean_psi - tdp43kd_mean_psi) |> 
  ggplot(aes(x = junc_cat, y = psi_com)) +
  geom_boxplot(fill = "lightgrey") +
  theme_minimal() +
  geom_jitter(size = 0.8, aes(colour = junc_cat)) +
  theme(legend.position = "none") +
  labs(
    x = "Junction type",
    y = "Mean PSI/SRSF3KD - TDP43KD") +
  geom_hline(yintercept = 0) +
  scale_colour_brewer(palette = "Set3")

psi_cj |> 
  mutate(psi_com = srsf3kd_mean_psi - tdp43kd_mean_psi) |> 
  mutate(outlier = ifelse(((psi_com < 0.2) & (psi_com > -0.2)), FALSE, TRUE)) |> 
  ggplot(aes(x = paste_into_igv_junction, y = psi_com)) +
  geom_point(aes(colour = outlier)) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  ggrepel::geom_text_repel(aes(label = gene_name), size = 3, max.overlaps = 10) +
  labs(
    x = "Cyptic junctions",
    y = "Mean PSI/SRSF3KD - TDP43KD") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0.2, colour = 'grey') +
  geom_hline(yintercept = -0.2, colour = 'grey') +
  scale_colour_manual(values = c('darkgrey', 'orange'))

psi_cj |> 
  mutate(psi_com = srsf3kd_mean_psi - tdp43kd_mean_psi) |> 
  mutate(outlier = ifelse(((psi_com < 0.2) & (psi_com > -0.2)), FALSE, TRUE)) |> 
  ggplot(aes(x = tdp43kd_mean_psi, y = srsf3kd_mean_psi)) +
  geom_point(aes(colour = outlier)) +
  geom_abline(intercept = 0.2, colour = 'darkgrey') +
  geom_abline(intercept = -0.2, colour = 'darkgrey') +
  geom_abline() +
  theme_minimal() +
  scale_colour_manual(values = c('grey', 'orange')) +
  ggrepel::geom_text_repel(aes(label = gene_name), size = 3, max.overlaps = 5) +
  theme(legend.position = "none") +
  labs(
    x = "PSI/TDP43KD",
    y = "PSI/SRSF3KD"
  )





tdpkd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,tdp43kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & tdp43kd_mean_psi > 0.10) |> 
  arrange(-tdp43kd_mean_psi) |> filter(gene_name == "HDGFL2")

# check overlap between nygc-tdp-cortex dataset and srsf3kd
nygc_cj = nygc |>
  filter(tdp_path == "path") |>
  filter(disease_tissue == TRUE) |> 
  filter(grepl("Cortex", tissue_clean)) |> 
  pull(paste_into_igv_junction) |> 
  unique()

nygc_srsf3kd = srsf3kd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,srsf3kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & srsf3kd_mean_psi > 0.10) |> 
  arrange(-srsf3kd_mean_psi) |> 
  filter(paste_into_igv_junction %in% nygc_cj) |> 
  filter(junc_cat != "annotated") |> 
  distinct(paste_into_igv_junction, .keep_all = TRUE) |> 
  filter(junc_cat != "none")

counts = fread("~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/sj_tabs_parsed/srsf3kd_counts.aggregated.clean.annotated.bed")
srsf3_counts = counts |> 
  mutate(seq = paste0(V1, ":", V2, "-", V3), .keep= "unused") |> 
  filter(grepl("^SRSF3KD", V4)) |> 
  filter(!grepl("annotated", V7)) |> 
  filter(!grepl("none", V7)) |> 
  filter(seq %in% nygc_cj) 

mcol_srsf3 = data.frame(stringr::str_split_fixed(srsf3_counts$V7, "[|]", 3)) 
names(mcol_srsf3) = c("gene", "junc_cat", "n_tdpkd_dataset")
srsf3_counts = cbind(srsf3_counts, mcol_srsf3)

srsf3_counts = srsf3_counts |> 
  select(-V7) |> 
  group_by(seq) |> 
  mutate(mean_counts = mean(V5)) |> 
  ungroup() |> 
  distinct(seq, .keep_all = TRUE) |> 
  mutate(pass_psi = ifelse((seq %in% nygc_srsf3kd$paste_into_igv_junction), TRUE, FALSE))

srsf3kd_psi = srsf3kd |> 
  select(paste_into_igv_junction, srsf3kd_mean_psi, control_mean_psi)

srsf3_counts_psi = srsf3_counts |> 
  left_join(srsf3kd_psi, by = c("seq" = "paste_into_igv_junction")) |> 
  filter(!is.na(srsf3kd_mean_psi)) |> 
  mutate(delta_psi = log2(srsf3kd_mean_psi/control_mean_psi)) |> 
  select(-V4, -V5, -V6) |> 
  unique()
  

srsf3_counts_psi |> 
  ggplot(aes(x = delta_psi, y = mean_counts)) +
  geom_point(aes(colour = junc_cat, shape = pass_psi)) +
  theme_minimal() +
  ggrepel::geom_text_repel(aes(label = gene), size = 3, max.overlaps = 10) +
  labs(
    x = "Mean PSI_Log2(SRSF3KD/Control)",
    y = "Counts",
    colour = "Pass PSI test?"
  ) +
  # theme(
  #   axis.text.x = element_blank()
  # ) +
  # scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.5)) 
  # scale_colour_manual(values = c('black', 'red')) 
  

# Are counts correlated with mean_psi?
srsf3_sg1 = fread("~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/majiq/psi_voila_tsv_single/SRSF3KD_rep1.Aligned.sorted.out_parsed.csv")
srsf3_sg2 = fread("~/cluster/first_weeks/splicing_comparisons/tdp_srsf3/majiq/psi_voila_tsv_single/SRSF3KD_rep2.Aligned.sorted.out_parsed.csv")

srsf3_sg1_psi = srsf3_sg1 |> 
  select(gene_name, mean_psi_per_lsv_junction, paste_into_igv_junction) |> 
  filter(paste_into_igv_junction %in% srsf3_counts_psi$seq) |> 
  mutate(mean_psi_rep1 = mean_psi_per_lsv_junction, .keep = "unused")

srsf3_sg2_psi = srsf3_sg2 |> 
  select(gene_name, mean_psi_per_lsv_junction, paste_into_igv_junction) |> 
  filter(paste_into_igv_junction %in% srsf3_counts_psi$seq) |> 
  mutate(mean_psi_rep2 = mean_psi_per_lsv_junction, .keep = "unused")

srsf3_sg_psi = srsf3_counts |> 
  left_join(srsf3_sg1_psi, by = c("seq" = "paste_into_igv_junction", "gene" ="gene_name"), relationship = "many-to-many") |> 
  left_join(srsf3_sg2_psi, by = c("seq" = "paste_into_igv_junction", "gene" ="gene_name"), relationship = "many-to-many") |> 
  mutate(mean_single_psi = (mean_psi_rep1 + mean_psi_rep2)/2) |> 
  select(-V4, -V5, -V6)  |> 
  group_by(seq, junc_cat) |> 
  mutate(mean_sg_psi = mean(mean_single_psi), .keep = "unused") |> 
  ungroup() |> 
  distinct(seq, junc_cat, .keep_all = TRUE) |> 
  replace(is.na(srsf3_sg_psi), 0)

srsf3_sg_psi |> 
  ggplot(aes(y = mean_counts, x = mean_sg_psi, colour = pass_psi)) +
  geom_point(aes(shape = ifelse((mean_sg_psi > 0), "Yes", "No"))) +
  theme_minimal() +
  ggrepel::geom_text_repel(aes(label = gene), size = 3, max.overlaps = 5) +
  scale_colour_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5)) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.5)) +
  labs(
    y = "Counts",
    x = "Mean single PSI",
    colour = "Pass PSI test?",
    shape = "Have single PSI?"
  ) + 
  scale_shape_manual(values = c(4, 16))


# plot numbers of cryptic junctions in TDP43KD, SRSF3KD, and double KD groups
tdpkd_cj = tdpkd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,tdp43kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & tdp43kd_mean_psi > 0.10) |> 
  arrange(-tdp43kd_mean_psi) |> 
  filter(junc_cat != 'annotated') |> 
  select(gene_name, paste_into_igv_junction, junc_cat) |> 
  mutate(tdpKD = TRUE) |> 
  distinct(paste_into_igv_junction, junc_cat, .keep_all = TRUE)

srsf3kd_cj = srsf3kd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,srsf3kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & srsf3kd_mean_psi > 0.10) |> 
  arrange(-srsf3kd_mean_psi) |> 
  filter(junc_cat != 'annotated') |> 
  select(gene_name, paste_into_igv_junction, junc_cat) |> 
  mutate(srsf3KD = TRUE) |> 
  distinct(paste_into_igv_junction, junc_cat, .keep_all = TRUE)


tdpsrsf3kd_cj = tdpsrsf3kd |> 
  dplyr::select(gene_name,paste_into_igv_junction,junc_cat,probability_changing,control_mean_psi,tdp43kdsrsf3kd_mean_psi,strand) |>  
  filter(control_mean_psi < 0.05 & tdp43kdsrsf3kd_mean_psi > 0.10) |> 
  arrange(-tdp43kdsrsf3kd_mean_psi) |> 
  filter(junc_cat != 'annotated') |> 
  select(gene_name, paste_into_igv_junction, junc_cat) |> 
  mutate(doubleKD = TRUE) |> 
  distinct(paste_into_igv_junction, junc_cat, .keep_all = TRUE)

nygc_srsf3kd_cj = nygc_srsf3kd |> 
  select(paste_into_igv_junction, junc_cat) |> 
  mutate(nygc_srsf3kd = TRUE)
 
overlap_cj = full_join(tdpkd_cj, srsf3kd_cj, by = c("paste_into_igv_junction", "junc_cat"), relationship = "many-to-many") |> 
  full_join(tdpsrsf3kd_cj, by = c("paste_into_igv_junction", "junc_cat"), relationship = "many-to-many") |> 
  full_join(nygc_srsf3kd_cj, by = c("paste_into_igv_junction", "junc_cat"), relationship = "many-to-many") |> 
  mutate(gene_id = ifelse(is.na(gene_name.x), gene_name.y, gene_name.x), .keep = "unused") |> 
  mutate(gene_id = ifelse(is.na(gene_id), gene_name, gene_id), .keep = "unused") |> 
  filter(junc_cat != 'none') |> 
  replace(is.na(overlap_cj), FALSE) 


overlap_cj |> 
  select(tdpKD, srsf3KD, doubleKD, nygc_srsf3kd) |> 
  eulerr::euler() |> 
  plot(quantities = TRUE, main = "Number of cyptic junctions", 
       legend = list(labels = c("TDP43KD", "SRSF3KD", "DoubleKD", "NYGC_SRSF3KD")))

cj_longer = overlap_cj |> 
  pivot_longer(cols = ends_with("KD"),
               names_to = "condition",
               values_to = "value") |> 
  filter(value == TRUE)

cj_longer |> 
  ggplot(aes(x = condition, fill = junc_cat)) +
  geom_bar() +
  theme_minimal() +
  labs(
    title = "Number of cryptic junctions",
    fill = "Junction type"
  ) +
  scale_fill_brewer(palette = "Set3")
















# Select cryptic events from NYGC dataset
  # load in dataset
nygc = read_parquet("/Users/Ewann/splicing_comparison/data/nygc/selective_cryptic_psi_in_nygc.parquet")

gencode = as.data.table(import.bed("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.bed12"))

cluster_mount_point = "/Users/Ewann/cluster"
tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v42.tx2gene.csv")
tx2gene <- fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

# generate a seq2gene/tx table
gencodes = gencode |> 
  left_join(tx2gene, by = c("name" = "TXNAME")) |> 
  left_join(gentab, by = c("GENEID" = "Gene.ID")) |>  
  select(seqnames, start, end, Gene.Symbol, name) |> 
  mutate(gene = Gene.Symbol, .keep = "unused") |> 
  mutate(txname = name, .keep = "unused") |> 
  GRanges() |> 
  unique()


# select tdp-pathology patients, get rid of other columns
nygc_tdp = nygc |>
  filter(tdp_path == "path") |>
  filter(disease_tissue == TRUE) |> 
  select(paste_into_igv_junction, psi, tissue_clean, disease) 

# create a GRange list for psi_tdp
gene_tdp = data.frame(str_split_fixed(nygc_tdp$paste_into_igv_junction, "[:-]", 3)) 
names(gene_tdp) <- c("seqnames", "start", "end")
psi_tdp = cbind(nygc_tdp, gene_tdp) |> select(-paste_into_igv_junction) |> GRanges()


# get gene names for psi_tdp
overlaps_gene_tdp = findOverlaps(psi_tdp, gencodes)
uniqueHits = overlaps_gene_tdp |> 
  as.data.frame() |> 
  distinct(queryHits, .keep_all = TRUE)
overlaps_tdp = psi_tdp[uniqueHits$queryHits]
mcols(overlaps_tdp)$gene = mcols(gencodes[uniqueHits$subjectHits])$gene


# get annotated intronics - location + gene names
gtf_path = "/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.gtf.gz"
gtf_obj =  GenomicFeatures::makeTxDbFromGFF(gtf_path,format = 'gtf')

intronics = GenomicFeatures::intronsByTranscript(gtf_obj,use.names = TRUE)
intronics = unlist(intronics) 
intronics$transcript_id = names(intronics)
intronics_gene = intronics |> 
  as.data.table() |> 
  left_join(tx2gene, by = c("transcript_id" = "TXNAME")) |> 
  left_join(gentab, by = c("GENEID" = "Gene.ID")) |> 
  select(-Biotype, -Source) |> 
  GRanges()

# subsetting selective genes from NYGC dataset
intronics_gene = intronics_gene[which(intronics_gene$Gene.Symbol %in% overlaps_tdp$gene)]


