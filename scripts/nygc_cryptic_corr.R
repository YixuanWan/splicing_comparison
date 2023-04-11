library(arrow)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(stringr)


# load in dataset
nygc = read_parquet("/Users/Ewann/splicing_comparison/data/nygc/selective_cryptic_psi_in_nygc.parquet")
gencode = data.table::as.data.table(import.bed("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.bed12"))

gentab = read.table("/Users/Ewann/splicing_comparison/data/GRCh38/gencode.v42.annotation.gtf.genes.tab", sep = "\t", header = TRUE)

tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v42.tx2gene.csv")
tx2gene <- data.table::fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")

gencodes = gencode |> 
  left_join(tx2gene, by = c("name" = "TXNAME")) |> 
  left_join(gentab, by = c("GENEID" = "Gene.ID")) |> 
  # mutate(seq = paste0(seqnames, ":", start, "-", end)) |>   
  # select(seq, Gene.Symbol) |>
  select(seqnames, start, end, Gene.Symbol) |> 
  mutate(gene = Gene.Symbol, .keep = "unused")


nygc_gene = data.frame(str_split_fixed(nygc$paste_into_igv_junction, "[:-]", 3)) 
names(nygc_gene) <- c("seqnames", "start", "end")


gr_gencodes = GRanges(gencodes)
gr_nygc = GRanges(nygc_gene)
overlaps = findOverlaps(gr_gencodes, gr_nygc_gene)
overlapping_intervals = gr_gencodes[queryHits(overlaps)]

nygc_gene$start <- as.numeric(nygc_gene$start)
nygc_gene$end <- as.numeric(nygc_gene$end)

gene_overlap = as.data.table(overlapping_intervals) |> unique() |> right_join(nygc_gene) |> filter(!is.na(gene))

# psi_tdp =
  nygc |>
  cbind(nygc_seq) |> 
  filter(tdp_path == "path") |> 
  select(psi, seqnames, start, end) |> 
  left_join(gencodes) |> 
  filter(gene == "FUS")
  
  
  group_by(paste_into_igv_junction) |>
  pivot_wider(names_from = paste_into_igv_junction, 
              names_repair = "check_unique", 
              values_from = psi,
              values_fn = list) |> 
  purrr::map(as.data.table)

psi_wide = do.call(cbind, psi_tdp) |> setnames(names(psi_tdp))

cor_matrix = psi_wide |> cor() 
diag(cor_matrix) = 0

# Convert the correlation matrix to a long format data frame
cor_df = as.data.frame(as.table(cor_matrix))|> 
  arrange(desc(Freq)) 

# Identify the most correlated cryptics
top_cor = cor_df |>  
  filter(Freq > 0.5 | Freq < -0.5)

# Create the heatmap using ggplot2
top_cor |>   
ggplot(aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "purple") +
  labs(x = "", y = "", fill = "Correlation") +
  theme_minimal() +
  # theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# write out the results
write.csv(cor_df, "/Users/Ewann/splicing_comparison/results/nygc/cryptic_tdp_corr.csv", row.names = FALSE)
  