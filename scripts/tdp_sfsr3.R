library(GenomicFeatures)
library(annotatr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(GenomicRanges)


source("scripts/snapcount_sra_cryptics.R")

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


