library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tximport)

# function to run deseq2 analysis without writing results into a file
salmon_deseq <- function(cluster_mount_point,
                              tx2gene_path,
                              salmonpath,
                              metadatapath,
                              baseline,
                              contrast,
                              AbundanceEstimates,
                              celltype)
{
  
  # import files
  salmon_quant_directory = file.path(cluster_mount_point, salmonpath)
  metadata_filepath = file.path(cluster_mount_point, metadatapath)
  
  metadata_orig = read.csv(metadata_filepath,header=TRUE)
  column_name = "group"

  
  # identify baseline and contrast
  metadata = metadata_orig |>   
    filter(is.na(exclude_sample_downstream_analysis)) |> 
    unique(by = "sample_name") |> 
    mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline',
                                            !!as.symbol(column_name) == contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) |> 
    filter(!is.na(comparison_condition))
  
  
 
  # modifications specific to sy5y datasets
  if(celltype == 'sy5y'){
    files = unique(file.path(salmon_quant_directory,paste0(metadata$unit,"_", metadata$sample_name),"quant.sf"))
  } 
  else{
    files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf"))
  }
  
  names(files) = unique(metadata$sample_name)
  
  if(all(file.exists(files)) == FALSE) {
    stop("It seems that I cannot find those files...Please check if your directory is correct.")
  }
  
  # import tx2gene file
  tx2gene = file.path(cluster_mount_point, tx2gene_path)
  tx2gene <- data.table::fread(tx2gene,header=FALSE)  
  colnames(tx2gene) = c("TXNAME", "GENEID")
  
  # perform transcription counts
  txi.tx <- tximport(files, 
                     type="salmon", 
                     tx2gene=tx2gene,
                     # countsFromAbundance = AbundanceEstimates,
                     countsFromAbundance = 'no',
                     ignoreTxVersion = TRUE,
                     ignoreAfterBar = TRUE,
                     txOut = TRUE)
  
  # txi.sum <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance = AbundanceEstimates)
  txi.sum <- summarizeToGene(txi.tx, tx2gene, countsFromAbundance = 'no')
  
  dds = DESeqDataSetFromTximport(txi.sum,
                               colData = metadata,
                               design = ~ comparison_condition)

  dds = DESeq(dds)
  
  return(dds)
}

# ==================================================================================================

# run deseq with sy5y 
sy5y_deseq = salmon_deseq("/Users/Ewann/cluster",
                          "vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv",
                          "first_weeks/TDP_CHX_CLONES_GLIA/salmon",
                          "first_weeks/TDP_CHX_CLONES_GLIA/sample_sheet.csv",
                          baseline = "no_dox",
                          contrast = "dox_075",
                          "lengthScaledTPM",
                          celltype = "sy5y")



# run deseq with dz
dz_deseq = salmon_deseq("/Users/Ewann/cluster",
                        "vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv",
                        "first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/salmon_quant",
                        "first_weeks/TDP_curves_UPF1_GLIA/DZ_curves/sample_sheet_dz_curves.csv",
                        baseline = "DZ_curves_0",
                        contrast = "DZ_curves_1",
                        "lengthScaledTPM",
                        celltype = "dz")

# humphrey i3
humphrey_deseq = salmon_deseq("/Users/Ewann/cluster",
                              "vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv",
                              "TDP43_RNA/i3_Cortical_long_reads/new_data_nanopore_plus_sr/salmon",
                              "first_weeks/splicing_comparisons/symbams/humphrey_i3cortical/majiq_sample_sheet.csv",
                              baseline = "control",
                              contrast = "tdp43KD",
                              "lengthScaledTPM",
                              celltype = "i3")
# sedigghi i3
sedigghi_deseq = salmon_deseq("/Users/Ewann/cluster",
                              "vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv",
                              "first_weeks/WARD_BAMS_NEW/salmon",
                              "first_weeks/WARD_BAMS_NEW/ward_i3_newer_longer_bams_sample_sheet.csv",
                              baseline = "CTL",
                              contrast = "TDP43_KD",
                              "lengthScaledTPM",
                              celltype = "sedigghi_i3")
# brown i3
brown_deseq = salmon_deseq("/Users/Ewann/cluster",
                      "vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv",
                      "alb_projects/data/ward_bams/salmon",
                      "alb_projects/data/ward_bams/ward_bams.csv",
                      baseline = "control",
                      contrast = "tdpKD",
                      "lengthScaledTPM",
                      celltype = "i3")


# get normalised counts
sy5y_counts = DESeq2::counts(sy5y_deseq, normalized = TRUE)  
dz_counts = DESeq2::counts(dz_deseq, normalized = TRUE)
humphrey_counts = DESeq2::counts(humphrey_deseq, normalized = TRUE)
sedigghi_counts = DESeq2::counts(sedigghi_deseq, normalized = TRUE)
brown_counts = DESeq2::counts(brown_deseq, normalized = TRUE)

# modify the table
melt_table <- function(file){
  df = file |> as.data.frame() |> 
            tibble::rownames_to_column('gene') |> 
            mutate(ensmbleID = gsub("\\..*", "", gene)) |>  
            select(-gene) |> 
            janitor::clean_names() |>
            as.data.table() |> 
            melt(id.vars ="ensmble_id")
  return(df)
}

sy5y_long = melt_table(sy5y_counts)
dz_long = melt_table(dz_counts)
humphrey_long = melt_table(humphrey_counts)
sedigghi_long = melt_table(sedigghi_counts)
brown_long = melt_table(brown_counts)

# extract sample names for identifying celltypes
sy5y_id = sy5y_long$variable |> unique()
dz_id = dz_long$variable |> unique()
humphrey_id = humphrey_long$variable |> unique()
sedigghi_id = sedigghi_long$variable |> unique()
brown_id = brown_long$variable |> unique()

# =============================
# table for plots - all groups
counts_long = rbind(sy5y_long, dz_long, humphrey_long, sedigghi_long, brown_long) |> 
  mutate(celltype = case_when(variable %in% dz_id ~ 'SK-N-DZ',
                              variable %in% sy5y_id ~ 'SH-SY-5Y',
                              variable %in% humphrey_id ~ 'Humphrey-i3',
                              variable %in% sedigghi_id ~ 'Sedigghi-i3',
                              variable %in% brown_id ~ 'Brown-i3',
                              TRUE ~ NA_character_))  |> 
  filter(!is.na(celltype)) 

# plot tdp43 level - all groups
counts_long |> 
  filter(ensmble_id == "ENSG00000120948") |> 
  ggplot(aes(x = celltype, y = value)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = celltype)) +
  theme_minimal()

tdp = counts_long |> 
  filter(ensmble_id == "ENSG00000120948") 
