
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(dplyr,tidyr,tximport,rlang,DESeq2,data.table,annotables)
#directory set up 
cluster_mount_point = "/Users/Ewann/cluster/"
salmon_quant_directory = file.path("/Users/Ewann/cluster/first_weeks/splicing_comparisons/dz_sy5y/salmon_quant")

# import files
metadata_filepath = file.path("/Users/Ewann/cluster/first_weeks/splicing_comparisons/rna_seq_snakemake/config/sample_sheet_dz_sy5y.csv")
tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv")
tx2gene <- data.table::fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")


# ============================ section 1: import data ===========================

#(1) First read in the metadata
metadata_orig = read.csv(metadata_filepath,header=TRUE)
column_name = "group"
baseline = "dz_ctl"
contrast = "sy5y_ctl"

controls_name = "dz_ctl"
contrast_name = "sy5y_ctl"
#What folder and what  will my outputs be written to? #CHANGE THIS
outputdir = "/Users/Ewann/splicing_comparison/results/deseq_dz_sy5y"
exp = paste0(contrast_name,"|",controls_name)
output_path = paste0(outputdir,"/",exp)

#Make a new minimal metadata file containing the sample for this comparison

metadata = metadata_orig %>%  
    filter(is.na(exclude_sample_downstream_analysis)) %>% 
    dplyr::select(sample_name, !!(column_name)) %>% #!!(column_name) - this is called 'unquoting
    mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline', #here we need as symbol due to the logic of case_when
                                            !!as.symbol(column_name) == contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) %>% 
    filter(!is.na(comparison_condition))   


#(2) Generate a vector of the wanted file names.
files = unique(file.path(salmon_quant_directory,metadata$sample_name,"quant.sf")) 
names(files) = unique(metadata$sample_name)


#(3) To check if all the files exist
if(all(file.exists(files)) == FALSE) {
    stop("It seems that I cannot find those files...Please check if your directory is correct.")
}

#import transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages
#output matrix: average transcript length, weighted by sample-specific transcript abundance estimates 
txi.tx <- tximport(files, 
                   type="salmon", 
                   tx2gene=tx2gene,
                   ignoreTxVersion = TRUE,
                   ignoreAfterBar = TRUE,
                   txOut = TRUE)

# ========================================== section 3b: WRITE out transcript counts (optional) =============================================================

#For downstream analysis it can sometimes
#be helpful to have a table of all the transcript counts
transcript_counts = as.data.frame(txi.tx$abundance) %>% 
    tibble::rownames_to_column('TXNAME') %>% 
    left_join(tx2gene)


fwrite(transcript_counts,paste0(output_path,".transcript_counts.csv"))

txi.sum <- summarizeToGene(txi.tx, tx2gene)

# ========================================== section 4: RUN A DEFAULT DESEQ 2 =============================================================


dds = DESeqDataSetFromTximport(txi.sum,
                               colData = metadata,
                               design = ~ comparison_condition) 


# 'Note that the tximport-to-DESeq2 approach uses estimated gene counts from the transcript abundance quantifiers, 
# but not normalized counts' -- <Deseq2 vignette> (just a note - Deseq() wraps the normalization step inside)
# perform the Deseq function
dds = DESeq(dds)

results_table = results(dds) %>%  
    as.data.frame() %>% 
    tibble::rownames_to_column('ensgene') %>%   
    #remove fractions in ensgene
    mutate(ensgene = gsub("\\..*","",ensgene))  %>% 
    left_join(annotables::grch38 %>% dplyr::select(ensgene,symbol)) %>% 
    dplyr::rename(gene_name = symbol) %>% 
    mutate(experiment = "humphrey_i3cortical") %>% 
    unique() %>% as.data.table()

# Now, write everything out

output_path=paste0(outputdir,"/",exp)

fwrite(results_table,paste0(output_path,".DESEQ2_results.csv"))

#Also save the dds object to later capture
#more metadata
saveRDS(dds, file.path(paste0(output_path,".DESEQ2_object.RDS")))

