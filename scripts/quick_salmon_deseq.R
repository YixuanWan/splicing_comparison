
if (!require("pacman")) install.packages("pacman")
library(pacman)
pacman::p_load(dplyr,tidyr,tximport,rlang,DESeq2,data.table,annotables)
#directory set up #CHANGE THIS - Make sure you have the cluster mounted as a drive
cluster_mount_point = "/Users/Ewann/cluster/"
salmon_quant_directory = file.path("/Users/Ewann/splicing_comparison/humphreycorticali3/salmon")
#This is a sample sheet containing metadata about the samples
metadata_filepath = file.path("/Users/Ewann/splicing_comparison/humphreycorticali3/majiq_sample_sheet.csv")
#This is a flat file containing the mapping between transcripts and genes
#call head() on it to see what it looks like
tx2gene = file.path(cluster_mount_point,"vyplab_reference_genomes/annotation/human/GRCh38/gencode.v40.tx2gene.csv")
#read in the tx2gene file
tx2gene <- data.table::fread(tx2gene,header=FALSE)
colnames(tx2gene) = c("TXNAME", "GENEID")
#BONUS What gene has the second most transcripts? What does it do? 
#check out its expression on https://www.proteinatlas.org/
    #ENSG00000147168 -> IL2RG, encodes interleukin 2 receptor subunit gamma, expressed in T-cells, mainly for adaptive immune response
# ============================ section 1: import data ===========================

#(1) First read in the metadata
metadata_orig = read.csv(metadata_filepath,header=TRUE)
#call head() on this, what does it contain? -> sample name, fastq1&2, group, experiment, 
#Which column from the metadata should I read from?
column_name = "group"
#Which values in that column refer to baseline and contrast conditions?
baseline = "control"
contrast = "tdp43KD"
#What am I going to call those on the final output file?
controls_name = "control"
contrast_name = "humphrey_i3cortical"
#What folder and what  will my outputs be written to? #CHANGE THIS
outputdir = "/Users/Ewann/splicing_comparison/humphreycorticali3"
exp = paste0(contrast_name,"|",controls_name)
output_path = paste0(outputdir,"/",exp)

#Make a new minimal metadata file containing the sample for this comparison
#heads up! there's some kinda advance R code below turning a variable (column_name) into 
#quotation in order to be able to select from a data.frame using what's stored in the variable

metadata = metadata_orig %>%  
    filter(is.na(exclude_sample_downstream_analysis)) %>% 
    dplyr::select(sample_name, !!(column_name)) %>% #!!(column_name) - this is called 'unquoting
    mutate(comparison_condition = case_when(!!as.symbol(column_name) == baseline ~ 'baseline', #here we need as symbol due to the logic of case_when
                                            !!as.symbol(column_name) == contrast ~ 'contrast',
                                            TRUE ~ NA_character_)) %>% 
    filter(!is.na(comparison_condition))   

#Why does this part ^ exist? What did the whole metadata_orig look like?
    #1. selecting the target experiments for analysis, and add a comparison_condition column
    #2. Incorporating unit into sample_name and deleted unit column

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

