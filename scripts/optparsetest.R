library(optparse)
library(data.table)
library(dplyr)
library(stringr)


option_list = list(
  make_option(c("-s", "--salmon_quant"), type = "character", default = NULL, help = "file path for the folder containing salmon_quant output", metavar = "character"),
  make_option(c("-m", "--samplesheet"), type = "character", default = NULL, help = "file path for the sample sheet", metavar = "character"),
  make_option(c("-g", "--tx2gene"), type = "character", default = NULL, help = "file path for reference tx2gene file", metavar = "character"),
  make_option(c("-c", "--column_name"), type = "character", default = NULL, help = "column to read from the sample sheet", metavar = "character"),
  make_option(c("-b", "--baseline"), type = "character", default = NULL, help = "the value refer to baseline condition", metavar = "character"),
  make_option(c("-x", "--contrast"), type = "character", default = NULL, help = "the value refer to contrast condition", metavar = "character"),
  make_option(c("-e", "--baseline_name"), type = "character", default = NULL, help = "baseline name on output files", metavar = "character"),
  make_option(c("-t", "--contrast_name"), type = "character", default = NULL, help = "contrast name on output files", metavar = "character"),
  make_option(c("-o", "--outputdir"), type = "character", default = NULL, help = "directory for DESeq2 outputs", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


