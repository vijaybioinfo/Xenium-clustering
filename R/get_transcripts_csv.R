#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-24
# Version 1.0.0
## This script will only process cellranger output and create the transcript.csv file needed. 
# ---

######################
# Clustering: Seurat: Spatial analysis # Step 1: get_transcripts_csv
######################
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)
options(future.globals.maxSize= 13312*1024^2)

library(optparse)
optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
    optparse::make_option(
    opt_str = c("-v", "--verbose"), default = TRUE,
    help = "Verbose: Show progress."
  )
)

# Getting arguments from command line
opt <- optparse::parse_args(optparse::OptionParser(option_list = optlist))

resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/fcastaneda/bin/quality_control/utilities.R",
  "/home/kmlanderos/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/clustering/R/plotting.R",  # variable_features_report
  #"/home/fcastaneda/bin/clustering/R/plotting.R" #plots of clustering pipeline
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R", # stats_summary_table
  "/home/ciro/scripts/figease/figease.R"
)
for(i in resources){ source(i) }

get_transcripts_csv<-function(cellranger_path){ 
    cat ("transcripts.csv file not present \n"); 
    cat ("Converting transcripts.parquet to csv as required for Seurat \n"); 
    # cat ("Creating transcripts files \n"); 
      ## Thanks cellranger https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/example-code 13 Nov 2024
      # cellranger_path <- cellranger_out
        library(arrow)
        # Path to your parquet file, edit path to where parquet file saved
        PATH <- paste0(cellranger_path, '/transcripts.parquet')
        # Edit path and output name for new file
        OUTPUT <- gsub('\\.parquet$', '.csv', PATH)
        # Specify chunk size
        CHUNK_SIZE <- 1e6

        cat("Read in the parquet file \n")
        # This reads in the data as an arrow table
        parquet_file <- arrow::read_parquet(PATH, as_data_frame = FALSE)
        # To read in the table as a tibble, set data frame to true:
        # parquet_file <- arrow::read_parquet(PATH, as_data_frame = TRUE)
        #Optional: convert parquet data frame to CSV.

        cat("Writting: convert parquet data frame to CSV. \n") 
        start <- 0
        while(start < parquet_file$num_rows) {
          end <- min(start + CHUNK_SIZE, parquet_file$num_rows)
          chunk <- as.data.frame(parquet_file$Slice(start, end - start))
          data.table::fwrite(chunk, OUTPUT, append = start != 0)
          start <- end
        }
        if(require('R.utils', quietly = TRUE)) {
          R.utils::gzip(OUTPUT)
        }

}

# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/spatial_pilot_5k/scripts/devel/config.yaml"
  opt$verbose <- TRUE
}

config = yaml::read_yaml(opt$yaml)
if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)
outdir_sp = "." #This need to change to five the user an option to run it without the snakemake
config = yaml::read_yaml(opt$yaml)
cellranger_out= config$input_expression
setwd(outdir_sp)
cat("Working in:", getwd(), "\n")

if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
if(opt$verbose) cat("Current directory: ", getwd() , "\n") ;

## Chekinf if parquet files are already in csv as required for Seurat
if(!dir.exists(config$input_expression)) stop("Cellranger input and yaml input expression doesn't match") # probably this is unneessary in the snakemake
files_cellranger <- list.files(cellranger_out)
if(!any(grepl("transcripts.csv.gz", files_cellranger))) { get_transcripts_csv(cellranger_out) } 
writeLines("transcripts.csv.gz from parquet file processed", ".transcripts_done.txt")


if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
