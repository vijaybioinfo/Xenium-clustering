#!/usr/bin/R

# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# Taking as a base Ciro's single-cell clustering pipeline (https://github.com/vijaybioinfo/clustering)
# Date: 2025-06-24
# Version 1.0.0
    # Seurat object creation, first QC's, metadata filtering, normalization and high variable genes (HVG slection), PC and UMAP
# ---

######################
# Clustering: Seurat: Spatial analysis #
######################
options(future.globals.maxSize= 13312*1024^2)
lib4.3path = "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"
if(grepl("4.3", getRversion()) && file.exists(lib4.3path))
  .libPaths(new = lib4.3path)

library(optparse)
optlist <- list(
  optparse::make_option(
    opt_str = c("-y", "--yaml"), type = "character",
    help = "Configuration file: Instructions in YAML format."
  ),
  optparse::make_option(
    opt_str = c("-p", "--percent"), type = "numeric",
    help = "Percentage of variance."
  ),
  optparse::make_option(
    opt_str = c("-n", "--n_comp"), type = "numeric", default = 50,
    help = "Total number of components to explore."
  ),
  optparse::make_option(
    opt_str = c("-c", "--chosen_comp"), type = "numeric",
    help = "Chosen components."
  ),
  optparse::make_option(
    opt_str = c("-r", "--prefix"), type = "character",
    help = "Prefix for output. Name of the files."
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

#source("/home/ciro/scripts/figease/figease.R")
library(ggplot2)
library(reshape2)
library(dplyr)
library(Seurat)
library(parallel)
library(Seurat)
library(plotly)

# opt parameters have the priority
if(interactive()){ # Example/manually
  opt$yaml = "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/spatial_pilot_5k/scripts/asthma/config_SaHe_5k_run7_donor_specific_cd4_cd8_spa.yaml"
  opt$verbose <- TRUE
  opt$percent <- 30
  opt$n_comp <- 40
  opt$chosen_comp <- 20
  opt$prefix <- "init_mean0.01_pct30_pc30"
}

config = yaml::read_yaml(opt$yaml)
cellranger_out= config$input_expression
pcts= opt$percent
pc= opt$chosen_comp
opt$prefix <- paste0(".object_", opt$prefix, ".rds")
if(interactive()) outdir_sp= paste0(config$output_dir, "/", config$project_name)
outdir_sp = "."
cat("Working in:", getwd(), "\n")

if(opt$verbose) cat('Date and time:\n') ; st.time <- timestamp();
## Chekinf if parquet files are already in csv as required for Seurat
if(!dir.exists(outdir_sp)) dir.create(outdir_sp, recursive=TRUE)
setwd(outdir_sp)

# Load the Xenium data

#################################################
################# Preprocessing #################
#################################################

    cat("------- Reading object -------- \n")
    xenium.obj <- LoadXenium(cellranger_out)
    
    cat("------- Overal QCs -------- \n")

    dir.create("Section_qc")
    pdf("Section_qc/Tissue_section_complete.pdf", width=10, height=5)
    print(ImageFeaturePlot(xenium.obj, features = c("nCount_Xenium"), border.size=NA, max.cutoff=400, axes=TRUE))
    dev.off()

    # rm(xenium.obj)
    # xenium.obj<-xenium.obj_pbt
    # rm(xenium.obj_pbt)

    pdf(paste0("Section_qc/QCs_prefiltered_", ncol(xenium.obj), "_cells.pdf"))
    # xenium.obj <- PercentageFeatureSet(xenium.obj, pattern = "^MT-", col.name = "percent.mt")
    aver<-VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size=0.001, alpha=0.2, combine= FALSE)
    print(aver[[1]] + geom_boxplot(width=0.2))
    print(aver[[2]] + geom_boxplot(width=0.2))
    dev.off()

    pdf("Section_qc/UMI_per_gene_distribution.pdf", width=10, height=5)
    umi_per_genes<-rowSums(xenium.obj@assays$Xenium) %>% melt() %>% arrange(value)
    p<-ggplot(umi_per_genes, aes(x = log10(value))) + geom_density(color=4) + labs(x= "log10(UMI count per gene)") + theme_classic()
    print(p)
    dev.off()

    write.csv(umi_per_genes, "Section_qc/UMIs_per_gene.csv")

    umi_cells<-melt(table(xenium.obj$nCount_Xenium))
    xliC<-round(mean(xenium.obj$nCount_Xenium))
    pdf("Section_qc/UMI_per_cell_distribution.pdf", width=10, height=5)
    p<-ggplot(umi_cells, aes(y=value, x=Var1)) +
      geom_bar(color="blue", fill="blue", stat="identity")+ labs(y= "Number of cells", x="#UMI") + theme_classic() + geom_vline(xintercept = xliC, linetype="dotted", color = "red", size=1.5) + labs(subtitle = paste("mean nCount: ", xliC))
    print(p)
    dev.off()

    genes_cells<-melt(table(xenium.obj$nFeature_Xenium))
    xliF<-round(mean(xenium.obj$nFeature_Xenium))
    xliC<-sum(xenium.obj$nFeature_Xenium > xliF)

    pdf("Section_qc/genes_per_cell_distribution.pdf", width=10, height=5)
    p<-ggplot(genes_cells, aes(y=value, x=Var1)) +
      geom_bar(color="blue", fill="blue", stat="identity")+ labs(y= "Number of cells", x="# genes with > 3 UMI") + theme_classic() + geom_vline(xintercept = xliF, linetype="dotted",  color = "red", size=1.5) + labs(subtitle = paste("mean nFeatures: ", xliF, ". \nCells with more than mean: ", xliC))
    print(p)
    dev.off()

    meta_data <-  xenium.obj@meta.data

    if(!is.null(config$metadata)) {
      if(opt$verbose) cat(" ------- Adding metadata -------  \n")
      addmetadataf <- if(length(config$metadata) > 0) config$metadata else "no_file"
      addmetadataf <- path.expand(addmetadataf)
      tvar <- file.exists(addmetadataf)
      if(any(tvar)){
        addmetadataf <- addmetadataf[tvar]
        if(opt$verbose) cat("Adding metadata from:", addmetadataf, sep = "\n")
        for(addmetadataf_i in addmetadataf){
          # addmetadataf_i <- addmetadataf[1]
          addannot <- if(grepl(".csv", addmetadataf_i)) remove.factors(readfile(addmetadataf_i, row.names=1)) else remove.factors(readfile(addmetadataf_i))
          print(dim(meta_data))
          print(dim(addannot))
          meta_data <- meta_data[rownames(meta_data) %in% rownames(addannot),] # keeping only the cells present in both metadatas. Thanks Donas!
          # partial matching problem addressed in /home/ciro/covid19/scripts/partial_matching.R
          meta_data <- joindf(x = meta_data, y = addannot)
        }
        print(dim(meta_data))
      }
    } else cat("No extra metadata given \n")

   ### Checking
    if(!is.null(config$metadata_donor)){
      if(opt$verbose) cat(" ------- Adding donor metadata -------  \n")
        tvar <- as.character(meta_data$orig.TMA_core)
        maxln <- max(sapply(tvar, length))
        if(opt$verbose) cat("Creating donor metadata with", maxln, "columns.\n")
        meta_donor <- data.table::rbindlist(lapply(tvar, function(x){
          if(length(x) < maxln && length(x) == 1) x <- rep(x, length.out = maxln)
          if(length(x) < maxln) x <- rep(paste0(x, collapse = "-"), length.out = maxln)
          as.data.frame(t(x))
        }))
        meta_donor <- remove.factors(
          data.frame(meta_donor, row.names = rownames(meta_data))
        )
        if(opt$verbose){
          print(head(meta_donor)); print(tail(meta_donor))
          cat("Checking created columns\n")
          print(lapply(meta_donor, table, useNA = 'always'))
        }

          rnname <- ifelse(grepl("~", config$metadata_donor), gsub(".*~", "", config$metadata_donor), 1)
          config$metadata_donor <- gsub("~.*", "", config$metadata_donor)
          tmp <- paste0(config$metadata_donor, "\n- Using row names:", rnname)
          if(opt$verbose) cat("Adding given metadata: ", tmp, "\n")
          tmp <- read.csv(config$metadata_donor, stringsAsFactors = FALSE, row.names = rnname)
          tmp[tmp == ""] <- NA; mdonor_extra <- tmp[meta_donor$V1, , drop = FALSE]
          rownames(mdonor_extra) <- rownames(meta_donor)
          if(opt$verbose) print(sapply(mdonor_extra, table, useNA = 'always'))
          # meta_donor <- joindf(meta_donor, mdonor_extra)
          meta_donor <- joindf(data.frame(meta_donor), mdonor_extra) #For JaSh
          colnames(meta_donor) <- paste0("orig.", colnames(meta_donor))
          if(opt$verbose) print(head(meta_donor)); print(tail(meta_donor))
          meta_data <- joindf(meta_data, meta_donor)
    }
        
    if(!is.null(config$filtering$subset)) {
      cat('  -------  Filters present  -------  \n')
      cat('----------------------- Filtering data --------------------------\n')
        filtereddata <- filters_complex(
          mdata = meta_data,
          filters = lapply(names(config$filtering$subset), function(x) c(x, config$filtering$subset[[x]]) ),
          verbose = opt$verbose
        )
        cat("Preserving:", nrow(filtereddata[[1]]), "/", nrow(meta_data), "samples/cells\n")
        meta_data <- filtereddata[[1]]; rm(filtereddata)
        xenium.obj <- subset(xenium.obj, cells =  rownames(meta_data))
        meta_data <- meta_data %>% select(!matches("seurat_clus|^Xenium_snn_res.|^SCT_snn_res.|umap_1|umap_2"))
        xenium.obj <- AddMetaData(
        object = xenium.obj,
        metadata = meta_data)
        if(opt$verbose) str(meta_data)
    } else cat("No filters present ")

    pdf(paste0("Section_qc/QCs_postfiltered_", ncol(xenium.obj), "_cells.pdf"))
    aver<-VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size=0.001, alpha=0.2, combine= FALSE)
    print(aver[[1]] + geom_boxplot(width=0.2))
    print(aver[[2]] + geom_boxplot(width=0.2))
    dev.off()

    pdf("Section_qc/QC_in_Tissue_plot.pdf", width=10, height=5)
    print(ImageFeaturePlot(xenium.obj, fov = "fov", features = c("nCount_Xenium"), border.size=NA, max.cutoff=400))
    print(ImageFeaturePlot(xenium.obj, fov = "fov", features = c("nFeature_Xenium"), border.size=NA))
    # ImageFeaturePlot(xenium.obj, fov = "fov", features = c("nFeature_Xenium"), border.size = 0.00000001)
    dev.off()


#################################################
################# CLUSTERING ####################
#################################################

  #### CLUSTERING 
    cat("------- Clustering analysis -------- \n")
 
   ### Filter highly variable features
  hvg_final <- rownames(xenium.obj)
  tvar <- config$variable_features$file

  if(tvar != "no_file"){
    if(file.exists(tvar)){
      if(opt$verbose) cat("WARINING: Subsetting features from file:\n", tvar, "\n")
      file_con <- file(description = tvar, open = "r")
      these_feats <- readLines(con = file_con); close(file_con)
      hvg_final <- hvg_final[hvg_final %in% these_feats]
    } else{
      stop("Error: file config$variable_features$file not exists")
    }
  }

    if(!grepl(config$norm, "sctransform")){
      if(opt$verbose) cat('\n@@@@@@@@@ Normalising (NormalizeData) ...\n'); timestamp()
      xenium.obj <- NormalizeData(
        object = xenium.obj,
        normalization.method = config$norm,
        scale.factor = 10000,
        verbose = opt$verbose
      )  
      if(opt$verbose) cat('\n@@@@@@@@@ Selecting features (FindVariableFeatures) ...\n'); timestamp()
      xenium.obj <- FindVariableFeatures(
        object = xenium.obj,
        selection.method = config$variable_features$method,
        nfeatures = config$variable_features$nfeatures,
        mean.cutoff = config$variable_features$mean.cutoff,
        dispersion.cutoff =  config$variable_features$dispersion.cutoff,
        verbose = opt$verbose
      )
      hvg_df <- HVFInfo(xenium.obj, selection.method = "vst")
      disp_n <- grep("standardized", colnames(hvg_df), value=TRUE); hvg_df <- hvg_df[order(-hvg_df[, disp_n]), ] # order
      # mean_pct_filters <- TRUE
      mean_pct_filters <- !isTRUE(config$variable_features$file_only)
      # hvg_final <- rownames(xenium.obj)
    }else{
      if(opt$verbose) cat('\n@@@@@@@@@ SCTransform - no regression\n'); timestamp()
      xenium.obj <- SCTransform(
        object = xenium.obj,
        variable.features.n = config$variable_features$nfeatures,
        vars.to.regress = NULL,
        conserve.memory = FALSE,
        return.only.var.genes = TRUE,
        verbose = opt$verbose
      )
      if(opt$verbose) cat("Attention: If SCTransform default hvg are chosen for the clustering. Ignoring pct loop")
      pcts <- "default"
    }; if(opt$verbose) cat("Default Assay:", DefaultAssay(object = xenium.obj), "\n")
  

    # pcts= opt$percent
    clustering_pct<-pcts
    verbose <- TRUE
    clustering_per_pcts_num<-ifelse(clustering_pct<10 || is.character(clustering_pct), paste0("0", clustering_pct), clustering_pct)

    if(!grepl(config$norm, "sctransform") && mean_pct_filters){
      if(verbose) cat("Taking", clustering_pct, "%\n")
      hvg_df_keep<-hvg_df
      hvg_df_keep$cumulative = cumsum(hvg_df_keep[, disp_n])
      hvg_df_keep$cumulative_pct = round(hvg_df_keep$cumulative / sum(hvg_df_keep[, disp_n], na.rm = TRUE) * 100, 2)
      passed <- hvg_df_keep$cumulative_pct <= clustering_pct
      if(verbose) cat("Number of features:", sum(passed), "of", nrow(hvg_df_keep), "\n")
      hvg_df_keep$passed<- hvg_df_keep$cumulative_pct <= clustering_pct & (rownames(hvg_df_keep) %in% hvg_final) ## The last condition could be  when variable "passed" is definded

      pdf(paste0(outdir_sp, "/cumulative_pct", clustering_pct, ".pdf"))
          print(ggplot(hvg_df_keep, aes(y=cumulative_pct, x=1:length(cumulative_pct), color=passed)) + geom_point() + scale_colour_manual(values = c('TRUE'="red", "FALSE"="black")) + geom_hline(yintercept = pcts, linetype='dotted', col = 'red') + labs(title="Cumulative variance of Highly Variable Genes", subtitle = paste0("Number of features: ", sum(passed), " of ", nrow(hvg_df_keep), "\n"), x='Number of genes, order: variance.standardized', y='Cumulative percentage: variance standardized'))
      dev.off()

      hvg_df_keep <- hvg_df_keep[passed, ]; #log_history <- paste0(log_history, "Pct")
    
      hvg_final2 <- intersect(hvg_final, rownames(hvg_df_keep))
      if(verbose) cat("Taking", length(hvg_final2), "/", length(VariableFeatures(xenium.obj)), "features\n")
      tmp <- length(intersect(VariableFeatures(xenium.obj), hvg_final2))
      if(verbose) cat("Overlap with initial selection:", tmp, "\n")
      VariableFeatures(object = xenium.obj, assay = NULL) <- hvg_final2
      if(verbose) cat("Final number:", length(VariableFeatures(xenium.obj)), "features\n")
    }

    if(!grepl(config$norm, "sctransform")){
      if(opt$verbose) cat('\n@@@@@@@@@ Scaling data...\n'); timestamp()
      xenium.obj <- ScaleData(
        object = xenium.obj,
        vars.to.regress = config$regress_var,
        block.size = 2000,
        verbose = opt$verbose
      )
    }

     # lapply(rev(pcs), function(clustering_per_pcs){
      # pc= opt$chosen_comp
      clustering_per_pcs<-pc
      clustering_per_pcs_num<-ifelse(clustering_per_pcs<10, paste0("0", clustering_per_pcs), clustering_per_pcs)

      if(opt$verbose){
        cat('\n@@@@@@@@@ Linear dimensional reduction\n')
        cat('Computing:', casefold(config$dim_reduction$base$type, upper = TRUE), '\n')
        cat('Components:', opt$n_comp, '\n')
      }; timestamp()
      cat("PC:", clustering_per_pcs_num,  " \n")

      xenium.obj <- RunPCA(xenium.obj, 
        npcs = opt$n_comp,
        features = VariableFeatures(xenium.obj), 
        nfeatures.print = 15,
        verbose = opt$verbose)
      
      if (grepl("harmony", config$dim_reduction$base$type, ignore.case = TRUE)) {
        xenium.obj <- harmony::RunHarmony(
        object = xenium.obj,
        group.by.vars = config$dim_reduction$base$batch,
        dims.use = 1:opt$chosen_comp,
        verbose = opt$verbose
        )
      }

      pdf(paste0("Elbow_", clustering_pct, "pct_", clustering_per_pcs, "pcs.pdf"), width=10, height=5)
      print(ElbowPlot(object = xenium.obj, ndim=clustering_per_pcs, reduction="pca") + geom_vline(xintercept = clustering_per_pcs, linetype="dotted", color = "red", size=1.5) + xlim(0,30) + ylim(1,5))
      dev.off()

      n_neis <- if(!is.null(config$dim_reduction$umap)) unique(config$dim_reduction$umap$n.neighbors)

      xenium.obj <- RunUMAP(xenium.obj, 
        reduction = config$dim_reduction$base$type,
        dims = 1:clustering_per_pcs, 
        n.neighbors = n_neis,
        min.dist = config$dim_reduction$umap$min.dist,
        verbose = opt$verbose); timestamp()

      if(opt$verbose) cat("Resolutions:", paste0(config$resolution, collapse = ", "), "\n")
      
      xenium.obj <- FindNeighbors(xenium.obj, 
      reduction = config$dim_reduction$base$type, dims = 1:clustering_per_pcs, 
      # compute.SNN = compute_snn,
      verbose = opt$verbose)

      xenium.obj <- FindClusters(xenium.obj,
        resolution = config$resolution,
        verbose = opt$verbose #, 
        # algorithm =  4 ### For leiden clustering !!! Modify after run
        )


      global_resolution <- paste0(clustering_per_pcts_num, "pct", clustering_per_pcs_num, "pc_qc")
      dir.create(global_resolution)

      xenium.obj@meta.data$cellname <- rownames(xenium.obj@meta.data)
      xenium.obj@meta.data = joindf(xenium.obj@meta.data,
              as.data.frame(xenium.obj@reductions$umap@cell.embeddings))

      filename_init <- opt$prefix   
      if(opt$verbose) cat("Saving init \n")
      init_object_name <- filename_init
      saveRDS(object = xenium.obj, file=init_object_name)

      if(opt$verbose) cat("Saving reductions \n")
      reductionsf <- gsub("init", "reductions", filename_init)
      saveRDS(object = xenium.obj@reductions, file=reductionsf)

      if(opt$verbose) cat("Saving metadata \n")
      metadataf <- gsub("init", "metadata", filename_init)
      saveRDS(object = xenium.obj@meta.data, file=metadataf)

      if(opt$verbose) cat("Saving graphs \n")
      graphsf <- gsub("init", "graphs", filename_init)
      saveRDS(object = xenium.obj@graphs, file=graphsf)

    rm(xenium.obj); rm(aver)


if(opt$verbose){
  cat('\n\n*******************************************************************\n')
  cat('Starting time:\n'); cat(st.time, '\n')
  cat('Finishing time:\n'); timestamp()
  cat('*******************************************************************\n')
  cat('SESSION INFO:\n'); print(sessionInfo()); cat("\n")
  cat('Pipeline finished successfully\n')
}
  
