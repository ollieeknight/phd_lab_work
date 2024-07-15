
plot_metrics_ribosomal <- function(seurat_object, RNA_cutoff, ribosomal_cutoff) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("Tidyverse is required but not installed")
  }
  suppressWarnings(
    ggplot(seurat_object@meta.data, aes(nCount_RNA, percent_ribosomal)) + 
      geom_bin2d(bins = 200) + 
      scale_fill_continuous(type = 'viridis') + 
      theme_classic() +
      scale_x_log10() + scale_y_log10() + 
      geom_vline(xintercept = RNA_cutoff) + 
      geom_hline(yintercept = ribosomal_cutoff)
  )
}

plot_metrics_mitochondrial <- function(seurat_object, RNA_cutoff, mitochondrial_cutoff) {
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("Tidyverse is required but not installed")
  }
  suppressWarnings(
    ggplot(seurat_object@meta.data, aes(nCount_RNA, percent_mitochondrial)) + 
      geom_bin2d(bins = 200) + 
      scale_fill_continuous(type = 'viridis') + 
      theme_classic() +
      scale_x_log10() + scale_y_log10() + 
      geom_vline(xintercept = RNA_cutoff) + 
      geom_hline(yintercept = mitochondrial_cutoff)
  )
}

map_to_tonsil_reference <- function(seurat_object) {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Azimuth is required but not installed")
  }
  if (ncol(seurat_object) > 100) {
    # Check if 'tonsil_reference' is loaded in the environment
    if (!exists("tonsil_reference", envir = .GlobalEnv)) {
      print('Loading tonsil_reference from RDS...')
      tonsil_reference <<- readRDS('~/group/work/data/tissue_progenitor/objects/DOGMA_objects.rds')
      tonsil_reference <<- tonsil_reference[['CITE_Romagnani']]
    }
    print(paste0('Mapping ', names(table(seurat_object$orig.ident)), ' to Romagnani tonsil reference'))
    print('Generating transfer anchors between datasets...')
    tryCatch({
      suppressMessages(suppressWarnings({
        anchors <- Seurat::FindTransferAnchors(
          reference = tonsil_reference, query = seurat_object,
          reference.assay = 'RNA', query.assay = 'RNA',
          reduction = 'pcaproject', reference.reduction = 'harmony_ref',
          dims = 1:30, verbose = FALSE
        )
      }))
      gc()
      
      print('Mapping query against reference...')
      suppressMessages(suppressWarnings({
        test <- Seurat::MapQuery(
          anchorset = anchors, query = seurat_object,
          reference = tonsil_reference, refdata = tonsil_reference$celltype,
          reference.reduction = 'harmony_ref', reference.dims = 1:30,
          query.dims = 1:30, reduction.model = 'umap_rna', verbose = FALSE
        )
      }))
      
      seurat_object$Romagnani_tonsil_reference_celltype <- test$predicted.id
      seurat_object$score_Romagnani_tonsil_reference_celltype <- test$predicted.id.score
      rm(anchors, test)
      gc()
    }, error = function(e) {
      # Code to handle the error
      print(paste("Error during processing:", e$message))
      seurat_object$Romagnani_tonsil_reference_celltype <- 'Unmapped'
      seurat_object$score_Romagnani_tonsil_reference_celltype <- 0
    })
  } else {
    print('Library has too few cells to annotate confidently...')
    seurat_object$Romagnani_tonsil_reference_celltype <- 'Unmapped'
    seurat_object$score_Romagnani_tonsil_reference_celltype <- 0
  }
  return(seurat_object)
}

map_to_tilp_reference <- function(seurat_object) {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Azimuth is required but not installed")
  }
  if (ncol(seurat_object) > 100) {
    # Check if 'tilp_reference' is loaded in the environment
    if (!exists("tilp_reference", envir = .GlobalEnv)) {
      tilp_reference <<- readRDS('~/group/work/data/NK_IBD/objects/GEX_Colonna_TILP.rds')
    }
    print(paste0('Mapping ', names(table(seurat_object$orig.ident)), ' to Colonna TILP reference'))
    print('Generating transfer anchors between datasets...')
    tryCatch({
      suppressMessages(suppressWarnings({
        anchors <- Seurat::FindTransferAnchors(
          reference = tilp_reference, query = seurat_object,
          reference.assay = 'RNA', query.assay = 'RNA',
          reduction = 'pcaproject', reference.reduction = 'harmony_ref',
          dims = 1:30, k.filter = NA, verbose = FALSE
        )
      }))
      gc()
      
      print('Mapping query against reference...')
      suppressMessages(suppressWarnings({
        test <- Seurat::MapQuery(
          anchorset = anchors, query = seurat_object,
          reference = tilp_reference, refdata = tilp_reference$celltype,
          reference.reduction = 'harmony_ref', reference.dims = 1:30,
          query.dims = 1:30, reduction.model = 'umap_rna', verbose = FALSE
        )
      }))
      
      seurat_object$Colonna_TILP_reference_celltype <- test$predicted.id
      seurat_object$score_Colonna_TILP_reference_celltype_score <- test$predicted.id.score
      rm(anchors, test)
      gc()
    }, error = function(e) {
      # Code to handle the error
      print(paste("Error during processing:", e$message))
      seurat_object$Colonna_TILP_reference_celltype <- 'Unmapped'
      seurat_object$score_Colonna_TILP_reference_celltype_score <- 0
    })
  } else {
    print('Library has too few cells to annotate confidently...')
    seurat_object$Colonna_TILP_reference_celltype <- 'Unmapped'
    seurat_object$score_Colonna_TILP_reference_celltype_score <- 0
  }
  return(seurat_object)
}

map_to_azimuth_pbmc_reference <- function(seurat_object, assay = 'RNA') {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Azimuth is required but not installed")
  }
  if (ncol(seurat_object) > 100) {
    print(paste0('Mapping ', names(table(seurat_object$orig.ident)), ' to Azimuth PBMC reference'))
    tryCatch({
      suppressMessages(suppressWarnings({
        mapdata <- Azimuth::RunAzimuth(seurat_object, reference = 'pbmcref', assay = assay, verbose = FALSE)
      }))
      
      print('Adding predictions to object...')
      seurat_object$azimuth_pbmc_reference_celltype_l1 <- mapdata$predicted.celltype.l1
      seurat_object$score_azimuth_pbmc_reference_celltype_l1 <- mapdata$predicted.celltype.l1.score
      seurat_object$azimuth_pbmc_reference_celltype_l2 <- mapdata$predicted.celltype.l2
      seurat_object$score_azimuth_pbmc_reference_celltype_l2 <- mapdata$predicted.celltype.l2.score
      seurat_object$azimuth_pbmc_reference_celltype_l3 <- mapdata$predicted.celltype.l3
      seurat_object$score_azimuth_pbmc_reference_celltype_l3 <- mapdata$predicted.celltype.l3.score
      gc()
    }, error = function(e) {
      # Code to handle the error
      print(paste("Error during processing:", e$message))
      seurat_object$azimuth_pbmc_reference_celltype_l1 <- 'Unmapped'
      seurat_object$score_azimuth_pbmc_reference_celltype_l1 <- 0
      seurat_object$azimuth_pbmc_reference_celltype_l2 <- 'Unmapped'
      seurat_object$score_azimuth_pbmc_reference_celltype_l2 <- 0
      seurat_object$azimuth_pbmc_reference_celltype_l3 <- 'Unmapped'
      seurat_object$score_azimuth_pbmc_reference_celltype_l3 <- 0
    })
  } else {
    print('Library has too few cells to annotate confidently...')
    seurat_object$azimuth_pbmc_reference_celltype_l1 <- 'Unmapped'
    seurat_object$score_azimuth_pbmc_reference_celltype_l1 <- 0
    seurat_object$azimuth_pbmc_reference_celltype_l2 <- 'Unmapped'
    seurat_object$score_azimuth_pbmc_reference_celltype_l2 <- 0
    seurat_object$azimuth_pbmc_reference_celltype_l3 <- 'Unmapped'
    seurat_object$score_azimuth_pbmc_reference_celltype_l3 <- 0
  }
  return(seurat_object)
}

map_to_v4_pbmc_reference <- function(seurat_object, assay = 'RNA') {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Azimuth is required but not installed")
  }
  if (ncol(seurat_object) > 100) {
    print(paste0('Mapping ', names(table(seurat_object$orig.ident)), ' to 10x 3prime v4 PBMC reference'))
    tryCatch({
      suppressMessages(suppressWarnings({
        mapdata <- Azimuth::RunAzimuth(seurat_object, reference = '~/group/work/ref/seurat/v4_reference/', assay = assay, verbose = FALSE)
      }))
      
      print('Adding predictions to object...')
      seurat_object$v4_10x_pbmc_reference_celltype_l1 <- mapdata$predicted.celltype_l1
      seurat_object$score_v4_10x_pbmc_reference_celltype_l1 <- mapdata$predicted.celltype_l1.score
      seurat_object$v4_10x_pbmc_reference_celltype_l2 <- mapdata$predicted.celltype_l2
      seurat_object$score_v4_10x_pbmc_reference_celltype_l2 <- mapdata$predicted.celltype_l2.score
      gc()
    }, error = function(e) {
      # Code to handle the error
      print(paste("Error during processing:", e$message))
      seurat_object$v4_10x_pbmc_reference_celltype_l1 <- 'Unmapped'
      seurat_object$score_v4_10x_pbmc_reference_celltype_l1 <- 0
      seurat_object$v4_10x_pbmc_reference_celltype_l2 <- 'Unmapped'
      seurat_object$score_v4_10x_pbmc_reference_celltype_l2 <- 0
    })
  } else {
    print('Library has too few cells to annotate confidently...')
    seurat_object$v4_10x_pbmc_reference_celltype_l1 <- 'Unmapped'
    seurat_object$score_v4_10x_pbmc_reference_celltype_l1 <- 0
    seurat_object$v4_10x_pbmc_reference_celltype_l2 <- 'Unmapped'
    seurat_object$score_v4_10x_pbmc_reference_celltype_l2 <- 0
  }
  return(seurat_object)
}


map_to_zhang_bmmc_reference <- function(seurat_object, assay = 'RNA') {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Azimuth is required but not installed")
  }
  if (ncol(seurat_object) > 100) {
    print(paste0('Mapping ', names(table(seurat_object$orig.ident)), ' to Zhang BMMC reference'))
    tryCatch({
      suppressMessages(suppressWarnings({
        mapdata <- Azimuth::RunAzimuth(seurat_object, reference = '~/group/work/data/datasets/BMMC_atlas/data/ref/', assay = 'RNA', verbose = FALSE)
      }))
      
      print('Adding predictions to object...')
      seurat_object$zhang_reference_celltype_l1 <- mapdata$predicted.Level1
      seurat_object$score_zhang_reference_celltype_l1 <- mapdata$predicted.Level1.score
      seurat_object$zhang_reference_celltype_l2 <- mapdata$predicted.Level2
      seurat_object$score_zhang_reference_celltype_l2 <- mapdata$predicted.Level2.score
      seurat_object$zhang_reference_celltype_l3 <- mapdata$predicted.Level3R
      seurat_object$score_zhang_reference_celltype_l3 <- mapdata$predicted.Level3R.score
      seurat_object$zhang_reference_celltype_l4 <- mapdata$predicted.Level3M
      seurat_object$score_zhang_reference_celltype_l4 <- mapdata$predicted.Level3M.score
      gc()
    }, error = function(e) {
      # Code to handle the error
      print(paste("Error during processing:", e$message))
      seurat_object$zhang_reference_celltype_l1 <- 'Unmapped'
      seurat_object$score_zhang_reference_celltype_l1 <- 0
      seurat_object$zhang_reference_celltype_l2 <- 'Unmapped'
      seurat_object$score_zhang_reference_celltype_l2 <- 0
      seurat_object$zhang_reference_celltype_l3 <- 'Unmapped'
      seurat_object$score_zhang_reference_celltype_l3 <- 0
      seurat_object$zhang_reference_celltype_l4 <- 'Unmapped'
      seurat_object$score_zhang_reference_celltype_l4 <- 0
    })
  } else {
    print('Library has too few cells to annotate confidently...')
    seurat_object$zhang_reference_celltype_l1 <- 'Unmapped'
    seurat_object$score_zhang_reference_celltype_l1 <- 0
    seurat_object$zhang_reference_celltype_l2 <- 'Unmapped'
    seurat_object$score_zhang_reference_celltype_l2 <- 0
    seurat_object$zhang_reference_celltype_l3 <- 'Unmapped'
    seurat_object$score_zhang_reference_celltype_l3 <- 0
    seurat_object$zhang_reference_celltype_l4 <- 'Unmapped'
    seurat_object$score_zhang_reference_celltype_l4 <- 0
  }
  return(seurat_object)
}

library(Seurat)

merge_seurat_GEX_objects <- function(seurat_object_list, chunk_size = 5, assays = NULL) {
  # Check if assays parameter is provided and is a non-empty vector
  if (is.null(assays) || length(assays) == 0) {
    stop("Assays parameter must be provided as a non-empty vector.")
  }
  
  merged_data_list <- list()
  
  while (length(seurat_object_list) >= chunk_size) {
    chunk <- seurat_object_list[1:chunk_size]
    merged_chunk <- Reduce(merge, chunk)
    
    for (assay in assays) {
      merged_chunk <- JoinLayers(merged_chunk, assay = assay)
    }
    
    merged_data_list <- append(merged_data_list, list(merged_chunk))
    
    seurat_object_list <- seurat_object_list[-(1:chunk_size)]
    
    # Print progress info
    cat(paste("Merged", chunk_size, "objects. Remaining:", length(seurat_object_list), "\n"))
    
    # Clean up memory
    gc()
  }
  
  if (length(seurat_object_list) > 0) {
    remaining_merged <- Reduce(merge, seurat_object_list)
    
    for (assay in assays) {
      remaining_merged <- JoinLayers(remaining_merged, assay = assay)
    }
    
    merged_data_list <- append(merged_data_list, list(remaining_merged))
    
    # Print final progress info
    cat(paste("Merged", length(seurat_object_list), "objects. All objects processed.\n"))
    
    # Clean up memory
    gc()
  }
  
  final_merged_data <- merged_data_list[[1]]
  if (length(merged_data_list) > 1) {
    for (i in 2:length(merged_data_list)) {
      final_merged_data <- merge(final_merged_data, merged_data_list[[i]])
      
      for (assay in assays) {
        final_merged_data <- JoinLayers(final_merged_data, assay = assay)
      }
      
      # Print progress info
      cat(paste("Merged additional chunk. Remaining chunks:", length(merged_data_list) - i, "\n"))
      
      # Clean up memory
      gc()
    }
  }
  
  for (assay in assays) {
    final_merged_data <- JoinLayers(final_merged_data, assay = assay, layer = c('counts', 'data', 'scale.data'))
  }
  
  # Clean up memory
  gc()
  
  return(final_merged_data)
}

# Define the function to merge and join layers for multiple assays
merge_seurat_ATAC_libraries <- function(seurat_object_list, libraries = 5) {

  merged_data_list <- list()
  total_libraries <- length(seurat_object_list)
  
  while (length(seurat_object_list) >= libraries) {
    chunk <- seurat_object_list[1:libraries]
    merged_chunk <- Reduce(merge, chunk)
    
    merged_data_list <- append(merged_data_list, list(merged_chunk))
    
    seurat_object_list <- seurat_object_list[-(1:libraries)]
    
    # Print progress info
    cat(paste("Merged", libraries, "libraries. Remaining:", length(seurat_object_list), "\n"))
    
    # Clean up memory
    gc()
  }
  
  if (length(seurat_object_list) > 0) {
    remaining_merged <- Reduce(merge, seurat_object_list)

    merged_data_list <- append(merged_data_list, list(remaining_merged))
    
    # Print final progress info
    cat(paste("Merged", length(seurat_object_list), "libraries. All libraries processed.\n"))
    
    # Clean up memory
    gc()
  }
  
  final_merged_data <- merged_data_list[[1]]
  if (length(merged_data_list) > 1) {
    for (i in 2:length(merged_data_list)) {
      final_merged_data <- merge(final_merged_data, merged_data_list[[i]])

      # Print progress info
      cat(paste("Merged additional superlibrary. Remaining superlibraries:", length(merged_data_list) - i, "\n"))
      
      # Clean up memory
      gc()
    }
  }
  
  # Clean up memory
  gc()
  
  return(final_merged_data)
}


library(Seurat)

compare_overlaps <- function(libraries, top_n = 10) {
  get_overlap <- function(obj1, obj2) {
    overlap <- intersect(colnames(obj1), colnames(obj2))
    return(length(overlap))
  }
  
  n <- length(libraries)
  overlap_matrix <- matrix(0, n, n)
  rownames(overlap_matrix) <- colnames(overlap_matrix) <- names(libraries)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      overlap_matrix[i, j] <- get_overlap(libraries[[i]], libraries[[j]])
      overlap_matrix[j, i] <- overlap_matrix[i, j]
    }
  }
  
  overlap_df <- as.data.frame(as.table(overlap_matrix))
  colnames(overlap_df) <- c("seurat_object_1", "seurat_object_2", "common_cells")
  
  overlap_df <- overlap_df[overlap_df$seurat_object_1 != overlap_df$seurat_object_2 & overlap_df$common_cells > 0,]
  
  overlap_df$Pair <- apply(overlap_df[, c("seurat_object_1", "seurat_object_2")], 1, function(x) paste(sort(x), collapse = "_"))
  overlap_df <- overlap_df[!duplicated(overlap_df$Pair),]
  overlap_df$Pair <- NULL
  
  overlap_df <- overlap_df[order(-overlap_df$common_cells),]
  
  top_n_overlaps <- head(overlap_df, top_n)
  
  return(top_n_overlaps)
}
