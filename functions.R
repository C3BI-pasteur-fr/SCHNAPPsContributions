# modified version of Read10X_h5 from Seurat to include the sample name in the barcode
#' Read HDF5 File
#' 
#' Reads data from an HDF5 file and returns a sparse matrix or a list of sparse matrices.
#' modified version of Read10X_h5 from Seurat to include the sample name in the barcode
#' 
#' @param filename The path to the HDF5 file.
#' @param use.names Logical, indicating whether to use feature names or IDs.
#' @param unique.features Logical, indicating whether to make feature names unique.
#' @param sampleName A character string specifying the sample name to be appended to column names.
#' 
#' @return A sparse matrix or a list of sparse matrices containing the data from the HDF5 file.
#' 
#' @export readH5
#' 
readH5 <- function(filename, use.names = TRUE, unique.features = TRUE, sampleName = "1") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    if (use.names) {
      feature_slot <- "features/name"
    }
    else {
      feature_slot <- "features/id"
    }
  }  else {
    if (use.names) {
      feature_slot <- "gene_names"
    }
    else {
      feature_slot <- "genes"
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    features <- infile[[paste0(genome, "/", feature_slot)]][]
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1, p = indptr[],
      x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, USE.NAMES = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  
  if (length(x = output) == 1) {
    if (class(output[[genome]]) == "list") {
      for (le in 1:length(output[[genome]])) {
        colnames(output[[genome]][[le]]) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output[[genome]][[le]])), "-", sampleName)
      }
    }else{
      colnames(output[[genome]]) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output[[genome]])), "-", sampleName)
    }
    return(output[[genome]])
  }
  else {
    # colnames(output) <- paste0(sub("(.*)-(.*)", "\\1", colnames(output)), "-", sampleName)
    return(output)
  }
}

