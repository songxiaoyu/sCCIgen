#' Convert a data.frame to config.yml
#'
#' @param x a data.frame with two columns (parameters and value)
#' @param parameter_file name for the config file. Default = config.yml
#'
#' @returns a config.yml file
#' @export
#'
#' @examples
#' df <- data.frame(
#'     parameters = c("expression_data_file", "num_simulated_cells"),
#'     value = c("fake1_expr.Rdata", 10000))
#' df_to_config(df)
df_to_config <- function(x, parameter_file = "config.yml") {
  # add a new column to input df
  x$yaml <- paste0(x$parameters, ": ", x$value)

  # collapse parameters in a single line
  x_yaml <- paste0(x$yaml, collapse = "\n    ")
  x_yaml <- paste0("default:\n    ", x_yaml)

  # export parameters to config yml
  cat(x_yaml, file = parameter_file, sep = "")
}


#' Convert sCCIgen output to a Giotto object
#'
#' @param counts_file path to .tsv file with gene x cell counts matrix.
#' @param metadata_file path to .tsv file with metadata. The file should
#' contain the following columns: Cell, annotation, x.loc, y.loc, and region.
#'
#' @returns a Giotto object
#' @export
#'
#' @examples
#' giotto_object <- sCCIgen_to_Giotto()
sCCIgen_to_Giotto <- function(counts_file, metadata_file) {

  x_expression <- read.delim(counts_file,
                             row.names = 1)

  x_meta <- read.delim(metadata_file)

  x_spatlocs <- x_meta[,c("Cell", "x.loc", "y.loc")]
  colnames(x_spatlocs) <- c("cell_ID", "sdimx", "sdimy")

  x_meta <- x_meta[,c("Cell", "annotation", "region")]
  colnames(x_meta)[1] <- "cell_ID"

  x <- Giotto::createGiottoObject(expression = x_expression,
                                  spatial_locs = x_spatlocs)

  x <- Giotto::addCellMetadata(x,
                               new_metadata = x_meta[-1])

  return(x)
}


#' Convert sCCIgen output to a Seurat object
#'
#' @param counts_file path to .tsv file with gene x cell counts matrix.
#' @param metadata_file path to .tsv file with metadata. The file should
#' contain the following columns: Cell, annotation, x.loc, y.loc, and region.
#'
#' @returns a Seurat object
#' @export
#'
#' @examples
#' seurat_object <- sCCIgen_to_Seurat()
sCCIgen_to_Seurat <- function(counts_file, metadata_file) {

  x_expression <- read.delim(counts_file,
                             row.names = 1)

  x_meta <- read.delim(metadata_file)

  x_spatlocs <- x_meta[,c("x.loc", "y.loc")]
  rownames(x_spatlocs) <- x_meta$Cell
  colnames(x_spatlocs) <- c("x", "y")

  x_metadata <- x_meta[,c("annotation", "region")]
  rownames(x_metadata) <- x_meta$Cell

  x <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(as.matrix(x_expression), sparse = TRUE),
    assay = "Spatial")

  x <- Seurat::AddMetaData(x,
                           metadata = x_metadata,
                           col.name = names(x_metadata))

  x@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = x_spatlocs
  )

  return(x)
}

#' Convert sCCIgen output to a SpatialExperiment object
#'
#' @param counts_file path to .tsv file with gene x cell counts matrix.
#' @param metadata_file path to .tsv file with metadata. The file should
#' contain the following columns: Cell, annotation, x.loc, y.loc, and region.
#'
#' @returns a SpatialExperiment object
#' @export
#'
#' @examples
#' spe_object <- sCCIgen_to_SpatialExperiment()
sCCIgen_to_SpatialExperiment <- function(counts_file, metadata_file) {

  x_expression <- read.delim(counts_file,
                             row.names = 1)

  x_meta <- read.delim(metadata_file)

  x_spatlocs <- x_meta[,c("x.loc", "y.loc")]
  rownames(x_spatlocs) <- x_meta$Cell
  colnames(x_spatlocs) <- c("x", "y")

  x_metadata <- x_meta[,c("annotation", "region")]
  rownames(x_metadata) <- x_meta$Cell

  x <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = Matrix::Matrix(as.matrix(x_expression), sparse = TRUE)),
    colData = x_metadata,
    spatialCoords = as.matrix(x_spatlocs))

  return(x)
}
