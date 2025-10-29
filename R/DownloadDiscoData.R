#' Download data from DISCO database using the filtered metadata
#'
#' @param metadata The resulting metadata object produced by FilterDiscoMetadata function
#' @param output_dir The directory where the downloaded files will be stored
#' @return A character vector containing sample IDs that failed to download
#' @examples
#' # Download the data from project 'GSE174748' and store it in the 'disco_data' directory
#' metadata = FilterDiscoMetadata(
#'   project = "GSE174748"
#' )
#' failed = DownloadDiscoData(metadata, output_dir = "disco_data")
#' @importFrom progress progress_bar
#' @import Matrix
#' @export
DownloadDiscoData <- function(metadata, output_dir = "DISCOtmp") {

  tryCatch({
    if (!dir.exists(output_dir)) {
      message("Create output directory")
      dir.create(output_dir)
    }
  }, error = function(e){
    stop("The output directory cannot be created.")
  })

  samples = metadata$sample_metadata
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nrow(samples))
  cell_type_list = list()
  failed_list = c()

  message("Start downloading")

  for (i in 1:nrow(samples)) {
    sample_id <- samples$sample_id[i]
    project_id <- samples$project_id[i]
    output_file <- file.path(output_dir, paste0(sample_id, ".rds"))

    h5_url <- paste0(
      getOption("disco_url"),
      "download/getRawH5/",
      project_id, "/",
      sample_id
    )

    tmp <- tempfile(fileext = ".h5")

    success <- tryCatch({
      # Download
      utils::download.file(h5_url, tmp, mode = "wb", quiet = TRUE)
      # Read
      rna <- Seurat::Read10X_h5(tmp)

      # Cell type info
      cell <- tryCatch({
        read.csv(
          paste0(getOption("disco_url"), "toolkit/getCellTypeSample?sampleId=", sample_id),
          sep = "\t"
        )
      }, error = function(e) {
        warning(paste("Cell type info failed for", sample_id))
        NULL
      })

      if (!is.null(cell)) {
        rownames(cell) = cell$cell_id
        cell = cell[, c(3, 6)]
      }

      # Apply filters
      if (!is.null(metadata$filter$cell_type)) {
        cell = cell[which(cell$cell_type %in% metadata[["cell_type_metadata"]][["cell_type"]]), , drop = FALSE]
        if (metadata[["filter"]][["cell_type_confidence"]] == "high") {
          cell = cell[which(cell$cell_type_score >= 0.8), , drop = FALSE]
        } else if (metadata[["filter"]][["cell_type_confidence"]] == "medium") {
          cell = cell[which(cell$cell_type_score >= 0.6), , drop = FALSE]
        }
        rna = rna[, rownames(cell), drop = FALSE]
      }

      saveRDS(rna, output_file)
      cell_type_list[[i]] <- cell
      closeAllConnections()
      TRUE
    }, error = function(e) {
      warning(paste("Failed:", sample_id, "Reason:", e$message))
      failed_list <<- c(failed_list, sample_id)
      FALSE
    })

    pb$tick()
  }
  if (length(cell_type_list) > 0) {
    cell_type_list = do.call(rbind, cell_type_list)
    saveRDS(cell_type_list, file.path(output_dir, "cell_type.rds"))
  }

  message("Download complete")
  if (length(failed_list) > 0) {
    message("Some samples failed to download: ", paste(failed_list, collapse = ", "))
  }

  return(failed_list)
}
