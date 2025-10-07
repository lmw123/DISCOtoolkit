#' Download data from DISCO database using the filtered metadata
#'
#' @param metadata The resulting metadata object produced by FilterDiscoMetadata function
#' @param output_dir The directory where the downloaded files will be stored
#' @return NULL
#' @examples
#' # Download the data from project 'GSE174748' and store it in the 'disco_data' directory
#' metadata = FilterDiscoMetadata(
#'   project = "GSE174748"
#' )
#' DownloadDiscoData(metadata, output_dir = "disco_data")
#' @importFrom progress progress_bar
#' @import Matrix
#' @export
DownloadDiscoData <- function(metadata, output_dir = "DISCOtmp") {

  tryCatch({
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }
  }, error = function(e){
    stop("The output directory cannot be created.")
  })


  samples = metadata$sample_metadata

  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nrow(samples))

  cell_type_list = list()
  message("Start downloading")
  for (i in 1:nrow(samples)) {
    output_file = paste0(output_dir, "/", samples$sample_id[i], ".h5")

    h5_url      <- paste0(
      getOption("disco_url"),
      "download/getRawH5/",
      samples$project_id[i], "/",
      samples$sample_id[i]
    )
    local_h5    <- file.path(output_dir, paste0(samples$sample_id[i], ".h5"))

    # download .h5 if it doesn't already exist

    download.file(h5_url, destfile = local_h5, mode = "wb")

    rna <- Seurat::Read10X_h5(local_h5)

    # 删除刚刚下载的 .h5 文件
    if (file.exists(local_h5)) {
      file.remove(local_h5)
    }

    cell = read.csv(paste0(getOption("disco_url"), "toolkit/getCellTypeSample?sampleId=", samples$sample_id[i]), sep = "\t")
    rownames(cell) = cell$cell_id
    cell = cell[,c(3,6)]
    if (!is.null(metadata$filter$cell_type)) {
      cell = cell[which(cell$cell_type %in% metadata[["cell_type_metadata"]][["cell_type"]]),,drop=F]

      if (metadata[["filter"]][["cell_type_confidence"]] == "high") {
        cell = cell[which(cell$cell_type_score >= 0.8),, drop=F]
      } else if  (metadata[["filter"]][["cell_type_confidence"]] == "medium") {
        cell = cell[which(cell$cell_type_score >= 0.6),, drop=F]
      }

      rna = rna[,rownames(cell),drop=F]
    }

    saveRDS(rna, output_file)
    cell_type_list[[i]] = cell
    closeAllConnections()
    pb$tick()
  }

  cell_type_list = do.call(rbind, cell_type_list)
  saveRDS(cell_type_list, paste0(output_dir, "/cell_type.rds"))
  message("Download complete")
}

