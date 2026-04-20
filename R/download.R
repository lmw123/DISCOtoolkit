#' Download DISCO samples
#'
#' Downloads H5 expression files for all samples in a metadata table (as
#' returned by \code{\link{FilterDiscoMetadata}}).  Already-downloaded files
#' are skipped.
#'
#' @param metadata A \code{data.frame} with a \code{sample_id} column.
#'   Typically the output of \code{FilterDiscoMetadata()}.
#' @param output_dir Character. Directory to save H5 files (default
#'   \code{"DISCOtmp"}).
#'
#' @return A \code{data.frame} with columns \code{sample_id}, \code{status}
#'   (\code{"ok"} or \code{"failed"}), and \code{path} (local file path, or
#'   \code{NA} on failure).
#' @export
#'
#' @examples
#' \dontrun{
#' meta <- FilterDiscoMetadata(tissue = "blood", disease = "COVID-19", limit = 5)
#' res  <- DownloadDiscoData(meta, output_dir = "DISCOtmp")
#' res[res$status == "ok", ]
#' }
DownloadDiscoData <- function(metadata, output_dir = "DISCOtmp") {
  if (!is.data.frame(metadata) || nrow(metadata) == 0L)
    cli::cli_abort("{.arg metadata} must be a non-empty data.frame.")

  sample_col <- intersect(c("sample_id", "sampleId"), names(metadata))[1]
  if (is.na(sample_col))
    cli::cli_abort("metadata must contain a {.field sample_id} column.")

  sample_ids <- unique(metadata[[sample_col]])
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  cli::cli_inform("Downloading {length(sample_ids)} sample(s) into {.path {output_dir}} ...")

  pb <- progress::progress_bar$new(
    format = "[:bar] :current/:total :percent  ETA :eta",
    total  = length(sample_ids), clear = FALSE
  )

  results <- lapply(sample_ids, function(sid) {
    pb$tick()
    h5_path <- file.path(output_dir, paste0(sid, ".h5"))
    if (file.exists(h5_path)) {
      return(data.frame(sample_id = sid, status = "ok", path = h5_path,
                        stringsAsFactors = FALSE))
    }
    tryCatch({
      raw <- .disco_get(paste0("/sample/downloadH5/", sid))
      writeBin(raw, h5_path)
      data.frame(sample_id = sid, status = "ok", path = h5_path,
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      cli::cli_warn("Failed to download {.val {sid}}: {conditionMessage(e)}")
      data.frame(sample_id = sid, status = "failed", path = NA_character_,
                 stringsAsFactors = FALSE)
    })
  })

  out <- do.call(rbind, results)
  n_ok   <- sum(out$status == "ok")
  n_fail <- sum(out$status == "failed")
  cli::cli_inform("Done. {n_ok} succeeded, {n_fail} failed.")
  out
}
