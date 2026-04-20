#' Predict cell types for clusters using CELLiD
#'
#' Annotates clusters by comparing their average expression against the DISCO
#' CELLiD reference matrix using Spearman correlation.  By default the
#' computation runs on the DISCO server (identical to the
#' \href{https://www.immunesinglecell.com/tool/cellid}{online CELLiD tool}).
#' Set \code{local = TRUE} to download the reference once (~280 MB), cache it,
#' and score locally — useful for repeated runs or offline use.
#'
#' @param rna A Seurat object (clusters must already be defined), or a
#'   numeric matrix/data.frame where rows are genes and columns are clusters
#'   (pre-averaged expression values).
#' @param atlas Character. Reference atlas (default \code{"all"}).  Any atlas
#'   name shown on the website is accepted; \code{"all"} uses every atlas.
#' @param n_predict Integer. Number of top predictions per cluster (1 or 2,
#'   default 1).
#' @param local Logical. When \code{TRUE}, download the reference matrix and
#'   score locally (same algorithm as the server).  The reference is cached in
#'   \code{ref_path} so subsequent calls skip the download.
#'   Default \code{FALSE} (server-side scoring).
#' @param ref_path Character. Cache directory for the downloaded reference
#'   (used only when \code{local = TRUE}, default \code{"~/.disco_cache"}).
#' @param verbose Logical. Print progress (default \code{TRUE}).
#'
#' @return A \code{data.frame} with columns \code{cluster},
#'   \code{predicted_cell_type}, \code{atlas}, and \code{score}.
#'   When \code{n_predict = 2} a \code{rank} column is also included.
#' @export
#'
#' @examples
#' \dontrun{
#' # Server-side scoring (default)
#' predictions <- CELLiDCluster(seu, atlas = "blood")
#'
#' # Local scoring — downloads reference once (~280 MB), then cached
#' predictions <- CELLiDCluster(seu, atlas = "blood", local = TRUE)
#'
#' cluster_map <- setNames(predictions$predicted_cell_type, predictions$cluster)
#' seu$cell_type <- setNames(cluster_map[as.character(Seurat::Idents(seu))], colnames(seu))
#' }
CELLiDCluster <- function(rna, atlas = "all", n_predict = 1L,
                            local = FALSE, ref_path = "~/.disco_cache",
                            verbose = TRUE) {
  n_predict <- min(as.integer(n_predict), 2L)
  avg <- .cellid_avg_matrix(rna)

  if (local) {
    .cellid_local(avg, atlas, n_predict, ref_path, verbose)
  } else {
    .cellid_remote(avg, atlas, n_predict, verbose)
  }
}


# --------------------------------------------------------------------------
# Build average expression matrix from Seurat or pre-averaged matrix
# --------------------------------------------------------------------------
.cellid_avg_matrix <- function(rna) {
  if (inherits(rna, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE))
      cli::cli_abort("Package {.pkg Seurat} is required.")
    avg <- Seurat::AverageExpression(rna, assays = "RNA", layer = "counts",
                                     verbose = FALSE)$RNA
    # AverageExpression prefixes numeric cluster names with "g" (e.g. 0 -> g0); strip it
    colnames(avg) <- gsub("^g(\\d+)$", "\\1", colnames(avg))
    return(round(avg, 2))
  }
  if (is.matrix(rna) || is.data.frame(rna)) return(rna)
  cli::cli_abort("{.arg rna} must be a Seurat object or a numeric matrix.")
}


# --------------------------------------------------------------------------
# Remote (server-side) scoring — identical to the website
# --------------------------------------------------------------------------
.cellid_remote <- function(avg, atlas, n_predict, verbose) {
  lines      <- apply(avg, 1, function(row) paste(round(row, 2), collapse = "\t"))
  input_text <- paste(paste0(rownames(avg), "\t", lines), collapse = "\n")

  if (verbose)
    cli::cli_inform("Sending {nrow(avg)} genes x {ncol(avg)} cluster(s) to DISCO CELLiD ...")

  res <- .disco_post("/cellid/run", list(input = input_text, tissue = atlas))
  if (is.null(res) || length(res) == 0)
    cli::cli_abort("CELLiD returned no results. Check your input.")

  if (!is.data.frame(res))
    res <- as.data.frame(do.call(rbind, lapply(res, as.list)))

  .cellid_parse_remote(res, n_predict)
}

.cellid_parse_remote <- function(res, n_predict) {
  rows <- lapply(seq_len(nrow(res)), function(i) {
    r   <- res[i, ]
    out <- data.frame(
      cluster             = as.character(r$id),
      predicted_cell_type = sub("--.*$", "", as.character(r$ct1)),
      atlas               = sub("^[^-]*--", "", as.character(r$ct1)),
      score               = suppressWarnings(as.numeric(r$cor1)),
      rank                = 1L,
      stringsAsFactors    = FALSE
    )
    if (n_predict >= 2L) {
      out <- rbind(out, data.frame(
        cluster             = as.character(r$id),
        predicted_cell_type = sub("--.*$", "", as.character(r$ct2)),
        atlas               = sub("^[^-]*--", "", as.character(r$ct2)),
        score               = suppressWarnings(as.numeric(r$cor2)),
        rank                = 2L,
        stringsAsFactors    = FALSE
      ))
    }
    out
  })
  result <- do.call(rbind, rows)
  if (n_predict == 1L) result$rank <- NULL
  rownames(result) <- NULL
  result
}


# --------------------------------------------------------------------------
# Local scoring — same algorithm as CELLiD.R on the server
# Download ref.rds once, cache, then run 2-stage Spearman correlation
# --------------------------------------------------------------------------
.cellid_load_ref <- function(ref_path, verbose) {
  ref_path   <- path.expand(ref_path)
  dir.create(ref_path, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(ref_path, "cellid_ref.rds")

  if (file.exists(cache_file)) {
    if (verbose) cli::cli_inform("Loading cached reference from {.path {cache_file}}")
    return(readRDS(cache_file))
  }

  if (verbose) cli::cli_inform("Downloading CELLiD reference (~280 MB, cached after first use) ...")
  raw <- .disco_get_raw("/cellid/getReference")
  writeBin(raw, cache_file)
  if (verbose) cli::cli_inform("Reference saved to {.path {cache_file}}")
  readRDS(cache_file)
}

.cellid_local <- function(avg, atlas, n_predict, ref_path, verbose) {
  ref <- .cellid_load_ref(ref_path, verbose)

  # filter by atlas
  if (!is.null(atlas) && atlas != "all") {
    keep <- which(sub(".*--", "", colnames(ref)) == atlas)
    if (length(keep) == 0L)
      cli::cli_abort("No reference profiles found for atlas {.val {atlas}}.")
    ref <- ref[, keep, drop = FALSE]
  }

  # intersect genes (uppercase to match CELLiD.R)
  query_genes <- toupper(rownames(avg))
  rownames(avg) <- query_genes
  ref_genes   <- toupper(rownames(ref))
  rownames(ref) <- ref_genes
  common <- intersect(query_genes, ref_genes)

  if (length(common) < 500L)
    cli::cli_abort("Only {length(common)} genes overlap with the reference (need >= 500). Check gene names.")

  avg <- avg[common, , drop = FALSE]
  ref <- ref[common, , drop = FALSE]

  if (verbose)
    cli::cli_inform("Scoring {ncol(avg)} cluster(s) against {ncol(ref)} reference profiles ...")

  # Stage 1: Spearman correlation across all ref profiles → top-5 candidates
  # Stage 2: re-score using top-2000 variable genes among candidates
  results <- lapply(colnames(avg), function(cl) {
    q <- as.numeric(avg[, cl])

    stage1 <- apply(ref, 2, function(r) {
      stats::cor(r, q, method = "spearman", use = "complete.obs")
    })
    top5_idx <- order(stage1, decreasing = TRUE)[seq_len(min(5L, length(stage1)))]
    ref_top  <- ref[, top5_idx, drop = FALSE]

    # top-2000 variable genes among candidates
    row_vars <- apply(ref_top, 1, stats::var, na.rm = TRUE)
    var_idx  <- order(row_vars, decreasing = TRUE)[seq_len(min(2000L, nrow(ref_top)))]
    ref_var  <- ref_top[var_idx, , drop = FALSE]
    q_var    <- q[var_idx]

    stage2 <- apply(ref_var, 2, function(r) {
      stats::cor(r, q_var, method = "spearman", use = "complete.obs")
    })
    ord <- order(stage2, decreasing = TRUE)[seq_len(min(n_predict, length(stage2)))]

    data.frame(
      cluster             = cl,
      predicted_cell_type = sub("--.*$", "", names(stage2)[ord]),
      atlas               = sub("^[^-]*--", "", names(stage2)[ord]),
      score               = round(stage2[ord], 3),
      rank                = seq_along(ord),
      stringsAsFactors    = FALSE
    )
  })

  result <- do.call(rbind, results)
  if (n_predict == 1L) result$rank <- NULL
  rownames(result) <- NULL
  result
}


# Internal: download raw bytes (for binary files)
.disco_get_raw <- function(path) {
  req  <- httr2::request(.disco_url(path)) |>
    httr2::req_timeout(3600L) |>
    httr2::req_headers("User-Agent" = paste0("DISCOtoolkit/", utils::packageVersion("DISCOtoolkit")))
  resp <- httr2::req_perform(req)
  httr2::resp_body_raw(resp)
}
