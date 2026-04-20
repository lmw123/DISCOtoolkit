#' Gene set enrichment against DISCO signatures
#'
#' Tests a gene list against two complementary DISCO references:
#'
#' \describe{
#'   \item{\strong{Cell type}}{Hypergeometric test against cell-type DEG
#'     markers derived from DISCO atlases.  Answers: "which cell type does
#'     my gene list represent?"}
#'   \item{\strong{Phenotype}}{Fisher exact test against phenotype DEG
#'     signatures (e.g. "COVID-19 vs control for T cell"), identical to the
#'     \href{https://www.immunesinglecell.com/tool/sc_enrichment}{online
#'     sc_enrichment tool}.  Answers: "which disease/condition is my gene
#'     list associated with?"}
#' }
#'
#' @param input Either a character vector of gene symbols, or a
#'   \code{data.frame} with a \code{gene} column and an optional \code{logFC}
#'   column.  When \code{logFC} is present only genes with \code{logFC > 0}
#'   are used for cell-type enrichment; the full signed values are passed to
#'   the phenotype test (matching website behaviour).
#' @param min_markers Integer. Minimum markers a cell type must have to be
#'   included in the cell-type test (default 5).
#' @param n_universe Integer. Gene universe size for the hypergeometric test
#'   (default 20000).
#'
#' @return A \code{data.frame} sorted by p-value with columns:
#'   \code{type} (\code{"Cell Type"} or \code{"Phenotype"}), \code{name},
#'   \code{overlap}, \code{overlap_genes}, \code{n_geneset}, \code{pval},
#'   \code{pval_adj}, \code{odds_ratio} (phenotype only, \code{NA} for
#'   cell-type rows).
#' @export
#'
#' @examples
#' \dontrun{
#' # Gene list only
#' res <- CELLiDEnrichment(c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"))
#' head(res)
#' #         type      name overlap ...  pval pval_adj odds_ratio
#' # 1  Cell Type      Treg       5 ...  ...      ...         NA
#' # 2  Phenotype COVID-19…       3 ...  ...      ...        8.2
#'
#' # DEG table with logFC
#' res <- CELLiDEnrichment(deg_df)   # deg_df has columns: gene, logFC
#' }
CELLiDEnrichment <- function(input,
                              min_markers = 5L,
                              n_universe  = 20000L) {
  # --- parse input ---
  if (is.character(input)) {
    genes_up   <- unique(input)
    genes_full <- unique(input)   # no logFC: same set for both tests
    has_fc     <- FALSE
  } else if (is.data.frame(input)) {
    gene_col <- if ("gene" %in% names(input)) "gene" else names(input)[1]
    fc_col   <- intersect(c("logFC", "logfc", "log2FC", "log2fc", "LFC"), names(input))
    has_fc   <- length(fc_col) > 0L
    if (has_fc) {
      genes_up   <- unique(as.character(input[[gene_col]][input[[fc_col[1]]] > 0]))
      genes_full <- input[, c(gene_col, fc_col[1])]   # keep full df for phenotype
    } else {
      genes_up   <- unique(as.character(input[[gene_col]]))
      genes_full <- genes_up
    }
  } else {
    cli::cli_abort("{.arg input} must be a character vector or data.frame.")
  }

  ct_result  <- .enrichment_cell_type(genes_up, min_markers, n_universe)
  phe_result <- .enrichment_phenotype(genes_full, has_fc)

  # Unify into a single data.frame
  ct_df <- if (nrow(ct_result) > 0L) {
    data.frame(
      type          = "Cell Type",
      name          = ct_result$cell_type,
      overlap       = ct_result$overlap,
      overlap_genes = ct_result$overlap_genes,
      n_geneset     = ct_result$n_markers,
      pval          = ct_result$pval,
      pval_adj      = ct_result$pval_adj,
      odds_ratio    = NA_real_,
      stringsAsFactors = FALSE
    )
  } else data.frame()

  phe_df <- if (nrow(phe_result) > 0L) {
    data.frame(
      type          = "Phenotype",
      name          = phe_result$name,
      overlap       = phe_result$overlap,
      overlap_genes = phe_result$overlap_genes,
      n_geneset     = phe_result$n_geneset,
      pval          = phe_result$pval,
      pval_adj      = stats::p.adjust(phe_result$pval, method = "BH"),
      odds_ratio    = phe_result$odds_ratio,
      stringsAsFactors = FALSE
    )
  } else data.frame()

  out <- rbind(ct_df, phe_df)
  out[order(out$pval), ]
}


# --------------------------------------------------------------------------
# Cell type enrichment — hypergeometric vs Cell Type DEG markers
# --------------------------------------------------------------------------
.enrichment_cell_type <- function(genes, min_markers, n_universe) {
  if (length(genes) == 0L) {
    cli::cli_warn("No genes for cell-type enrichment after filtering.")
    return(data.frame())
  }

  ct_data   <- .disco_get("/toolkit/getCellTypeDEGMarkers")
  if (!is.data.frame(ct_data)) ct_data <- as.data.frame(ct_data)

  gene_sets <- lapply(seq_len(nrow(ct_data)), function(i) {
    raw <- ct_data[["genes"]][i]
    if (is.null(raw) || is.na(raw) || raw == "") return(character())
    trimws(strsplit(raw, ";")[[1]])
  })
  names(gene_sets) <- ct_data[["cell_type"]]
  gene_sets <- Filter(function(g) length(g) >= min_markers, gene_sets)

  if (length(gene_sets) == 0L)
    cli::cli_abort("No cell types with >= {min_markers} markers found.")

  k    <- length(genes)
  rows <- lapply(names(gene_sets), function(ct) {
    ref_genes <- gene_sets[[ct]]
    overlap   <- intersect(genes, ref_genes)
    m         <- length(ref_genes)
    x         <- length(overlap)
    pval      <- stats::phyper(x - 1L, m, n_universe - m, k, lower.tail = FALSE)
    data.frame(
      cell_type     = ct,
      overlap       = x,
      overlap_genes = paste(overlap, collapse = ";"),
      n_markers     = m,
      pval          = pval,
      stringsAsFactors = FALSE
    )
  })

  out          <- do.call(rbind, rows)
  out$pval_adj <- stats::p.adjust(out$pval, method = "BH")
  out[order(out$pval_adj, out$pval), ]
}


# --------------------------------------------------------------------------
# Phenotype enrichment — server-side Fisher test (same as website)
# --------------------------------------------------------------------------
.enrichment_phenotype <- function(genes_full, has_fc) {
  # Serialise input for /cellid/run_gene
  if (has_fc) {
    # data.frame with gene + logFC columns
    gene_col <- names(genes_full)[1]
    fc_col   <- names(genes_full)[2]
    lines    <- paste(toupper(genes_full[[gene_col]]),
                      genes_full[[fc_col]], sep = "\t")
    input_text <- paste(lines, collapse = "\n")
  } else {
    input_text <- paste(toupper(genes_full), collapse = "\n")
  }

  res <- tryCatch(
    .disco_post("/cellid/run_gene", list(input = input_text)),
    error = function(e) {
      cli::cli_warn("Phenotype enrichment failed: {conditionMessage(e)}")
      NULL
    }
  )

  if (is.null(res) || length(res) == 0) return(data.frame())

  # run_gene returns list of rows or a data.frame
  if (!is.data.frame(res))
    res <- as.data.frame(do.call(rbind, lapply(res, as.list)))

  # Normalise column names from the TSV output
  # expected: pval, or, name, gene, background, overlap, geneset
  colnames(res) <- tolower(colnames(res))
  rename_map <- c(or = "odds_ratio", gene = "overlap_genes",
                  background = "n_background", geneset = "n_geneset")
  for (old in names(rename_map)) {
    if (old %in% names(res)) names(res)[names(res) == old] <- rename_map[old]
  }

  for (col in c("pval", "odds_ratio", "overlap", "n_geneset", "n_background")) {
    if (col %in% names(res)) res[[col]] <- suppressWarnings(as.numeric(res[[col]]))
  }

  cols_order <- intersect(
    c("name", "pval", "odds_ratio", "overlap", "n_geneset", "n_background", "overlap_genes"),
    names(res)
  )
  res[order(res$pval), cols_order, drop = FALSE]
}
