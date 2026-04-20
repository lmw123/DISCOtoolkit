#' Filter and retrieve DISCO sample metadata
#'
#' @param sample_id Character vector. Filter by sample ID.
#' @param project_id Character vector. Filter by project ID.
#' @param tissue Character vector. Filter by tissue (e.g. \code{"blood"}, \code{"lung"}).
#' @param disease Character vector. Filter by disease (e.g. \code{"COVID-19"}, \code{"healthy"}).
#' @param platform Character vector. Sequencing platform.
#' @param sample_type Character vector. Sample type.
#' @param limit Integer. Maximum rows to return. \code{NULL} returns all.
#'
#' @return A \code{data.frame} of sample metadata.
#' @export
#'
#' @examples
#' \dontrun{
#' meta <- FilterDiscoMetadata(tissue = "blood", disease = "COVID-19")
#' meta <- FilterDiscoMetadata(tissue = c("blood", "liver"), limit = 500)
#' }
FilterDiscoMetadata <- function(
    sample_id   = NULL,
    project_id  = NULL,
    tissue      = NULL,
    disease     = NULL,
    platform    = NULL,
    sample_type = NULL,
    limit       = NULL) {

  body <- list(
    tissue     = as.list(tissue      %||% character()),
    disease    = as.list(disease     %||% character()),
    platform   = as.list(platform    %||% character()),
    sampleType = as.list(sample_type %||% character()),
    projectId  = as.list(project_id  %||% character())
  )
  data <- .disco_post("/repository/get_all_metadata", body)
  df   <- .as_df(data)

  if (!is.null(sample_id)) df <- df[df$sample_id %in% sample_id, ]
  if (!is.null(limit) && nrow(df) > limit) df <- df[seq_len(limit), ]
  df
}


#' List available filter values for a metadata field
#'
#' @param field Character. One of \code{"tissue"}, \code{"disease"},
#'   \code{"platform"}, \code{"sample_type"}, \code{"project_id"}.
#'
#' @return A \code{data.frame} with columns \code{value} and \code{count},
#'   sorted by count descending.
#' @export
#'
#' @examples
#' \dontrun{
#' ListMetadataItem("tissue")
#' ListMetadataItem("disease")
#' }
ListMetadataItem <- function(field) {
  data <- .disco_post("/repository/get_metadata_options", list(type = field))
  if (is.list(data) && !is.data.frame(data)) {
    data.frame(
      value = sapply(data, `[[`, 1),
      count = as.integer(sapply(data, `[[`, 2)),
      stringsAsFactors = FALSE
    )
  } else {
    .as_df(data)
  }
}


# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x


.disco_onto_env <- new.env(parent = emptyenv())
.disco_onto_env$cache <- NULL

.get_disco_ontology <- function() {
  if (!is.null(.disco_onto_env$cache)) return(.disco_onto_env$cache)
  data <- .disco_get("/cell_type/getCellOntology")
  onto <- if (is.data.frame(data)) data else as.data.frame(data)
  .disco_onto_env$cache <- onto
  onto
}

.collect_children <- function(ct, onto) {
  out      <- ct
  children <- onto$cell_name[onto$parent == ct]
  for (ch in children) out <- c(out, .collect_children(ch, onto))
  unique(out)
}
