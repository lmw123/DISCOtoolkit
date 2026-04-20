# Internal package state
.disco_env <- new.env(parent = emptyenv())
.disco_env$base_url <- "https://immunesinglecell.com/disco_v3_api"
.disco_env$timeout  <- 60L

.onAttach <- function(libname, pkgname) {
  stats <- tryCatch({
    resp <- httr2::req_perform(
      httr2::request(paste0(.disco_env$base_url, "/getStatistics")) |>
        httr2::req_timeout(5L)
    )
    httr2::resp_body_json(resp)
  }, error = function(e) NULL)

  ver  <- as.character(utils::packageVersion("DISCOtoolkit"))
  link <- cli::style_hyperlink(
    cli::col_cyan("immunesinglecell.com"),
    "https://immunesinglecell.com"
  )

  title <- paste0(
    cli::style_bold(cli::col_blue("DISCO")),
    cli::col_silver("toolkit"),
    "  ",
    cli::col_silver(paste0("v", ver))
  )

  if (!is.null(stats)) {
    fmt <- function(n) formatC(as.numeric(n), format = "d", big.mark = ",")
    stats_line <- paste0(
      "  ", cli::col_yellow("\u25cf"), " ", cli::style_bold(fmt(stats$sampleCount)), " samples",
      "  ", cli::col_yellow("\u25cf"), " ", cli::style_bold(fmt(stats$cellCount)),   " cells",
      "  ", cli::col_yellow("\u25cf"), " ", cli::style_bold(fmt(stats$atlasCount)),  " atlases"
    )
    update_line <- if (!is.null(stats$lastUpdate)) {
      paste0("  ", cli::col_silver(paste0("last update: ", stats$lastUpdate)))
    } else ""
    packageStartupMessage(title, "\n", stats_line, "\n", update_line, "\n  ", link)
  } else {
    packageStartupMessage(title, "\n  ", link)
  }
}

# --------------------------------------------------------------------------
# Internal HTTP helpers
# --------------------------------------------------------------------------

.disco_url <- function(path) {
  paste0(.disco_env$base_url, "/", sub("^/", "", path))
}

.disco_get <- function(path, ...) {
  req <- httr2::request(.disco_url(path)) |>
    httr2::req_timeout(.disco_env$timeout) |>
    httr2::req_headers("User-Agent" = paste0("DISCOtoolkit/", utils::packageVersion("DISCOtoolkit"))) |>
    httr2::req_url_query(...)
  resp <- httr2::req_perform(req)
  .parse_response(resp)
}

.disco_post <- function(path, body) {
  req <- httr2::request(.disco_url(path)) |>
    httr2::req_timeout(.disco_env$timeout) |>
    httr2::req_headers("User-Agent" = paste0("DISCOtoolkit/", utils::packageVersion("DISCOtoolkit"))) |>
    httr2::req_body_json(body)
  resp <- httr2::req_perform(req)
  .parse_response(resp)
}

.parse_response <- function(resp) {
  ct <- httr2::resp_content_type(resp)
  if (grepl("json", ct)) {
    httr2::resp_body_json(resp, simplifyVector = TRUE)
  } else {
    httr2::resp_body_raw(resp)
  }
}

.as_df <- function(data) {
  if (is.null(data) || length(data) == 0) return(data.frame())
  df <- if (is.data.frame(data)) data else as.data.frame(data)
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.data.frame(col) || is.list(col)) {
      df[[nm]] <- vapply(seq_len(nrow(df)), function(i) {
        x <- if (is.data.frame(col)) col[i, , drop = FALSE] else col[[i]]
        if (is.null(x) || length(x) == 0) NA_character_
        else jsonlite::toJSON(x, auto_unbox = TRUE, na = "null")
      }, character(1L))
    }
  }
  id_col <- intersect(c("sample_id", "sampleId"), names(df))[1]
  if (!is.na(id_col) && !anyDuplicated(df[[id_col]])) {
    rownames(df) <- df[[id_col]]
  }
  df
}

# --------------------------------------------------------------------------
# Internal helper — replaces missing-value NULL with default
# --------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
