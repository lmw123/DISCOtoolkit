#' Load DISCO example datasets
#'
#' Returns one of the built-in example gene lists used on the DISCO website.
#'
#' @param dataset Character. Which example to load:
#'   \describe{
#'     \item{\code{"deg"}}{A \code{data.frame} with columns \code{gene} and
#'       \code{logFC} — pancreatic acinar cell DEGs (type 2 diabetes atlas).
#'       Use with \code{CELLiDEnrichment()} for weighted enrichment.}
#'     \item{\code{"genes"}}{A character vector of ISG/interferon-stimulated
#'       genes.  Use with \code{CELLiDEnrichment()} for standard enrichment.}
#'   }
#'
#' @return A \code{data.frame} or character vector depending on \code{dataset}.
#' @export
#'
#' @examples
#' deg   <- disco_example("deg")
#' genes <- disco_example("genes")
#'
#' res <- CELLiDEnrichment(deg)
#' res <- CELLiDEnrichment(genes)
disco_example <- function(dataset = c("deg", "genes")) {
  dataset <- match.arg(dataset)

  if (dataset == "deg") {
    raw <- "PRSS1\t5.24\nCTRB2\t5.18\nCELA3A\t4.84\nCTRB1\t4.72\nREG1B\t4.70\nCLPS\t4.65\nCPB1\t4.51\nCPA1\t4.35\nPLA2G1B\t4.31\nREG3A\t4.27\nREG1A\t4.25\nCTRC\t4.21\nPNLIP\t4.06\nCELA3B\t4.00\nPRSS3\t3.85\nSYCN\t3.78\nCELA2A\t3.54\nCPA2\t3.50\nAMY2A\t3.45\nREG3G\t3.20"
    con <- textConnection(raw)
    on.exit(close(con))
    df  <- utils::read.table(con, sep = "\t", header = FALSE,
                             col.names = c("gene", "logFC"),
                             stringsAsFactors = FALSE)
    return(df)
  }

  # dataset == "genes"
  c("IFI6","IFITM1","ISG15","IFITM3","APOBEC3A","IL1R2","LY6E","CCL2",
    "S100A8","CD163","AREG","RNASE1","RNASE2","RETN","MT2A","THBS1",
    "CLU","LMNA","ISG20","S100A12","IFI44L","SIGLEC1","RNF213","CTSD","FOLR3")
}
