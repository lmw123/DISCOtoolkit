
#' Get metadata from DISCO and apply filters on it
#'
#' @param sample sample ID(s) of interest, can be a string or vector
#' @param project project ID(s) of interest, can be a string or vector
#' @param tissue tissue(s) of interest, can be a string or vector
#' @param disease disease(s) of interest, can be a string or vector
#' @param platform platform(s) of interest, can be a string or vector
#' @param sample_type sample type(s) of interest, can be a string or vector
#' @param cell_type cell type(s) of interest, can be a string or vector
#' @param cell_type_confidence  Filter the results based on the confidence of cell type prediction, which can be categorized as high, medium, or all
#' @param include_cell_type_children Whether the cell type is an exact match or includes its child elements
#' @param min_cell_per_sample Only samples with a cell count greater than the minimum cell count per sample will be retained.
#' @export
#' @importFrom jsonlite fromJSON
FilterDiscoMetadata <- function(sample_id = NULL,
                                project_id = NULL,
                                tissue = NULL,
                                disease = NULL,
                                platform = NULL,
                                sample_type = NULL,
                                cell_type = NULL,
                                cell_type_confidence = "medium",
                                include_cell_type_children = T,
                                min_cell_per_sample = 100) {

  filter_data = list(
    sample_metadata = NULL,
    cell_type_metadata = NULL,
    sample_count = NULL,
    cell_count = NULL,
    filter =  list(
      sample_id = sample_id, project_id = project_id, tissue = tissue, disease = disease, platform = platform,
      sample_type = sample_type, cell_type = cell_type, cell_type_confidence = cell_type_confidence,
      include_cell_type_children = include_cell_type_children, min_cell_per_sample = min_cell_per_sample
    )
  )

  message("Fetching sample metadata")
  metadata = read.csv(paste0(getOption("disco_url"), "toolkit/getSampleMetadata"), sep = "\t")
  rownames(metadata) = metadata$sample_id

  message("Filtering sample")
  if (is.null(sample_id) == F) {
    metadata = metadata[which(metadata$sample_id %in% sample_id),]
  }

  if (is.null(project_id) == F) {
    metadata = metadata[which(metadata$project_id %in% project_id),]
  }

  if (is.null(tissue) == F) {
    metadata = metadata[which(metadata$tissue %in% tissue),]
  }

  if (is.null(platform) == F) {
    metadata = metadata[which(metadata$platform %in% platform),]
  }

  if (is.null(disease) == F) {
    metadata = metadata[which(metadata$disease %in% disease),]
  }

  if (is.null(sample_type) == F) {
    metadata = metadata[which(metadata$sample_type %in% sample_type),]
  }

  if (nrow(metadata) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }

  message("Fetching cell type information")
  sample_ct_info = read.csv(paste0(getOption("disco_url"), "toolkit/getCellTypeSummary"), sep = "\t")


  if (is.null(cell_type)) {

    sample_ct_info = sample_ct_info[which(sample_ct_info$sample_id %in% metadata$sample_id),]
    sample_cell_count = aggregate(sample_ct_info$cell_number, by=list(sample=sample_ct_info$sample_id), FUN=sum)
    rownames(sample_cell_count) = sample_cell_count$sample
    metadata$cell_number = sample_cell_count[rownames(metadata),"x"]

    metadata = metadata[which(metadata$cell_number > min_cell_per_sample),]
    if (nrow(metadata) == 0) {
      stop("Sorry, no samples passed the applied filters.")
    }
    sample_ct_info = sample_ct_info[which(sample_ct_info$sample_id %in% metadata$sample_id),]

    metadata = metadata[,apply(metadata, 2, function(i){sum(is.na(i))}) < nrow(metadata)]


    filter_data[["sample_metadata"]] = metadata
    filter_data[["cell_type_metadata"]] = sample_ct_info
    filter_data[["sample_count"]] = nrow(metadata)
    filter_data[["cell_count"]] = sum(metadata$cell_number)
    message(paste0(filter_data[["sample_count"]], " samples and ", filter_data[["cell_count"]]," cells were found"))
    return(filter_data)
  }

  if (include_cell_type_children) {
    cell_type = GetCellTypeChildren(cell_type)
  }

  sample_ct_info = sample_ct_info[which(sample_ct_info$cell_type %in% cell_type),]
  sample_ct_info = sample_ct_info[which(sample_ct_info$sample_id %in% metadata$sample_id),]

  if (!(cell_type_confidence %in% c("high", "medium", "all"))) {
    stop("cell_type_confidence can only be high, medium, or all")
  }

  if (cell_type_confidence == "high") {
    sample_ct_info = sample_ct_info[which(sample_ct_info$cell_type_score >= 0.8),]
  } else if  (cell_type_confidence == "medium") {
    sample_ct_info = sample_ct_info[which(sample_ct_info$cell_type_score >= 0.6),]
  }


  if (nrow(sample_ct_info) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }

  metadata = metadata[which(metadata$sample_id %in% sample_ct_info$sample_id),]
  sample_cell_count = aggregate(sample_ct_info$cell_number, by=list(sample=sample_ct_info$sample_id), FUN=sum)

  rownames(sample_cell_count) = sample_cell_count$sample
  metadata$cell_number = sample_cell_count[rownames(metadata),"x"]

  metadata = metadata[which(metadata$cell_number > min_cell_per_sample),]
  if (nrow(metadata) == 0) {
    stop("Sorry, no samples passed the applied filters.")
  }
  sample_ct_info = sample_ct_info[which(sample_ct_info$sample_id %in% metadata$sample_id),]

  filter_data[["sample_metadata"]] = metadata
  filter_data[["cell_type_metadata"]] = sample_ct_info
  filter_data[["sample_count"]] = nrow(metadata)
  filter_data[["cell_count"]] = sum(metadata$cell_number)

  message(paste0(filter_data[["sample_count"]], " samples and ", filter_data[["cell_count"]]," cells were found"))

  return(filter_data)
}


#' Check children of input cell type
#'
#' @param cell_type The cell type of the input, which may be a string or a list of string
#' @param cell_ontology The cell type ontology. If not specified, the cell type ontology will be retrieved from the DISCO database.
#' @return List of cell type
#' @examples
#' # Get children of B cell
#' GetCellTypeChildren(cell_type = c("B cell"))
#'
#' # Get children of B cell and Macrophage
#' GetCellTypeChildren(cell_type = c("B cell", "Macrophage"))
#' @export
GetCellTypeChildren <- function(cell_type, cell_ontology = NULL) {
  if(is.null(cell_ontology)) {
    cell_ontology = GetJson(
      url = paste0(getOption("disco_url"), "toolkit/getCellOntology"),
      info.msg = "Fetching ontology from DISCO database",
      error.msg = "Failed to retrieve ontology Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance."
    )
  }

  children = c(cell_type)

  while(length(cell_type) > 0) {
    children = append(children, cell_ontology$cell_name[which(cell_ontology$parent %in% cell_type)])
    cell_type = cell_ontology$cell_name[which(cell_ontology$parent %in% cell_type)]
  }

  children = unique(children)
  return(children)
}



#' Check children of input cell type
#'
#' @param term A partial or complete name of a cell type
#' @param cell_ontology The cell type ontology. If not specified, the cell type ontology will be retrieved from the DISCO database.
#' @return List of cell type
#' @examples
#' # Retrieve cell types whose names contain the term 'blast'
#' GetCellTypeChildren(cell_type = c("blast"))
#' @export
FindCellType <- function(term = "", cell_ontology = NULL) {
  if(is.null(cell_ontology)) {
    cell_ontology = GetJson(
      url = paste0(getOption("disco_url"), "toolkit/getCellOntology"),
      info.msg = "Retrieving ontology from DISCO database",
      error.msg = "Failed to retrieve ontology Please try again. If the issue persists, please contact us at li_mengwei@immunol.a-star.edu.sg for assistance."
    )
  }
  cell_type = cell_ontology$cell_name
  cell_type = grep(term, cell_type, ignore.case = T, value = T)
  return(cell_type)
}


#' List all potential values of a metadata field
#'
#' @export
ListMetadataItem <- function(field){
  message("Fetching sample metadata")
  metadata = read.csv(paste0(getOption("disco_url"), "toolkit/getSampleMetadata"), sep = "\t")
  if (field %in% colnames(metadata)) {
    return(unique(metadata[,field]))
  } else {
    stop(paste0("DISCO data don't contain ", field, " field"))
  }
}


