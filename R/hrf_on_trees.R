#' @include hrf_functions.R
NULL

#' Prepare branches for HRF calculation
#' @description
#' Takes a 'phylo' object and extracts the relevant branches for each HRF calculation.
#'
#' @param tree a phylogeny of class 'phylo'
#'
#' @return an object of class tbl or data_frame
#' @export
#'
prep_branch_lengths <- function(tree) {

  edges_df <- dplyr::data_frame(id = 1:length(tree$edge[, 1]),
                      from = tree$edge[, 1],
                      to = tree$edge[, 2],
                      code = purrr::map2_chr(from, to, function(x, y) paste(x, y, sep = "-")),
                      this_branch = tree$edge.length)

  edges_df <- dplyr::mutate(edges_df,
                    descend = purrr::map(to, function(x) which(tree$edge[, 1] == x)),
                    ancest = purrr::map_int(from, function(x) match(x, tree$edge[,2])),
                    sib = purrr::map2_int(from, id, function(x, y) setdiff(which(tree$edge[, 1] == x), y)))

  edges_df <- dplyr::transmute(edges_df,
                        code,
                        from,
                        to,
                        this,
                        des = purrr::map(descend, function(x) tree$edge.length[x]),
                        anc = tree$edge.length[ancest],
                        sib = tree$edge.length[sib])

  return(edges_df)
}

#' Calculate HRF for all branches on a tree
#'
#' @description
#' It assumes:
#' no uncertainty on species tree ((A,B),C), one sample per species, branch lengths in coalescent units.
#'
#'
#' @param edges a tbl (or data_frame)
#' @param mutation numeric, population-wide mutation rate
#' @param mode string
#'
#' @return an object of class tbl or data_frame
#' @export
#' @seealso pr_hemiplasy(), pr_homoplasy(), single_hrf()
#'
tree_hrf <- function(edges, mutation = 0.01, mode = "standard") {
  model_avail <- c("standard", "minimal", "strict")
  if(!mode %in% model_avail) warning(paste("Requested model", mode,"in tree_hrf() not recognized."))

  if(is.null(edges$this_branch) | is.null(edges$descendants) | is.null(edges$ancestor) | is.null(edges$sibling)){
    warning("At least on variable is missing in the 'edges' data frame.")
  }
  mutate(edges,
         hrf = purrr::pmap_dbl(list(this_branch, descendants, ancestor, sibling), single_hrf, mutation_rate = mutation, model = mode))
}

#' Convert to 'phylo4' class
#'
#' @param tree
#' @param df
#'
#' @return an object of class phylo4
#' @export
#'
to_phylo4d <- function(tree, df) {
  if (!requireNamespace("phylobase", quietly = TRUE)) {
    stop("phylobase needed for this function to work. Please install it.",
         call. = FALSE)
  }

  earliest <- dplyr::filter(df, is.na(anc))
  edgeno <- dplyr::slice(earliest, 1)
  anc_node <- dplyr::data_frame(code = paste("0", edgeno$from, sep = "-"), to = edgeno$from, this = Inf, des = list(earliest$this), hrf = 0)
  new_df <- dplyr::bind_rows(anc_node, df)
  legacy_df <- as.data.frame(new_df)
  rownames(legacy_df) <- new_df$to

  phylobase::phylo4d(tree, all.data = legacy_df)
}

#' Convert to 'treedata' structure
#'
#' @param tree
#' @param df
#'
#' @return an object of class treedata
#' @export
#'
to_treedata <- function(tree, df) {
   ph <- ape::as.phylo(tree)
   dfout <- dplyr::mutate(df, node = to)
   new("treedata", phylo = ph, data = dfout)
}

