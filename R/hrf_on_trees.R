#' @include hrf_functions.R
NULL


#' Correlation between number of substitutions and coalescent length units
#'
#' @description Basic correlation, to be used in the generation of approximate
#' branch lengths in coalescent units (not estimated by coalescent tree methods)
#'
#' This is an implementation of Meng Wu's idea to correct MP-EST branches.
#'
#' @param mpest_branches vector of coalescent lengths
#' @param concate_branches vector of substitutions
#'
#' @return a vector with slope and intercept
#'
branch_length_correlation <- function(mpest_branches, concate_branches){
  fit <- lm(mpest_branches ~ concate_branches)
  intercept <- coef(fit)[1]
  slope <- coef(fit)[2]
  return(c(slope,intercept))
}

#' Change values of branch lengths at tips
#'
#' @description Because estimating branch length at tips (terminal lineages) is
#'   impossible for coalescent methods, some implementations (like MP-EST)
#'   return tips of arbitrary length. These lengths don't have biological
#'   meaning but can affect the value of HRF. This function changes the values
#'   of the tips for another, equally arbitrary length. 'newval' can take a
#'   single value that is applied to all tips, a vector of equal length to the
#'   number of tips, or a phylogeny (object class 'phylo') with branch lengths
#'   in substitution rate (will be used to extrapolate coalescent branch length)
#'
#' @param branches numeric vector
#' @param is_tip logical vector
#' @param newval numeric or class 'phylo'
#'
#' @return a numeric vector
#'
fix_mpest_tips <- function(branches, is_tip, newval){
  if(is.null(newval)) return(branches)

  if(is.numeric(newval)){
    l <- length(newval)
    t <- sum(is_tip)
    if( l != 1 & l != t) {warning(paste("Replacement value for tips expected to be length 1 or", t, "\n", sep=" "))}
    branches[is_tip] <- newval
  }
  else if(is(newval, 'phylo')){
    concat_notip <- newval$edge.length[!is_tip]
    correlat <- branch_length_correlation(branches[!is_tip], concat_notip)
    branches[is_tip] <- correlat[1] * newval$edge.length[is_tip] + correlat[2]
  }
  else { warning("Class of object for tip replacement not recognized. No changes made to lengths.")}
  return(branches)
}

inverse_cf <- function(vec){
  epsilon <- 1e-4
  y <- -log((3/2)*(1-vec))
  ifelse(y < 0, epsilon, y)
}

#' Use concordance factors as approximate branch lengths
#'
#' @description Under the assumption that all gene tree discordance is due to
#' incomplete lineage sorting, this function estimates length for internal branches of a tree
#' by taking the proportion of discordant trees (i.e. 1 minus the concordance factor)
#' and solving for branch length t = -ln(3⁄2* (1-CF)).
#' It will return 0 for CF < 1/3, and
#' Infinite for CF == 1 (which isn't realistic but essentially eliminates the possibility of hemiplasy)
#'
#' @param tree a phylogeny of class 'phylo'
#' @param concordance_factors string or numeric vector of length equal to number of branches in tree
#'
#' @return numeric vector of estimated branches
#' @export
#'
convert_from_cf <- function(tree, concordance_factors){

  if(is.numeric(concordance_factors)){
    l <- length(concordance_factors)
    t <- length(tree$edge.length)
    if( l != 1 & l != t) {warning(paste("Replacement object for concordance factors expected to be length 1 or", t, "\n", sep=" "))}
    branches <- inverse_cf( concordance_factors )
  }
  else if(is.character(concordance_factors)){
    branches <- inverse_cf( with(tree, concordance_factors))
  }
  else { warning("Class of object for Concordance Factor transformation not recognized.")}

  if(sum(is.infinite(branches))> 0){warning("At least one branch length inferred from concordance factors is infinite.")}

  return(branches)
}


#' Prepare branches for HRF calculation
#' @description Take a 'phylo' object and extract the relevant branches for each
#'   HRF calculation.
#'
#'   Optionally, adjust the tip length values (if 'tip_replace' is not NULL).
#'   Because estimating branch length at tips is impossible for coalescent
#'   methods, some implementations (like MP-EST) return tips of arbitrary
#'   length. These lengths don't have biological meaning but can affect the
#'   value of HRF. This function changes the values of the tips for another,
#'   equally arbitrary length.
#'
#'   If 'tip_replace' is a single value, it's applied to all tips. If
#'   'tip_replace' is a vector, it will be passed directly to all tips (it needs
#'   to be of length equal to number of tips).
#'
#'   If 'tip_replace' is a phylogeny (an object of class 'phylo'), the
#'   replacement values will be extrapolated from the correlation of the two
#'   sets of branch lengths (the ones given in 'tree' and the ones in
#'   'tip_replace'), excluding the tips. This method assumes that the nodes in
#'   the 'tip_replace' phylogeny are in the same order as in 'tree'. While not
#'   mandatory, it is strongly recommended that the topologies of the two trees
#'   are the identical.
#'
#'   If 'concordance_factors' is not NULL, it'll modify branch lengths using
#'   concordance factors. If 'concordance_factors' is a character, it'll be taken as the
#'   name of the variable in the data frame that corresponds to the concordance factors. If it's
#'   a numeric vector, it'll be taken as the values of concordance factors (in the same order as
#'   the corresponding branches in the data frame).
#'
#'
#' @param tree a phylogeny of class 'phylo'
#' @param tip_replace NULL(default), numeric of length 1, or numeric vector of length equal to number of tips,
#' or object of class 'phylo'.
#' @param concordance_factors NULL(default), character, numeric vector of length equal to number of branches in tree
#'
#' @return a tibble
#' @export
#'
prep_branch_lengths <- function(tree, tip_replace = NULL, concordance_factors = NULL) {
  if (!is(tree, 'phylo')) {
    warning("'tree' is not of class 'phylo', prep_branch_lengths() might give bad results.")
  }

  skeleton <- dplyr::data_frame(id = 1:length(tree$edge[, 1]), from = tree$edge[, 1], to = tree$edge[, 2])

  skeleton <- dplyr::mutate(skeleton,
                            code = purrr::map2_chr(from, to, function(x, y) paste(x, y, sep = "-")),
                            descend = purrr::map(to, function(x) which(tree$edge[, 1] == x)),
                            ancest = purrr::map_int(from, function(x) match(x, tree$edge[,2])),
                            sib = purrr::map2_int(from, id, function(x, y) setdiff(which(tree$edge[, 1] == x), y)))

  is_tip <-  sapply(skeleton$descend, length) == 0
  branches <- fix_mpest_tips(tree$edge.length, is_tip, tip_replace)

  if(!is.null(concordance_factors)){
    branches <- convert_from_cf(tree, concordance_factors)
  }

  edges_df <- dplyr::transmute(skeleton,
                               code,
                               from,
                               to,
                               this_branch = branches,
                               descendants = purrr::map(descend, function(x) branches[x]),
                               ancestor = branches[ancest],
                               sibling = branches[sib])

  return(edges_df)
}

#' Calculate HRF for all branches on a tree
#'
#' @description
#' It assumes:
#' no uncertainty on species tree (A,(B,C)), one sample per species, branch lengths in coalescent units.
#'
#' If 'mutation' is a numeric of length 1, forward (0->1) and reverse (1->0) mutation rates are the same. If
#' it's of length 2, the first element is the forward- and the second is the reverse mutation rate.
#'
#' 'mode' can be one of three:
#' 'minimal' -- Reversals are not allowed. Ancestral (ABC) time is assumed to be infinite.
#' 'standard' -- Reversals allowed. Ancestral (ABC) time is assumed to be at least 4. This is equivalent to assuming that
#' coalescence in the ABC ancestor is very likely (most coalescence happens before this time).
#' 'strict' -- Reversals allowed. No changes to ABC time. If this time is very short, probabilities of reversals
#' and of hemiplasy are affected. The ABC branch (where the first 0->1 would happen in a reversal scenario) and the
#' AC branch (where the only 0->1 happens in hemiplasy) are very short, lowering the probability of events on them.
#'
#'
#'
#' @param edges a tibble (or data_frame)
#' @param mutation numeric (forward, reverse population-wide mutation rate)
#' @param mode character
#'
#' @return a tibble
#' @export
#' @seealso \code{\link{pr_homoplasy}}, \code{\link{pr_hemiplasy}}, \code{\link{calc_hrf}}
#'
tree_hrf <- function(edges, mutation = 0.01, mode = "standard") {
  model_avail <- c("standard", "minimal", "strict")
  if(!mode %in% model_avail) warning(paste("Requested model", mode,"in tree_hrf() not recognized."))

  if(is.null(edges$this_branch) | is.null(edges$descendants) | is.null(edges$ancestor) | is.null(edges$sibling)){
    warning("At least one variable is missing in the 'edges' data frame. You might need to run prep_branch_lengths() first.")
  }
 return(
      dplyr::mutate(edges,
                    hrf = purrr::pmap_dbl(list(this_branch, descendants, ancestor, sibling), calc_hrf, mutation_rate = mutation, model = mode))
    )
}

#' Convert to 'treedata'
#' @description Take a data frame and 'phylo' object and turn them into
#' a 'treedata' object (from 'treeio' package)
#' @param tree an object of class "phylo"
#' @param df accesory data to be added
#'
#' @return an object of class treedata
#' @export
#'
to_treedata <- function(tree, df) {
  if(!requireNamespace('treeio', quietly = TRUE)) {
    message("The to_treedata() function needs the 'treeio' packaged. Please install it.")
    return(NA)
  }
  ph <- ape::as.phylo(tree)
  dfout <- dplyr::mutate(df, node = to)
  new("treedata", phylo = ph, data = dfout)
}

#' Convert to 'phylo4' class
#'
#' @inheritParams to_treedata
#'
#' @return an object of class phylo4
#' @export
#'
to_phylo4d <- function(tree, df) {

  earliest <- dplyr::filter(df, is.na(ancestor))
  edgeno <- dplyr::slice(earliest, 1)
  anc_node <- dplyr::data_frame(code = paste("0", edgeno$from, sep = "-"), to = edgeno$from, this_branch = Inf, descendants = list(earliest$this_branch), hrf = 0)
  new_df <- dplyr::bind_rows(df, anc_node)
  legacy_df <- as.data.frame(new_df)
  rownames(legacy_df) <- new_df$to

  phylobase::phylo4d(tree, all.data = legacy_df)
}
