% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrf_on_trees.R
\name{prep_branch_lengths}
\alias{prep_branch_lengths}
\title{Prepare branches for HRF calculation}
\usage{
prep_branch_lengths(tree, tip_replace = NULL,
  concordance_factors = NULL)
}
\arguments{
\item{tree}{a phylogeny of class 'phylo'}

\item{tip_replace}{NULL(default), numeric of length 1, or numeric vector of length equal to number of tips,
or object of class 'phylo'.}

\item{concordance_factors}{NULL(default), character, numeric vector of length equal to number of branches in tree}
}
\value{
a tibble
}
\description{
Take a 'phylo' object and extract the relevant branches for each
  HRF calculation.

  Optionally, adjust the tip length values (if 'tip_replace' is not NULL).
  Because estimating branch length at tips is impossible for coalescent
  methods, some implementations (like MP-EST) return tips of arbitrary
  length. These lengths don't have biological meaning but can affect the
  value of HRF. This function changes the values of the tips for another,
  equally arbitrary length.

  If 'tip_replace' is a single value, it's applied to all tips. If
  'tip_replace' is a vector, it will be passed directly to all tips (it needs
  to be of length equal to number of tips).

  If 'tip_replace' is a phylogeny (an object of class 'phylo'), the
  replacement values will be extrapolated from the correlation of the two
  sets of branch lengths (the ones given in 'tree' and the ones in
  'tip_replace'), excluding the tips. This method assumes that the nodes in
  the 'tip_replace' phylogeny are in the same order as in 'tree'. While not
  mandatory, it is strongly recommended that the topologies of the two trees
  are the identical.

  If 'concordance_factors' is not NULL, it'll modify branch lengths using
  concordance factors. If 'concordance_factors' is a character, it'll be taken as the
  name of the variable in the data frame that corresponds to the concordance factors. If it's
  a numeric vector, it'll be taken as the values of concordance factors (in the same order as
  the corresponding branches in the data frame).
}
