#' Calculate the probability of at least one mutation on a branch
#'
#' @param mu mutation rate in coalescent units, 2N*per_generation_mutation_rate
#' @param t branch length, in coalescent units
#' @return Numeric, the probability of at least one mutation.
#' @examples
#' pr_mutation(0.001, 1)

pr_mutation <- function(mu, t) {
  (1 - exp(-mu * t))
}

#' Calculate the length of a branch that ends with coalescence of two samples (i.e., coal rate = 1)
#'
#' @param t0 Numeric. Time at start of coalescent process (i.e., constant length to be added to terminal branches; in coal units)
#' @param tm Numberic. Maximum time to coalescence (ABC ancestor branch length)
#' @return Numeric, a branch length in coalescent units

branch_k2 <- function(t0 = 0, tm = Inf) {
  coalcdf <- 1 - exp(-tm)
  if (tm == Inf) {
    return(t0 + 1)
  }
  (1 + t0 - exp(-tm) * (1 + t0 + tm)) / coalcdf
}

#' Calculate branch lengths in genealogy that involves ILS
#'
#' @description
#' Estimates average length of any branch in the genealogy of the sample (a,b,c), given that
#' the species tree is ((A,B),C) and there was no coalescence of the (a,b) sample in the AB ancestral species.
#' In other words, coalescence of the three samples (a,b,c) happens in the ABC ancestor.
#'
#' @inheritParams branch_k2
#' @param branch String. Branch to be returned: first ('short') or second ('long') coalescent events, mean of those two ('avg'), internal ('inter'), or ancestral branch ('anc').
#' @return Numeric, the probability of at least one mutation.
#' @examples
#' branch_k3(1, 10, 'short')

branch_k3 <- function(t0 = 0, tm = Inf, branch = "avg") {

  coalcdf <- (1/2) * (2 + exp(-3 * tm) - 3 * exp(-tm))
  shortb <- (t0 + (1/12) * exp(-3 * tm) * (5 - 9 * exp(2 * tm) + 4 * exp(3 * tm) + 6 * tm)) / coalcdf
  longb <- (t0 + (1/6) * exp(-3 * tm) * (1 + 8 * exp(3 * tm) + 3 * tm - 9 * exp(2 * tm) * (1 + tm))) / coalcdf

  if (tm == Inf) {
    shortb <- t0 + 1/3
    longb <-  t0 + 4/3
  }

  out <- switch(branch,
                short = shortb,
                long = longb,
                inter = longb - shortb,
                anc = tm - longb,
                avg = (2/3) * shortb + (1/3) * longb,
                NA)

  if (is.na(out)) print("In 'branch_k3(): Value of 'branch' parameter not recognized.")

  return(out)
}

#' Calculate probability of two mutations on three-species tree
#' @description
#' One of two types of double-mutation events will be returned, depending on the value of \code{where}.
#' Mutation and back-mutation are assumed to have the same rate.
#'
#' "ac"--> product of the following probabilities: mutation on the [a] and [c] branches, no mutations on [b], ancestral, or internal branches.
#' "banc" --> product of the following probabilities: mutation on the [abc] ancestor and on the [b] branch, no mutation on [a],[c] or internal branches.
#'
#' @param mu_list List with "anc", "a", "b", "c", and "ab" entries, each of which is a mutation probability.
#' @param where String. "ac" or "banc".
#'
#' @return Numeric. Probability
#' @seealso homoplasy_mutations()
pr_double_hit <- function(mu_list, where) {
  out <- NA
  if (where == "ac") out <-  mu_list$a * mu_list$c * (1 - mu_list$anc) * (1 - mu_list$b) * (1 - mu_list$ab)
  if (where == "banc") out <- mu_list$b * mu_list$anc * (1 - mu_list$a) * (1 - mu_list$c) * (1 - mu_list$ab)
  if (is.na(out)) print("In pr_double_hit(): Parameter value for 'where' not recognized")
  return(out)
}

#' Get probabilities of mutation along all branches
#' @description
#' Given a species tree, calculate the probability of mutation on each of the branches of the gene tree.
#'
#' @inheritParams pr_hemiplasy
#' @param coal logical. Did coalescence in the AB ancestor happen?
#'
#' @return List of size 5
homoplasy_mutations <- function(this, desc_a, desc_b, sib, anc, mu_rate, coal = T, mode) {
  reversals <- TRUE
  if (mode == "minimal")
    reversals <- FALSE

  if (coal) {
    # Return probabilities for genealogies assuming coalescence of AB in this branch
    # t for the ab branch is the sum of ancestral ab branch up to the start of the ancestral population
    # and the branch in the ancestor
    out <- list(anc = ifelse(reversals, pr_mutation(mu = mu_rate, t = anc - branch_k2(t0 = 0, tm = anc)), 0),
                a = pr_mutation(mu = mu_rate, t = branch_k2(t0 = desc_a, tm = this)),
                b = pr_mutation(mu = mu_rate, t = branch_k2(t0 = desc_b, tm = this)),
                c = pr_mutation(mu = mu_rate, t = branch_k2(t0 = sib, tm = anc)),
                ab = pr_mutation(mu = mu_rate, t = branch_k2(t0 = 0, tm = anc) + (this - branch_k2(t0 = 0, tm = this))))
  } else {
    # Return probabilities for genealogies assuming NO coalescence of AB
    # A,B coalesce in the ancestral branch with C, so here we average probabilities of the three possible gene trees
    out <- list(anc = ifelse(reversals, pr_mutation(mu = mu_rate, t = branch_k3(t0 = 0, tm = anc, branch = "anc")), 0),
                a = pr_mutation(mu = mu_rate, t = branch_k3(t0 = desc_a + this, tm = anc, branch = "avg")),
                b = pr_mutation(mu = mu_rate, t = branch_k3(t0 = desc_b + this, tm = anc, branch = "avg")),
                c = pr_mutation(mu = mu_rate, t = branch_k3(t0 = sib, tm = anc, branch = "avg")),
                ab = pr_mutation(mu = mu_rate, t = branch_k3(t0 = 0, tm = anc, branch = "inter")))
  }
  return(out)
}

