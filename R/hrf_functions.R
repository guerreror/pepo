#' @include mutation_functions.R
NULL

#' Probability of incongruence due to hemiplasy in a three-species clade.
#'
#' @description
#' It assumes:
#' no uncertainty on species tree (A,(B,C)), one sample per species, observed states of
#' binary trait were a=1, b=0, and c=1, branch lengths in coalescent units.
#'
#' Roughly, it computes the product of the probabilities of events necessary for the observed
#' incongruence. Specifically, the probability of samples a and c coalescing (which includes incomplete
#' lineage sorting in the internal BC branch), and one mutation on the ac branch (which implies no mutations
#' on other branches of the tree).
#'
#' See documentation of tree_hrf() for description of the three models implemented.
#'
#' @param this Numeric. Branch length of focal internal branch.
#' @param desc_a,desc_b Lengths of branches descending from focal branch.
#' @param sib Length of sibling branch (leading to tip node C).
#' @param anc Length of ancestral branch (the ABC ancestor)
#' @param mu_rate Numeric of length 2. Population-wide mutation rate (forward, reverse)
#' @param mode String. One of the three models implemented: "standard" (default), "minimal", or "strict".
#'
#' @return Numeric, a probability
#' @seealso \code{\link{pr_homoplasy}}, \code{\link{single_hrf}}, \code{\link{tree_hrf}}
#'
pr_hemiplasy <- function(this, desc_a, desc_b, sib, anc, mu_rate, mode= "standard") {

   treeb <- list(
      anc = ifelse(mode == "minimal", 0, nu_k3_anc(t3 = anc, mu = mu_rate[1])),
      a = nu_k3_short(t0 = desc_a + this, t3 = anc, mu = mu_rate),
      b = nu_k3_long(t0 = desc_b + this, t3 = anc, mu = mu_rate),
      c = nu_k3_short(t0 = sib, t3 = anc, mu = mu_rate),
      ac = nu_k3_inter(t3 = anc, mu = mu_rate))

    pr_no_coal <- exp(-this)

    (1/3) * pr_no_coal * treeb$ac[1] * (1-treeb$a[2]) * (1-treeb$c[2]) * (1-treeb$b[1]) * (1-treeb$anc[1])
}

#' Probability of incongruence due to homoplasy in a three-species clade.
#'
#' @description
#' It assumes:
#' no uncertainty on species tree (A,(B,C)), one sample per species, observed states of
#' binary trait were a=1, b=0, and c=1, branch lengths in coalescent units.
#'
#' Roughly, it computes the product of the probabilities of events necessary for the observed
#' incongruence. Specifically, mutations along the A and C branches (implying no mutations
#' in other branches) or --depending on \code{mode}- mutations along the ABC ancestor and B.
#'
#' See documentation of tree_hrf() for description of the three models implemented.
#'
#' @inheritParams pr_hemiplasy
#'
#' @return Numeric, a probability
#' @seealso \code{\link{pr_hemiplasy}}, \code{\link{single_hrf}}, \code{\link{tree_hrf}}
#'
pr_homoplasy <- function(this, desc_a, desc_b, sib, anc, mu_rate, mode = "standard") {
    pr_coal <- 1 - exp(-this)

    mu_noco <- homoplasy_mutations(this, desc_a, desc_b, sib, anc, mu_rate, coal = F, mode)
    mu_coal <- homoplasy_mutations(this, desc_a, desc_b, sib, anc, mu_rate, coal = T, mode)

    hit_ac <- pr_coal * pr_double_hit(mu_coal, "ac") +
      (1 - pr_coal) * ((1/3)*pr_double_hit(mu_noco$alpha, "ac") +
                       (1/3)*pr_double_hit(mu_noco$beta, "ac") +
                       (1/3)*pr_double_hit(mu_noco$gamma, "ac"))

    if (mode != "minimal") {
        hit_banc <- pr_coal * pr_double_hit(mu_coal, "banc") +
          (1 - pr_coal) * ((1/3)*pr_double_hit(mu_noco$alpha, "banc") +
                           (1/3)*pr_double_hit(mu_noco$beta, "banc") +
                           (1/3)*pr_double_hit(mu_noco$gamma, "banc"))
        return(hit_ac + hit_banc)
    } else {
        return(hit_ac)
    }
}

#' Calculate the Hemiplasy Risk Factor in a three-species clade.
#'
#' @description
#' HRF is the expected fraction of incongruency due to hemiplasy (vs. homoplasy). It assumes:
#' no uncertainty on species tree (A,(B,C)), one sample per species, branch lengths in coalescent units.
#'
#' Calls \code{pr_homoplasy()} and \code{pr_hemiplasy()} for the given values. The total probability of
#' each scenario is the sum of probabilities for two incongruencies:
#' when A and C share the derived state, and when B and C share the derived state.
#'
#' If 'mutation_rate' is a numeric of length 1, forward (0->1) and reverse (1->0) mutation rates are the same.
#' If 'mutation_rate' is of length 2, the first element is the forward- and the second is the reverse mutation rate.
#'
#' See documentation for tree_hrf() for description of the three models implemented.
#'
#' @param this_branch Numeric. Branch length of focal internal branch.
#' @param descendants Numeric vector of size 2. Lengths of branches descending from focal branch.
#' @param sibling Numeric. Length of sibling branch (leading to tip node A).
#' @param ancestor Numeric. Length of ancestral branch (the ABC ancestor)
#' @param mutation_rate Numeric (forward, reverse population-wide mutation rate)
#' @param model String. One of the three models implemented: "standard" (default), "minimal", or "strict".
#'
#'
#' @return Numeric, a ratio
#' @seealso \code{\link{pr_homoplasy}}, \code{\link{pr_hemiplasy}}, \code{\link{tree_hrf}}
calc_hrf <- function(this_branch, descendants, ancestor, sibling, mutation_rate, model = "standard", pepo=F) {
    if (length(descendants) != 2) return(NA) # this branch doesn't have two descendants, HRF can't be computed
    if (is.na(this_branch)) return(NA) # this branch doesn't have a positive length. Happens at root. No HRF calculated.
    if (is.na(ancestor)) ancestor <- Inf # if ancestral branch is undefined, it is assumed to be infinite (ancestor may be root).

    if (model == "minimal") ancestor <- Inf
    if (model == "standard") ancestor <- max(ancestor, 4)  # 97% of all genealogies of k=3 have coalesced by time 4

    a <- descendants[1]
    b <- descendants[2]

    mu <- fix_mu_vals(mutation_rate)

    pe_a <- pr_hemiplasy(this = this_branch, desc_a = a, desc_b = b, sib = sibling, anc = ancestor, mu_rate = mu, mode = model)
    pe_b <- pr_hemiplasy(this = this_branch, desc_a = b, desc_b = a, sib = sibling, anc = ancestor, mu_rate = mu, mode = model)
    pe <- pe_a + pe_b


    po_a <- pr_homoplasy(this = this_branch, desc_a = a, desc_b = b, sib = sibling, anc = ancestor, mu_rate = mu, mode = model)
    po_b <- pr_homoplasy(this = this_branch, desc_a = b, desc_b = a, sib = sibling, anc = ancestor, mu_rate = mu, mode = model)
    po <- po_a + po_b
    if(pepo)    return(pe_a/po_a)

    return(pe/(pe+po))
}
