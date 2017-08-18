#' @include mutation_functions.R
NULL

#' Probability of incongruence due to hemiplasy in a three-species clade.
#'
#' @description
#' It assumes:
#' no uncertainty on species tree ((A,B),C), one sample per species, observed states of
#' binary trait were a=1, b=0, and c=1, branch lengths in coalescent units.
#'
#' Roughly, it computes the product of the probabilities of events necessary for the observed
#' incongruence. Specifically, the probability of samples a and c coalescing (which includes incomplete
#' lineage sorting in the internal AB branch), and one mutation on the ac branch (which implies no mutations
#' on other branches of the tree).
#'
#' See documentation of tree_hrf() for description of the three models implemented.
#'
#' @param this Numeric. Branch length of focal internal branch.
#' @param desc_a,desc_b Lengths of branches descending from focal branch.
#' @param sib Length of sibling branch (leading to tip node C).
#' @param anc Length of ancestral branch (the ABC ancestor)
#' @param mu_rate Population-wide mutation rate
#' @param mode String. One of the three models implemented: "standard" (default), "minimal", or "strict".
#'
#' @return Numeric, a probability
#' @seealso \code{\link{pr_homoplasy}}, \code{\link{single_hrf}}, \code{\link{tree_hrf}}
#'
pr_hemiplasy <- function(this, desc_a, desc_b, sib, anc, mu_rate, mode= "standard") {

    if (mode == "minimal") {
      nomut_anc <- 1
    } else {
      nomut_anc <- 1 - pr_mutation(mu = mu_rate, t = branch_k3(t0 = 0, tm = anc, branch = "anc"))
    }

    pr_no_coal <- exp(-this)
    mut_ac <- pr_mutation(mu = mu_rate, t = branch_k3(t0 = 0, tm = anc, branch = "inter"))
    nomut_a <- 1 - pr_mutation(mu = mu_rate, t = branch_k3(t0 = desc_a + this, tm = anc, branch = "short"))
    nomut_c <- 1 - pr_mutation(mu = mu_rate, t = branch_k3(t0 = sib, tm = anc, branch = "short"))
    nomut_b <- 1 - pr_mutation(mu = mu_rate, t = branch_k3(t0 = desc_b + this, tm = anc, branch = "long"))

    (1/3) * pr_no_coal * mut_ac * nomut_a * nomut_c * nomut_b * nomut_anc
}

#' Probability of incongruence due to homoplasy in a three-species clade.
#'
#' @description
#' It assumes:
#' no uncertainty on species tree ((A,B),C), one sample per species, observed states of
#' binary trait were a=1, b=0, and c=1, branch lengths in coalescent units.
#'
#' Roughly, it computes the product of the probabilities of events necessary for the observed
#' incongruence. Specifically, mutations along the A and C branches (implying no mutations
#' in other branches) or --depending on \code{mode}- mutations along the ABC ancestor and B
#' (again, no other mutations).
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

    hit_ac <- pr_coal * pr_double_hit(mu_coal, "ac") + (1 - pr_coal) * pr_double_hit(mu_noco, "ac")

    if (mode != "minimal") {
        hit_banc <- pr_coal * pr_double_hit(mu_coal, "banc") + (1 - pr_coal) * pr_double_hit(mu_noco, "banc")
        return(hit_ac + hit_banc)
    } else {
        return(hit_ac)
    }
}

#' Ratio of hemiplasy to homoplasy in a three-species clade.
#'
#' @description
#' It assumes:
#' no uncertainty on species tree ((A,B),C), one sample per species, branch lengths in coalescent units.
#'
#' Calls \code{pr_homoplasy()} and \code{pr_hemiplasy()} for the given values. If \code{both_sides == TRUE}
#' (default), the total probability of each scenario is the sum of probabilities for two incongruencies:
#' when A and C share the derived state, and when B and C share the derived state.
#' Otherwise, it only calculates the first case (assuming the first value in the
#' \code{descendants} vector as A).
#'
#' See documentation for tree_hrf() for description of the three models implemented.
#'
#' @param this_branch Numeric. Branch length of focal internal branch.
#' @param descendants Numeric vector of size 2. Lengths of branches descending from focal branch.
#' @param sibling Numeric. Length of sibling branch (leading to tip node C).
#' @param ancestor Numeric. Length of ancestral branch (the ABC ancestor)
#' @param mutation_rate Numeric. Population-wide mutation rate
#' @param model String. One of the three models implemented: "standard" (default), "minimal", or "strict".
#'
#' @return Numeric, a ratio
#' @seealso \code{\link{pr_homoplasy}}, \code{\link{pr_hemiplasy}}, \code{\link{tree_hrf}}
single_hrf <- function(this_branch, descendants, ancestor, sibling, mutation_rate, model = "standard") {
    if (length(descendants) != 2) return(NA) # this branch doesn't have two descendants, HRF can't be computed
    if (is.na(this_branch)) return(NA) # this branch doesn't have a positive length. Happens at root. No HRF calculated.
    if (is.na(ancestor)) ancestor <- Inf # if ancestral branch is undefined, it is assumed to be infinite (ancestor may be root).

    if (model == "minimal") ancestor <- Inf
    if (model == "standard") ancestor <- max(ancestor, 4)  # 97% of all genealogies of k=3 have coalesced by time 4

    a <- descendants[1]
    b <- descendants[2]

    pe_a <- pr_hemiplasy(this = this_branch, desc_a = a, desc_b = b, sib = sibling, anc = ancestor, mu_rate = mutation_rate, mode = model)
    pe_b <- pr_hemiplasy(this = this_branch, desc_a = b, desc_b = a, sib = sibling, anc = ancestor, mu_rate = mutation_rate, mode = model)
    pe <- pe_a + pe_b


    po_a <- pr_homoplasy(this = this_branch, desc_a = a, desc_b = b, sib = sibling, anc = ancestor, mu_rate = mutation_rate, mode = model)
    po_b <- pr_homoplasy(this = this_branch, desc_a = b, desc_b = a, sib = sibling, anc = ancestor, mu_rate = mutation_rate, mode = model)
    po <- po_a + po_b

    return(pe/po)
}
