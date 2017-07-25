#' Calculates the probability of at least one mutation on a branch
#'
#' @param mu mutation rate in coalescent units, 2N*per_generation_mutation_rate
#' @param t branch length, in coalescent units
#' @return Numeric, the probability of at least one mutation.
#' @examples
#' pr_mutation(0.001, 1)

pr_mutation <- function(mu, t) {
    (1 - exp(-mu * t))
}

#' Calculates the length of a branch that ends with coalescence of two samples (i.e., coal rate = 1)
#'
#' @param t0 time to beginning of ancestral lineage (i.e., a constant branch length to add) in coalescent units
#' @param tm maximum time to coalescence (total time of ancestral lineage) in coalescent units
#' @return Numeric, a branch length in coalescent units

branch_k2 <- function(t0 = 0, tm = Inf) {
    coalcdf <- 1 - exp(-tm)
    if (tm == Inf) {
        return(t0 + 1)
    }
    (1 + t0 - exp(-tm) * (1 + t0 + tm)) / coalcdf
}

#' Calculates the length of a branch that ends in a population where three samples coalesce
#'
#'
#' @param t0 time to ancestor, or time to speciation event (i.e., a constant branch length, in coalescent units)
#' @param tm maximum time to coalescence (ancestral branch length)
#' @param whichBranch specifies if the time to the first ('short'), second ('long'), internal ('inter'), ancestral ('anc'), mean ('avg') coalescent event.
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

pr_double_hit <- function(mu_vec, where) {
  out <- NA
  if (where == "ac") out <-  mu_vec$a * mu_vec$c * (1 - mu_vec$anc) * (1 - mu_vec$b) * (1 - mu_vec$ab)
  if (where == "banc") out <- mu_vec$b * mu_vec$anc * (1 - mu_vec$a) * (1 - mu_vec$c) * (1 - mu_vec$ab)
  if (is.na(out)) print("In pr_double_hit(): Parameter value for 'where' not recognized")
  return(out)
}

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

pr_hemiplasy <- function(this, desc_a, desc_b, sib, anc, mu_rate, mode) {

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

tree_hrf <- function(edges, mutation = 0.01, mode = "standard") {
    edges %>% mutate(hrf = pmap_dbl(list(.$this_branch, .$descendants, .$ancestor, .$sibling), single_hrf, mutation_rate = mutation, model = mode))
}

prep_branch_lengths <- function(tree) {

    edges <- data_frame(id = 1:length(tree$edge[, 1]),
                        from = tree$edge[, 1],
                        to = tree$edge[, 2],
                        code = map2_chr(from, to, function(x, y) paste(x, y, sep = "-")),
                        this_branch = tree$edge.length)

    edges <- edges %>% mutate(descend = map(to, function(x) which(tree$edge[, 1] == x)), ancest = map_int(from, function(x) match(x, tree$edge[,
        2])), sib = map2_int(from, id, function(x, y) setdiff(which(tree$edge[, 1] == x), y)), des = map(descend, function(x) tree$edge.length[x]),
        anc = tree$edge.length[ancest], sib = tree$edge.length[sib]) %>% transmute(code, from, to, this, des, anc, sib)

    return(edges)
}

to_phylo4d <- function(tree, df) {
  if (!requireNamespace("phylo4d", quietly = TRUE)) {
    stop("phylo4d needed for this function to work. Please install it.",
         call. = FALSE)
  }

    earliest <- df %>% filter(is.na(anc))
    edgeno <- earliest %>% slice(1) %>% .$from
    anc_node <- data_frame(code = paste("0", edgeno, sep = "-"), to = edgeno, this = Inf, des = list(earliest$this), hrf = 0)
    new_df <- bind_rows(anc_node, df)
    legacy_df <- as.data.frame(new_df)
    rownames(legacy_df) <- new_df$to

    phylo4d(tree, all.data = legacy_df)
}

to_treedata <- function(tree, df) {
    require(ggtree)
    ph <- as.phylo(tree)
    dfout <- df %>% mutate(node = to)
    new("treedata", phylo = ph, data = dfout)
}

mini_pepo <- function(mu, t1, t2) {
    3 * exp(-t2) / ((5 + 3 * t1) * (4 + 3 * t1 + 3 * t2) * mu)
}

critical_mutation <- function(t1, t2, x = 1) {
    3 * exp(-t2) / ((5 + 3 * t1) * (4 + 3 * t1 + 3 * t2) * x)
}

mini_hrf <- function(mu, t1, t2, rev = F) {
    t3 <- Inf
    model <- "minimal"
    if (rev) {
      t3 <- 4
      model <- "standard"
    }

    pe <- pr_hemiplasy(t2, t1, t1, t1 + t2, t3, mu, mode = model)
    po <- pr_homoplasy(t2, t1, t1, t1 + t2, t3, mu, mode = model)

    pe / po
}

coarse_log_likelihood <- function(hrf, level = 2) {
    if (is.na(hrf)) {
        return(NA)
    }
    if (hrf > level) {
        return("hemiplasy")
    }
    if (hrf < (1 / level)) {
        return("homoplasy")
    }
    return("ns")
}
