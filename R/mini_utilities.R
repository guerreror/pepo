critical_mutation <- function(t1, t2, x = 1) {
  3 * exp(-t2) / ((5 + 3 * t1) * (4 + 3 * t1 + 3 * t2) * x)
}

mini_pepo <- function(mu, t1, t2, rev = F) {
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
