# Impute the combined parental ALT dosage from child genotypes and sibling IBD.
#
# For a singleton, E[G_parent | g] = g + 2p. For a sibling pair, the
# conditional expectation depends on the exact IBD0/1/2 state. The formulas
# match impute_parental_genotypes_from_siblings.qmd.

.impute_parental_alt_dosage <- function(g1, g2, ibd, p_alt) {
  if (length(g1) != length(g2) ||
      length(g1) != length(ibd)) {
    stop("g1, g2, and ibd must have equal lengths.")
  }
  if (anyNA(g1) || anyNA(g2) || anyNA(ibd) ||
      any(!(g1 %in% 0:2)) ||
      any(!(g2 %in% 0:2)) ||
      any(!(ibd %in% 0:2)) ||
      length(p_alt) != 1L ||
      !is.finite(p_alt) ||
      p_alt < 0 ||
      p_alt > 1) {
    stop("Invalid sibling genotype, IBD, or allele-frequency input.")
  }

  incompatible <- (
    ibd == 1L & abs(g1 - g2) == 2L
  ) | (
    ibd == 2L & g1 != g2
  )
  if (any(incompatible)) {
    stop(
      "Sibling genotypes are incompatible with the supplied IBD state."
    )
  }

  parent <- rep(NA_real_, length(g1))

  use_ibd0 <- ibd == 0L
  parent[use_ibd0] <- g1[use_ibd0] + g2[use_ibd0]

  use_ibd2 <- ibd == 2L
  parent[use_ibd2] <- g1[use_ibd2] + 2 * p_alt

  use_double_het <- ibd == 1L & g1 == 1L & g2 == 1L
  parent[use_double_het] <- 1 + 2 * p_alt

  use_resolved_ibd1 <- ibd == 1L & !use_double_het
  shared_allele <- ifelse(
    g1[use_resolved_ibd1] == 0L |
      g2[use_resolved_ibd1] == 0L,
    0,
    1
  )
  parent[use_resolved_ibd1] <-
    g1[use_resolved_ibd1] +
    g2[use_resolved_ibd1] -
    shared_allele +
    p_alt

  parent
}

#' Impute combined parental ALT dosages.
#'
#' @param genotype Numeric n by p matrix of child ALT dosages in 0, 1, or 2.
#' @param ibd_state Numeric n_pair by p matrix of exact IBD states. Row names
#'   must identify the paired families.
#' @param pedigree Data frame with one row per genotype row and columns
#'   family_id and individual_id. Families may contain one or two children.
#'
#' @return A list containing parental_alt_dosage, p_alt, pair_family_id,
#'   singleton_rows, first_sibling_rows, and second_sibling_rows.
impute_parental_genotypes <- function(
    genotype,
    ibd_state,
    pedigree) {
  if (!is.matrix(genotype) ||
      !is.numeric(genotype) ||
      nrow(genotype) == 0L ||
      ncol(genotype) == 0L ||
      anyNA(genotype) ||
      any(!(genotype %in% 0:2))) {
    stop("genotype must be a complete numeric n by p dosage matrix.")
  }
  if (!is.data.frame(pedigree) ||
      nrow(pedigree) != nrow(genotype) ||
      !all(c("family_id", "individual_id") %in% names(pedigree))) {
    stop(
      "pedigree must align with genotype and contain family_id ",
      "and individual_id."
    )
  }

  family_id <- as.character(pedigree$family_id)
  individual_id <- as.character(pedigree$individual_id)
  if (anyNA(family_id) ||
      anyNA(individual_id) ||
      any(family_id == "") ||
      any(individual_id == "") ||
      anyDuplicated(individual_id)) {
    stop("Pedigree identifiers must be complete and individual IDs unique.")
  }

  family_levels <- unique(family_id)
  family_index <- match(family_id, family_levels)
  family_size <- tabulate(
    family_index,
    nbins = length(family_levels)
  )
  if (any(!(family_size %in% 1:2))) {
    stop("Every analysis family must contain one child or one sibling pair.")
  }

  family_rows <- split(
    seq_len(nrow(genotype)),
    factor(family_index, levels = seq_along(family_levels))
  )
  singleton_rows <- unlist(
    family_rows[family_size == 1L],
    use.names = FALSE
  )
  pair_rows <- family_rows[family_size == 2L]
  pair_family_id <- family_levels[family_size == 2L]
  n_pairs <- length(pair_rows)
  n_snps <- ncol(genotype)

  if (n_pairs == 0L) {
    if (!is.null(ibd_state) &&
        (!is.matrix(ibd_state) ||
          nrow(ibd_state) != 0L ||
          ncol(ibd_state) != n_snps)) {
      stop("With no pairs, ibd_state must be NULL or a 0 by p matrix.")
    }
  } else {
    if (!is.matrix(ibd_state) ||
        !is.numeric(ibd_state) ||
        nrow(ibd_state) != n_pairs ||
        ncol(ibd_state) != n_snps ||
        anyNA(ibd_state) ||
        any(!(ibd_state %in% 0:2))) {
      stop("ibd_state must be a complete numeric n_pair by p matrix.")
    }

    ibd_family_id <- rownames(ibd_state)
    if (is.null(ibd_family_id) ||
        anyDuplicated(ibd_family_id) ||
        anyNA(match(pair_family_id, ibd_family_id))) {
      stop("ibd_state row names must uniquely identify all paired families.")
    }
    ibd_state <- ibd_state[
      match(pair_family_id, ibd_family_id),
      ,
      drop = FALSE
    ]

    genotype_snp_id <- colnames(genotype)
    ibd_snp_id <- colnames(ibd_state)
    if (!is.null(genotype_snp_id) && !is.null(ibd_snp_id)) {
      if (anyDuplicated(genotype_snp_id) ||
          anyDuplicated(ibd_snp_id) ||
          anyNA(match(genotype_snp_id, ibd_snp_id))) {
        stop("Named genotype and IBD SNP columns must match uniquely.")
      }
      ibd_state <- ibd_state[
        ,
        match(genotype_snp_id, ibd_snp_id),
        drop = FALSE
      ]
    }
  }

  p_alt <- colMeans(genotype) / 2
  parent <- matrix(
    NA_real_,
    nrow = nrow(genotype),
    ncol = n_snps,
    dimnames = dimnames(genotype)
  )

  if (length(singleton_rows) > 0L) {
    parent[singleton_rows, ] <- sweep(
      genotype[singleton_rows, , drop = FALSE],
      MARGIN = 2L,
      STATS = 2 * p_alt,
      FUN = "+"
    )
  }

  first_sibling_rows <- integer()
  second_sibling_rows <- integer()
  if (n_pairs > 0L) {
    first_sibling_rows <- vapply(
      pair_rows,
      function(rows) rows[[1L]],
      integer(1L)
    )
    second_sibling_rows <- vapply(
      pair_rows,
      function(rows) rows[[2L]],
      integer(1L)
    )
    pair_parent <- matrix(
      NA_real_,
      nrow = n_pairs,
      ncol = n_snps,
      dimnames = list(pair_family_id, colnames(genotype))
    )

    for (snp in seq_len(n_snps)) {
      pair_parent[, snp] <- .impute_parental_alt_dosage(
        g1 = genotype[first_sibling_rows, snp],
        g2 = genotype[second_sibling_rows, snp],
        ibd = ibd_state[, snp],
        p_alt = p_alt[[snp]]
      )
    }
    parent[first_sibling_rows, ] <- pair_parent
    parent[second_sibling_rows, ] <- pair_parent
  }

  list(
    parental_alt_dosage = parent,
    p_alt = p_alt,
    pair_family_id = pair_family_id,
    singleton_rows = singleton_rows,
    first_sibling_rows = first_sibling_rows,
    second_sibling_rows = second_sibling_rows
  )
}
