# Regenerate the fixed chr18-22 family panel after download_reference.sh.
# Run from this folder: Rscript simulate_genotypes.R 18

chromosome_seeds <- c(
  "18" = 20261744L,
  "19" = 20261745L,
  "20" = 20261748L,
  "21" = 20261746L,
  "22" = 20261747L
)

simulate_chromosome <- function(chromosome, n_families = 5000L,
                               output_dir = "data/simulation") {
  .libPaths(unique(c(Sys.getenv("SIM1000G_R_LIB", unset = ".Rlib-sim1000G"), .libPaths())))
  stopifnot(requireNamespace("sim1000G", quietly = TRUE))
  data.table::setDTthreads(1L)
  reference_dir <- "data/simulation/reference"
  selection <- read.delim("data/simulation/snp_selection.tsv")
  selection <- selection[selection$chromosome == chromosome, ]
  seed <- unname(chromosome_seeds[as.character(chromosome)])
  stopifnot(nrow(selection) > 0L, !is.na(seed))
  output <- file.path(output_dir, sprintf("families.chr%d.rds", chromosome))
  if (file.exists(output)) stop("Output already exists: ", output)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # sim1000G expects chromosome, position, rate, and centimorgan columns.
  map <- read.table(file.path(reference_dir, sprintf("chr%d.map", chromosome)))
  map <- map[map$V1 == chromosome, ]
  map <- map[order(map$V4), ]
  converted <- data.frame(chromosome = map$V1, position = map$V4,
                          rate = NA_real_, centimorgan = map$V3)
  map_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(map_file), add = TRUE)
  write.table(converted, map_file, sep = "\t", quote = FALSE, row.names = FALSE)
  sim1000G::readGeneticMapFromFile(map_file)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  vcf <- sim1000G::readVCF(file.path(reference_dir, sprintf("chr%d.vcf.gz", chromosome)),
                          maxNumberOfVariants = nrow(selection), min_maf = 0.05, max_maf = 0.95)
  stopifnot(nrow(vcf$gt1) == nrow(selection))
  sim1000G::startSimulation(vcf, totalNumberOfIndividuals = 4L * n_families)
  sim1000G::setRecombinationModel("poisson")

  father <- mother <- sibling1 <- sibling2 <- integer(n_families)
  for (i in seq_len(n_families)) {
    father[i] <- sim1000G::SIM$addUnrelatedIndividual()
    mother[i] <- sim1000G::SIM$addUnrelatedIndividual()
    sibling1[i] <- sim1000G::SIM$mate(father[i], mother[i])
    sibling2[i] <- sim1000G::SIM$mate(father[i], mother[i])
    if (i == 1L || i %% 500L == 0L) message("chr", chromosome, ": family ", i, "/", n_families)
  }

  variant <- data.frame(
    chromosome = as.integer(vcf$vcf[[1L]]),
    chromosome_variant_index = seq_len(nrow(vcf$vcf)),
    snp_id = paste0(vcf$vcf[[1L]], ":", vcf$vcf[[2L]]),
    position = as.integer(vcf$vcf[[2L]]),
    effect_allele = as.character(vcf$vcf[[5L]]),
    other_allele = as.character(vcf$vcf[[4L]])
  )
  selected <- match(variant$position, selection$position)
  stopifnot(!anyNA(selected), !anyDuplicated(variant$snp_id))
  variant$block_id <- selection$block_id[selected]
  variant$block_variant_index <- selection$block_variant_index[selected]
  family_id <- sprintf("F%06d", seq_len(n_families))
  dosage <- function(index) {
    g <- sim1000G::SIM$gt1[index, , drop = FALSE] + sim1000G::SIM$gt2[index, , drop = FALSE]
    storage.mode(g) <- "integer"
    colnames(g) <- variant$snp_id
    g
  }

  # Compare inherited founder labels, which retain phase through recombination.
  q1 <- sim1000G::SIM$origin1[sibling1, , drop = FALSE]
  q2 <- sim1000G::SIM$origin2[sibling1, , drop = FALSE]
  r1 <- sim1000G::SIM$origin1[sibling2, , drop = FALSE]
  r2 <- sim1000G::SIM$origin2[sibling2, , drop = FALSE]
  f1 <- sim1000G::SIM$origin1[father, 1L]
  f2 <- sim1000G::SIM$origin2[father, 1L]
  q1_paternal <- q1 == f1 | q1 == f2
  r1_paternal <- r1 == f1 | r1 == f2
  paternal <- 1L * (ifelse(q1_paternal, q1, q2) == ifelse(r1_paternal, r1, r2))
  maternal <- 1L * (ifelse(q1_paternal, q2, q1) == ifelse(r1_paternal, r2, r1))
  phase <- 1L * (q1 == r1) + 2L * (q1 == r2) + 4L * (q2 == r1) + 8L * (q2 == r2)
  ibd <- paternal + maternal
  dimnames(ibd) <- dimnames(phase) <- dimnames(paternal) <- dimnames(maternal) <-
    list(family_id, variant$snp_id)
  result <- list(chromosome = chromosome, family_id = family_id, variant = variant,
    father_genotype = dosage(father), mother_genotype = dosage(mother),
    sibling1_genotype = dosage(sibling1), sibling2_genotype = dosage(sibling2),
    ibd_state = ibd, ibd_phase_code = phase, paternal_ibd = paternal, maternal_ibd = maternal)
  stopifnot(all(phase %in% c(0L, 1L, 2L, 4L, 8L, 6L, 9L)),
            !any(ibd == 2L & result$sibling1_genotype != result$sibling2_genotype),
            !any(ibd == 1L & abs(result$sibling1_genotype - result$sibling2_genotype) == 2L))
  saveRDS(result, output)
  sim <- sim1000G::SIM
  for (name in c("gt1", "gt2", "origin1", "origin2", "population_gt1",
                 "population_gt2", "haplodata", "pool")) sim[[name]] <- NULL
  sim$individuals_generated <- 0L
  invisible(output)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  chromosomes <- if (length(args)) as.integer(args) else 18:22
  stopifnot(all(chromosomes %in% 18:22))
  for (chr in chromosomes) simulate_chromosome(chr)
}
