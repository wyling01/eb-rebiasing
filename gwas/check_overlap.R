library(data.table)
library(ieugwasr)

run_howe_overlap <- function(
  opengwas_id = NULL,
  howe = NULL,
  method,
  procedure,
  pval = 1e-8,
  pop = "EUR",
  plink2 = Sys.getenv("PLINK2", unset = "plink2"),
  base_dir = "results/discovery_sets_bh",
  bfile = "ldref/1000G/EUR.chrpos",
  out_dir = "results/howe_overlap",
  ld_window_kb = 250,
  ld_window_r2 = 0.8
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ## 1. Get Howe/OpenGWAS significant SNPs
  ## Either pass a precomputed ieugwasr::tophits() output via `howe`,
  ## or pass `opengwas_id` and let the function call tophits().
  if (is.null(howe) && is.data.frame(opengwas_id)) {
    howe <- opengwas_id
    opengwas_id <- NULL
  }
  
  if (is.null(howe)) {
    if (is.null(opengwas_id)) {
      stop("Provide either `howe` as a tophits output or `opengwas_id`.")
    }
    
    howe <- tophits(
      id = opengwas_id,
      pval = pval,
      clump = 0,
      pop = pop
    )
  }

  setDT(howe)
  if (is.null(opengwas_id) && "id" %in% names(howe) && uniqueN(howe$id) == 1) {
    opengwas_id <- unique(howe$id)
  }
  
  chr_col <- intersect(c("chr", "chromosome", "CHR"), names(howe))[1]
  pos_col <- intersect(c("position", "pos", "bp", "BP"), names(howe))[1]
  
  if (is.na(chr_col) || is.na(pos_col)) {
    stop("`howe` must contain chromosome and position columns.")
  }
  
  howe[, chr_pos := paste0(get(chr_col), ":", get(pos_col))]
  howe_ids <- unique(howe$chr_pos)

  ## 2. Read your discovery SNPs as chr:pos
  disc_file <- file.path(
    base_dir,
    procedure,
    "chromosome_position",
    paste0(method, "_discoveries_chr_pos.txt")
  )

  disc <- fread(disc_file)
  disc[, chr_pos := paste0(chromosome, ":", pos)]
  disc_ids <- unique(disc$chr_pos)

  ## 3. Read local 1000G EUR reference panel
  bim <- fread(paste0(bfile, ".bim"), header = FALSE)
  setnames(bim, c("chr", "chr_pos", "cm", "bp", "a1", "a2"))

  ## 4. Restrict both sets to SNPs present in local reference
  disc_in_ref <- intersect(disc_ids, bim$chr_pos)
  howe_in_ref <- intersect(howe_ids, bim$chr_pos)

  exact_overlap <- intersect(disc_in_ref, howe_in_ref)

  cat("Discovery SNPs:", length(disc_ids), "\n")
  cat("Howe SNPs:", length(howe_ids), "\n")
  cat("Discovery in local 1000G EUR ref:", length(disc_in_ref), "\n")
  cat("Howe in local 1000G EUR ref:", length(howe_in_ref), "\n")
  cat("Exact overlap:", length(exact_overlap), "\n")

  ## 5. Write PLINK input SNP lists
  union_file <- file.path(out_dir, paste0("union_snps_", method, "_", procedure, ".txt"))
  disc_list_file <- file.path(out_dir, paste0("disc_snps_", method, "_", procedure, ".txt"))

  writeLines(unique(c(disc_in_ref, howe_in_ref)), union_file)
  writeLines(disc_in_ref, disc_list_file)

  ## 6. Run PLINK LD calculation
  out_prefix <- file.path(out_dir, paste0("howe_vs_", method, "_", procedure))

  plink_args <- c(
    "--bfile", bfile,
    "--extract", union_file,
    "--r2-unphased",
    "--ld-snp-list", disc_list_file,
    "--ld-window-kb", as.character(ld_window_kb),
    "--ld-window-r2", as.character(ld_window_r2),
    "--out", out_prefix
  )

  status <- system2(plink2, args = plink_args)

  if (status != 0) {
    stop("PLINK failed. Check log file: ", paste0(out_prefix, ".log"))
  }

  ## 7. Read PLINK LD output
  ld_file <- paste0(out_prefix, ".vcor")
  ld <- fread(ld_file)
  names(ld) <- sub("^#", "", names(ld))

  ## 8. Extract discovery-Howe LD pairs
  ## Expected orientation:
  ##   ID_A = discovery SNP, ID_B = Howe SNP
  ## But keep both orientations for robustness.
  ld_pairs <- rbindlist(list(
    ld[
      ID_A %in% disc_in_ref & ID_B %in% howe_in_ref,
      .(
        discovery_id = ID_A,
        howe_id = ID_B,
        r2 = UNPHASED_R2
      )
    ],
    ld[
      ID_B %in% disc_in_ref & ID_A %in% howe_in_ref,
      .(
        discovery_id = ID_B,
        howe_id = ID_A,
        r2 = UNPHASED_R2
      )
    ]
  ))

  ld_pairs <- unique(ld_pairs, by = c("discovery_id", "howe_id"))

  ## 9. Define matches
  ld_matched_discovery <- unique(
    ld_pairs[discovery_id != howe_id, discovery_id]
  )

  matched_discovery <- union(exact_overlap, ld_matched_discovery)
  unmatched_discovery <- setdiff(disc_in_ref, matched_discovery)
  disc_not_in_ref <- setdiff(disc_ids, disc_in_ref)

  ## 10. Summary table
  summary_table <- data.table(
    opengwas_id = ifelse(is.null(opengwas_id), NA_character_, opengwas_id),
    method = method,
    procedure = procedure,
    ld_window_kb = ld_window_kb,
    ld_window_r2 = ld_window_r2,
    category = c(
      "Discovery SNPs total",
      "Discovery SNPs in 1000G EUR ref",
      "Howe SNPs total",
      "Howe SNPs in 1000G EUR ref",
      "Discovery SNPs exact-matched to Howe",
      "Discovery SNPs non-exact LD-proxy matched to Howe",
      "Discovery SNPs exact or non-exact LD-proxy matched to Howe",
      "Discovery SNPs in ref but not matched to Howe",
      "Discovery SNPs not in 1000G EUR ref"
    ),
    count = c(
      length(disc_ids),
      length(disc_in_ref),
      length(howe_ids),
      length(howe_in_ref),
      length(exact_overlap),
      length(ld_matched_discovery),
      length(matched_discovery),
      length(unmatched_discovery),
      length(disc_not_in_ref)
    )
  )

  ## 11. Save outputs
  fwrite(summary_table, paste0(out_prefix, "_summary.tsv"), sep = "\t")
  fwrite(ld_pairs, paste0(out_prefix, "_ld_pairs.tsv"), sep = "\t")

  return(list(
    summary_table = summary_table,
    ld_pairs = ld_pairs,
    exact_overlap = exact_overlap,
    ld_matched_discovery = ld_matched_discovery,
    matched_discovery = matched_discovery,
    unmatched_discovery = unmatched_discovery,
    disc_not_in_ref = disc_not_in_ref,
    howe = howe
  ))
}
