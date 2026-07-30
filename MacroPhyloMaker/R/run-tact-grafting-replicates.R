run_tact_grafting_replicates <- function(
    n_replicates,
    seed_start = 1,
    out_prefix,
    ...
) {
  seeds <- seq.int(seed_start, length.out = n_replicates)

  out <- vector("list", n_replicates)

  for (i in seq_len(n_replicates)) {
    rep_prefix <- sprintf("%s_rep%03d", out_prefix, i)

    message(
      "\n===== TACT replicate ",
      i,
      " / ",
      n_replicates,
      " | seed = ",
      seeds[i],
      " =====\n"
    )

    out[[i]] <- run_tact_grafting(
      out_prefix = rep_prefix,
      seed = seeds[i],
      ...
    )
  }

  manifest <- do.call(
    rbind,
    lapply(seq_along(out), function(i) {
      paths <- out[[i]]$paths

      data.frame(
        replicate = i,
        seed = seeds[i],
        cleaned_tree = paths$cleaned_tree,
        cleaned_tree_biogeo_labels = paths$cleaned_tree_biogeo_labels,
        backbone_biogeo_labels = paths$backbone_biogeo_labels,
        validation = paths$validation,
        taxonomy_mismatch_summary = paths$taxonomy_mismatch_summary,
        tree_not_taxonomy = paths$tree_not_taxonomy,
        taxonomy_not_tree = paths$taxonomy_not_tree,
        raw_tact_output_tree = paths$raw_tact_output_tree,
        stringsAsFactors = FALSE
      )
    })
  )

  manifest_file <- paste0(out_prefix, "_replicate_manifest.tsv")

  utils::write.table(
    manifest,
    manifest_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  structure(
    list(
      results = out,
      manifest = manifest,
      manifest_file = manifest_file
    ),
    class = "tact_grafting_replicates"
  )
}