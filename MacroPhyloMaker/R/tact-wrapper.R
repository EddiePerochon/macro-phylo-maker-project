
# tact_wrapper.R
# MacroPhyloMaker TACT integration wrapper
#
# This file prepares inputs for TACT and calls either a Docker-based or
# system-installed TACT executable. It does not distribute or vendor TACT.

# -----------------------------------------------------------------------------
# Logging helpers
# -----------------------------------------------------------------------------

#' Create a TACT workflow logger
#'
#' Creates a small timestamped log file next to the requested output prefix.
#'
#' @param out_prefix Character. Output prefix used by the TACT wrapper.
#' @param file_prefix Character. Prefix for the log filename.
#'
#' @return A logger object with a `file` element.
#' @keywords internal
.tact_make_logger <- function(out_prefix, file_prefix = "tact_grafting") {
  log_dir <- file.path(dirname(normalizePath(out_prefix, mustWork = FALSE)), "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  logfile <- file.path(
    log_dir,
    paste0(file_prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  cat("TACT grafting log\n", file = logfile)
  cat("Started: ", format(Sys.time()), "\n\n", file = logfile, append = TRUE, sep = "")
  structure(list(file = logfile), class = "tact_logger")
}

#' Write a message to console and logger
#'
#' Uses a package-level `log_msg()` if present, otherwise prints to console and
#' appends to the logger file.
#'
#' @param ... Components of the message.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns the message string.
#' @keywords internal
.tact_msg <- function(..., logger = NULL) {
  msg <- paste0(...)
  if (exists("log_msg", mode = "function", inherits = TRUE)) {
    get("log_msg", mode = "function", inherits = TRUE)(logger, msg)
  } else {
    cat(msg, "\n", sep = "")
    if (!is.null(logger) && !is.null(logger$file)) {
      cat(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n",
        file = logger$file, append = TRUE, sep = ""
      )
    }
  }
  invisible(msg)
}

#' Write a section heading to the TACT log
#'
#' @param title Character section title.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns `title`.
#' @keywords internal
.tact_section <- function(title, logger = NULL) {
  if (exists("log_section", mode = "function", inherits = TRUE)) {
    get("log_section", mode = "function", inherits = TRUE)(logger, title)
  } else {
    line <- paste(rep("-", nchar(title) + 4), collapse = "")
    .tact_msg(line, logger = logger)
    .tact_msg("| ", title, " |", logger = logger)
    .tact_msg(line, logger = logger)
  }
  invisible(title)
}

#' Require R packages used by the TACT wrapper
#'
#' @param pkgs Character vector of package names.
#'
#' @return Invisibly returns `TRUE` if all packages are available.
#' @keywords internal
.tact_require <- function(pkgs = c("ape", "stringr", "phytools")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Required R package(s) not installed: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Taxon-name helpers
# -----------------------------------------------------------------------------

#' Extract genus-like prefixes from tip labels
#'
#' Assumes tree labels use `Genus_epithet` style. Labels without underscores are
#' returned unchanged.
#'
#' @param labels Character vector of tip labels.
#'
#' @return Character vector of genus prefixes.
#' @keywords internal
.tact_get_genus <- function(labels) {
  has_us <- grepl("_", labels, fixed = TRUE)
  out <- labels
  out[has_us] <- sub("_.*$", "", labels[has_us])
  out
}

#' Extract epithets from tip labels
#'
#' @param labels Character vector of labels.
#'
#' @return Character vector of everything after the first underscore, or `NA`.
#' @keywords internal
.tact_get_epithet <- function(labels) {
  ifelse(grepl("_", labels, fixed = TRUE), sub("^[^_]+_", "", labels), NA_character_)
}

#' Detect genus-only labels
#'
#' @param labels Character vector of labels.
#'
#' @return Logical vector.
#' @keywords internal
.tact_is_genus_only <- function(labels) {
  !grepl("_", labels, fixed = TRUE)
}

#' Detect code-like species epithets
#'
#' Detects labels such as `Eburopone_MG04` or `Recurvidris_TH01` where the
#' epithet is uppercase letters followed by optional digits.
#'
#' @param labels Character vector of labels.
#'
#' @return Logical vector.
#' @keywords internal
.tact_is_code_epithet <- function(labels) {
  ep <- .tact_get_epithet(labels)
  !is.na(ep) & grepl("^[A-Z]{1,5}[0-9]{0,4}$", ep)
}

#' Convert underscores to spaces
#'
#' @param x Character vector.
#'
#' @return Character vector with all underscores converted to spaces.
#' @keywords internal
.tact_to_space <- function(x) gsub("_", " ", x, fixed = TRUE)

#' Convert spaces to underscores
#'
#' @param x Character vector.
#' @return Character vector with spaces converted to underscores.
#' @keywords internal
.tact_to_underscore <- function(x) gsub(" ", "_", x, fixed = TRUE)

#' Make a label unique among used labels
#'
#' @param x Proposed label.
#' @param used Character vector of already-used labels.
#'
#' @return Unique label.
#' @keywords internal
.tact_make_unique <- function(x, used) {
  if (!(x %in% used)) return(x)
  i <- 1L
  cand <- paste0(x, "_dup", i)
  while (cand %in% used) {
    i <- i + 1L
    cand <- paste0(x, "_dup", i)
  }
  cand
}

# -----------------------------------------------------------------------------
# Taxonomy readers/builders
# -----------------------------------------------------------------------------

#' Read and normalize taxonomy for TACT
#'
#' Reads AntWiki-style, simple, or TACT-style CSV taxonomy into a standard table
#' used by the wrapper. The returned table includes both a space-separated
#' species label and an underscore-separated species label.
#'
#' @param taxonomy Character path to a taxonomy file, or a data frame.
#' @param taxonomy_format One of `"antwiki"`, `"simple"`, or `"tact_csv"`.
#'   `"antwiki"` expects columns `TaxonName`, `Genus`, and `Species`.
#'   `"simple"` accepts columns `genus` and `species`, or `genus.species`.
#'   `"tact_csv"` expects columns `Family`, `genus`, and `genus.species`.
#' @param family Character family name to assign when absent.
#'
#' @return A data frame with columns `Family`, `genus`, `species`, and
#'   `species_underscore`.
#' @export
read_tact_taxonomy <- function(taxonomy,
                               taxonomy_format = c("antwiki", "simple", "tact_csv"),
                               family = "Formicidae") {
  taxonomy_format <- match.arg(taxonomy_format)

  if (is.data.frame(taxonomy)) {
    tx <- taxonomy
  } else if (taxonomy_format == "antwiki") {
    tx <- utils::read.table(
      taxonomy,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  } else {
    tx <- utils::read.table(
      taxonomy,
      header = TRUE,
      sep = if (grepl("\\.csv$", taxonomy, ignore.case = TRUE)) "," else "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fill = TRUE
    )
  }

  if (taxonomy_format == "antwiki") {
    req <- c("TaxonName", "Genus", "Species")
    miss <- setdiff(req, names(tx))
    if (length(miss)) {
      stop(
        "AntWiki taxonomy is missing required column(s): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }

    is_species_rank <- vapply(strsplit(trimws(tx$TaxonName), "\\s+"), length, integer(1)) == 2L
    sp <- tx[is_species_rank, c("Genus", "Species"), drop = FALSE]
    sp <- sp[nchar(sp$Genus) > 0 & nchar(sp$Species) > 0, , drop = FALSE]

    out <- data.frame(
      Family = family,
      genus = sp$Genus,
      species = paste(sp$Genus, sp$Species),
      species_underscore = paste(sp$Genus, sp$Species, sep = "_"),
      stringsAsFactors = FALSE
    )
  } else if (taxonomy_format == "tact_csv") {
    req <- c("Family", "genus", "genus.species")
    miss <- setdiff(req, names(tx))
    if (length(miss)) {
      stop(
        "TACT CSV taxonomy is missing required column(s): ", paste(miss, collapse = ", "),
        call. = FALSE
      )
    }

    out <- data.frame(
      Family = tx$Family,
      genus = tx$genus,
      species = tx[["genus.species"]],
      species_underscore = .tact_to_underscore(tx[["genus.species"]]),
      stringsAsFactors = FALSE
    )
  } else {
    if (!"Family" %in% names(tx)) tx$Family <- family
    if (!"genus" %in% names(tx) && "Genus" %in% names(tx)) tx$genus <- tx$Genus

    if (!"species" %in% names(tx)) {
      if ("genus.species" %in% names(tx)) {
        tx$species <- tx[["genus.species"]]
      } else if (all(c("genus", "Species") %in% names(tx))) {
        tx$species <- paste(tx$genus, tx$Species)
      } else {
        stop("Simple taxonomy needs columns genus and species, or genus.species.", call. = FALSE)
      }
    }

    out <- data.frame(
      Family = tx$Family,
      genus = tx$genus,
      species = tx$species,
      species_underscore = .tact_to_underscore(tx$species),
      stringsAsFactors = FALSE
    )
  }

  out <- out[!is.na(out$genus) & !is.na(out$species_underscore), , drop = FALSE]
  out <- out[nchar(out$genus) > 0 & nchar(out$species_underscore) > 0, , drop = FALSE]
  out <- out[!duplicated(out$species_underscore), , drop = FALSE]
  out <- out[order(out$genus, out$species_underscore), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Write a TACT rank CSV
#'
#' Writes taxonomy as a rank-ordered CSV for `tact_build_taxonomic_tree`. The
#' species column uses spaces rather than underscores because DendroPy treats
#' unquoted Newick underscores as spaces when reading the backbone and taxonomy
#' tree. For example, `TACTTMPProcryptocerus01_carbonarius` is written as
#' `TACTTMPProcryptocerus01 carbonarius`.
#'
#' @param tax Data frame produced by [read_tact_taxonomy()] or by the wrapper's
#'   preprocessing steps.
#' @param file Character output CSV path.
#' @param family Character family name. Present for API consistency; the value
#'   in `tax$Family` is used.
#'
#' @return Invisibly returns `file`.
#' @export
write_tact_taxonomy_csv <- function(tax, file, family = "Formicidae") {
  species_label <- if ("species" %in% names(tax)) {
    tax$species
  } else {
    tax$species_underscore
  }

  species_label <- .tact_to_space(species_label)

  out <- data.frame(
    Family = tax$Family,
    genus = tax$genus,
    species = species_label,
    stringsAsFactors = FALSE
  )

  out <- out[
    !is.na(out$Family) & !is.na(out$genus) & !is.na(out$species),
    ,
    drop = FALSE
  ]
  out <- out[
    nchar(out$Family) > 0 & nchar(out$genus) > 0 & nchar(out$species) > 0,
    ,
    drop = FALSE
  ]
  out <- out[!duplicated(out$species), , drop = FALSE]
  out <- out[order(out$Family, out$genus, out$species), , drop = FALSE]

  utils::write.table(
    out,
    file = file,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )

  invisible(file)
}

#' Write a simple taxonomy Newick tree
#'
#' This function is retained for debugging and fallback use. The recommended
#' production workflow is to write a CSV with [write_tact_taxonomy_csv()] and
#' build the taxonomy tree using TACT's `tact_build_taxonomic_tree` command.
#'
#' @param tax Data frame with taxonomy.
#' @param file Character output Newick path.
#' @param family Character root label.
#'
#' @return Invisibly returns `file`.
#' @keywords internal
write_tact_taxonomy_newick <- function(tax, file, family = "Formicidae") {
  tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  genus_groups <- split(tax$species_underscore, tax$genus)
  genus_strings <- vapply(names(genus_groups), function(g) {
    tips <- genus_groups[[g]]
    if (length(tips) == 1L) {
      paste0("(", tips, ")", g)
    } else {
      paste0("(", paste(tips, collapse = ","), ")", g)
    }
  }, character(1))
  nwk <- paste0("(", paste(genus_strings, collapse = ","), ")", family, ";")
  writeLines(nwk, file)
  invisible(file)
}

#' Write the cleaned TACT tree from a TACT grafting result
#'
#' @param res Result object returned by `run_tact_grafting()`.
#' @param file Output Newick file path.
#'
#' @return Invisibly returns `file`.
#' @export
write_tact_result_tree <- function(res, file) {
  if (is.null(res$tree) || !inherits(res$tree, "phylo")) {
    stop("`res$tree` is not a valid phylo object.", call. = FALSE)
  }

  ape::write.tree(
    phy = res$tree,
    file = file
  )

  invisible(file)
}

# -----------------------------------------------------------------------------
# Exclusion utilities
# -----------------------------------------------------------------------------

#' Resolve taxa and MRCA clades to exclude from TACT grafting
#'
#' Exclusions are used in two ways: tips present in the backbone are protected
#' from receiving missing taxa during non-monophyletic genus splitting; taxa not
#' present in the backbone can be removed from the TACT taxonomy by
#' `.remove_excluded_taxa_not_in_backbone()`.
#'
#' @param tree An object of class `phylo`.
#' @param exclude_taxa Optional character vector of exact tip labels to protect.
#' @param exclude_mrca Optional character vector or list of character vectors.
#'   Each vector is resolved to an MRCA and all descendant tips are protected.
#' @param exclude_missing One of `"warn"`, `"error"`, or `"ignore"`, controlling
#'   behavior when an excluded tip is absent from the backbone.
#'
#' @return Character vector of backbone tip labels to protect.
#' @export
resolve_tact_exclusions <- function(tree,
                                    exclude_taxa = NULL,
                                    exclude_mrca = NULL,
                                    exclude_missing = c("warn", "error", "ignore")) {
  exclude_missing <- match.arg(exclude_missing)
  out <- character(0)

  if (!is.null(exclude_taxa) && length(exclude_taxa)) {
    exclude_taxa <- unique(.tact_to_underscore(exclude_taxa))
    missing <- setdiff(exclude_taxa, tree$tip.label)
    if (length(missing)) {
      msg <- paste("Excluded tip(s) not found in backbone:", paste(missing, collapse = ", "))
      if (exclude_missing == "error") stop(msg, call. = FALSE)
      if (exclude_missing == "warn") warning(msg, call. = FALSE)
    }
    out <- union(out, intersect(exclude_taxa, tree$tip.label))
  }

  if (!is.null(exclude_mrca) && length(exclude_mrca)) {
    sets <- if (is.list(exclude_mrca)) exclude_mrca else list(exclude_mrca)
    for (anchors in sets) {
      anchors <- .tact_to_underscore(anchors)
      missing <- setdiff(anchors, tree$tip.label)
      if (length(missing)) {
        stop("MRCA exclusion anchor(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
      }

      node <- if (length(anchors) == 1L) {
        match(anchors, tree$tip.label)
      } else {
        ape::getMRCA(tree, anchors)
      }

      if (is.null(node) || is.na(node)) {
        stop("Could not resolve MRCA exclusion for: ", paste(anchors, collapse = ", "), call. = FALSE)
      }

      desc <- if (node <= ape::Ntip(tree)) node else phytools::getDescendants(tree, node)
      desc <- desc[desc <= ape::Ntip(tree)]
      out <- union(out, tree$tip.label[desc])
    }
  }

  sort(unique(out))
}

#' Remove excluded taxa that are absent from the backbone from taxonomy
#'
#' @param tax TACT taxonomy data frame.
#' @param exclude_taxa Character vector of taxa requested for exclusion.
#' @param backbone_tips Character vector of current backbone tip labels.
#' @param logger Optional logger object.
#'
#' @return A list with `taxonomy` and `removed` elements.
#' @keywords internal
.remove_excluded_taxa_not_in_backbone <- function(tax, exclude_taxa, backbone_tips, logger = NULL) {
  empty <- data.frame(
    species_underscore = character(),
    species = character(),
    genus = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )

  if (is.null(exclude_taxa) || !length(exclude_taxa)) {
    return(list(taxonomy = tax, removed = empty))
  }

  exclude_taxa <- unique(.tact_to_underscore(exclude_taxa))
  absent_from_backbone <- setdiff(exclude_taxa, backbone_tips)
  if (!length(absent_from_backbone)) {
    return(list(taxonomy = tax, removed = empty))
  }

  hit <- tax$species_underscore %in% absent_from_backbone
  removed <- tax[hit, , drop = FALSE]
  removed_log <- data.frame(
    species_underscore = removed$species_underscore,
    species = removed$species,
    genus = removed$genus,
    reason = "excluded_taxon_absent_from_backbone_removed_from_tact_taxonomy",
    stringsAsFactors = FALSE
  )

  if (nrow(removed_log)) {
    .tact_msg(
      "Excluded taxa removed from TACT taxonomy because absent from backbone: ",
      paste(removed_log$species_underscore, collapse = ", "),
      logger = logger
    )
  }

  list(taxonomy = tax[!hit, , drop = FALSE], removed = removed_log)
}

# -----------------------------------------------------------------------------
# Monophyly diagnostic utilities
# -----------------------------------------------------------------------------

#' Detect non-monophyletic genera for TACT preprocessing
#'
#' Uses Eddie's `nonmono_genera()` helper if it is available. Otherwise it falls
#' back to `MonoPhy::AssessMonophyly()`.
#'
#' @param tree An object of class `phylo`.
#' @param exclude_tips Optional tips to drop before testing. For TACT splitting,
#'   this should usually be `character(0)` so excluded tips can still trigger
#'   pseudo-genus creation.
#' @param genus_sep Character genus separator passed to Eddie's helper.
#'
#' @return Character vector of non-monophyletic genera.
#' @export
#' @importFrom utils capture.output
#' @examples
#' \dontrun{
#' nonmono <- detect_nonmono_genera_for_tact(tree)
#' }
detect_nonmono_genera_for_tact <- function(tree, exclude_tips = character(0), genus_sep = "_") {
  tr <- tree
  if (length(exclude_tips)) {
    tr <- ape::drop.tip(tr, intersect(exclude_tips, tr$tip.label))
  }

  if (exists("nonmono_genera", mode = "function", inherits = TRUE)) {
    return(get("nonmono_genera", mode = "function", inherits = TRUE)(tr, genus_sep = genus_sep))
  }

  if (!requireNamespace("MonoPhy", quietly = TRUE)) {
    stop("Package MonoPhy is required unless Eddie's nonmono_genera() is available.", call. = FALSE)
  }

  genera <- .tact_get_genus(tr$tip.label)
  tab <- table(genera)
  multi <- names(tab)[tab >= 2]
  tax <- data.frame(species = tr$tip.label, genus = genera, stringsAsFactors = FALSE)
  ass <- MonoPhy::AssessMonophyly(tr, tax)
  res <- MonoPhy::GetResultMonophyly(ass)
  if (!"Taxa" %in% names(res)) return(character(0))

  status_col <- intersect(c("Monophyly", "Status", "monophyly"), names(res))[1]
  if (is.na(status_col)) return(character(0))

  sort(intersect(res$Taxa[!tolower(res[[status_col]]) %in% "monophyletic"], multi))
}

#' Return descendant tip labels from a node
#'
#' @param tree An object of class `phylo`.
#' @param node Numeric node number.
#'
#' @return Character vector of descendant tip labels.
#' @keywords internal
.tact_desc_tips <- function(tree, node) {
  if (node <= ape::Ntip(tree)) return(tree$tip.label[node])
  dd <- phytools::getDescendants(tree, node)
  tree$tip.label[dd[dd <= ape::Ntip(tree)]]
}

#' Partition a non-monophyletic genus by MRCA child branches
#'
#' @param tree An object of class `phylo`.
#' @param genus Character genus name.
#' @param exclude_tips Character vector of protected tips that should not receive
#'   missing species.
#'
#' @return A list of character vectors, each containing graftable tips for one
#'   temporary pseudo-genus.
#' @keywords internal
.partition_genus_by_mrca_children <- function(tree, genus, exclude_tips = character(0)) {
  all_tips <- grep(paste0("^", genus, "_"), tree$tip.label, value = TRUE)
  graftable <- setdiff(all_tips, exclude_tips)
  if (length(graftable) == 0L) return(list())
  if (length(graftable) == 1L) return(list(graftable))

  mrca <- ape::getMRCA(tree, graftable)
  children <- tree$edge[tree$edge[, 1] == mrca, 2]
  groups <- list()
  for (ch in children) {
    tips <- intersect(.tact_desc_tips(tree, ch), graftable)
    if (length(tips)) groups[[length(groups) + 1L]] <- tips
  }

  if (length(groups) <= 1L) groups <- as.list(graftable)
  groups
}

# -----------------------------------------------------------------------------
# Backbone preparation
# -----------------------------------------------------------------------------

#' Prepare backbone labels for TACT
#'
#' Handles genus-only tips, code-like epithets, and backbone tips that are absent
#' from the taxonomic table. Code-like tips are kept by default and added to the
#' taxonomy as known terminals so TACT does not treat them as missing.
#'
#' @param tree An object of class `phylo`.
#' @param tax Taxonomy data frame from [read_tact_taxonomy()].
#' @param genus_only One of `"replace_random_species"`, `"temporary_species"`,
#'   or `"keep"`.
#' @param species_code One of `"keep"`, `"temporary_species"`, or `"drop"`.
#' @param sample_with_replacement Logical. Whether genus-only tips can sample an
#'   already-used species if no unused valid species remain.
#' @param fallback_prefix Character prefix for temporary species labels.
#' @param seed Integer random seed.
#' @param logger Optional logger object.
#'
#' @return A list with `tree`, `taxonomy`, `label_map`, and `dropped` elements.
#' @export
prepare_tact_backbone_labels <- function(tree,
                                         tax,
                                         genus_only = c("replace_random_species", "temporary_species", "keep"),
                                         species_code = c("keep", "temporary_species", "drop"),
                                         sample_with_replacement = TRUE,
                                         fallback_prefix = "sp",
                                         seed = 1,
                                         logger = NULL) {
  genus_only <- match.arg(genus_only)
  species_code <- match.arg(species_code)
  set.seed(seed)

  labels <- tree$tip.label
  used <- labels
  map <- data.frame(old = character(), new = character(), reason = character(), stringsAsFactors = FALSE)
  dropped <- character(0)
  species_by_genus <- split(tax$species_underscore, tax$genus)

  idx_go <- which(.tact_is_genus_only(labels))
  if (length(idx_go)) {
    .tact_msg("Genus-only backbone tips detected: ", length(idx_go), logger = logger)
    for (i in idx_go) {
      g <- labels[i]
      if (genus_only == "keep") next

      if (genus_only == "replace_random_species") {
        pool <- setdiff(species_by_genus[[g]], used)
        if (length(pool) == 0L && isTRUE(sample_with_replacement)) pool <- species_by_genus[[g]]
        if (is.null(pool) || length(pool) == 0L) {
          proposal <- paste0(g, "_", fallback_prefix, sprintf("%03d", i))
          reason <- "genus_only_fallback"
        } else {
          proposal <- sample(pool, 1L)
          reason <- "genus_only_random_valid_species"
        }
      } else {
        proposal <- paste0(g, "_", fallback_prefix, sprintf("%03d", i))
        reason <- "genus_only_temporary_species"
      }

      new <- .tact_make_unique(proposal, used)
      map <- rbind(map, data.frame(old = labels[i], new = new, reason = reason, stringsAsFactors = FALSE))
      labels[i] <- new
      used <- c(used, new)
    }
  }

  idx_code <- which(.tact_is_code_epithet(labels))
  if (length(idx_code)) {
    .tact_msg("Species-code-like backbone tips detected: ", length(idx_code), logger = logger)
  }

  if (length(idx_code) && species_code == "drop") {
    dropped <- labels[idx_code]
    labels <- labels[-idx_code]
    tree <- ape::drop.tip(tree, dropped)
  } else if (length(idx_code) && species_code == "temporary_species") {
    for (i in idx_code) {
      g <- .tact_get_genus(labels[i])
      ep <- .tact_get_epithet(labels[i])
      new <- paste0(g, "_", fallback_prefix, ep)
      new <- .tact_make_unique(new, setdiff(labels, labels[i]))
      map <- rbind(map, data.frame(
        old = labels[i], new = new,
        reason = "species_code_temporary_species",
        stringsAsFactors = FALSE
      ))
      labels[i] <- new
    }
  }

  tree$tip.label <- labels

  missing_from_tax <- setdiff(tree$tip.label, tax$species_underscore)
  if (length(missing_from_tax)) {
    add <- data.frame(
      Family = tax$Family[1],
      genus = .tact_get_genus(missing_from_tax),
      species = .tact_to_space(missing_from_tax),
      species_underscore = missing_from_tax,
      stringsAsFactors = FALSE
    )
    tax <- rbind(tax, add)
    tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  }

  list(tree = tree, taxonomy = tax, label_map = map, dropped = dropped)
}

# -----------------------------------------------------------------------------
# Nonmonophyletic genus splitting
# -----------------------------------------------------------------------------

#' Split non-monophyletic genera into temporary pseudo-genera for TACT
#'
#' Non-monophyletic genera cannot be safely represented as a single taxonomic
#' unit for TACT grafting. This helper temporarily renames separate backbone
#' occurrences of a non-monophyletic genus as pseudo-genera such as
#' `TACTTMPNeoponera01`. Protected excluded tips are renamed to `TACTEXCL...`
#' pseudo-genera and receive no missing species. After TACT, temporary labels are
#' removed by [restore_tact_temp_names()].
#'
#' @param tree An object of class `phylo`.
#' @param tax Taxonomy data frame.
#' @param nonmono One of `"split"`, `"skip"`, or `"error"`.
#' @param nonmono_allocation Missing-species allocation among pseudo-genera:
#'   `"proportional"`, `"equal"`, or `"random"`.
#' @param exclude_tips Character vector of backbone tips to protect from grafting.
#' @param seed Integer random seed.
#' @param logger Optional logger object.
#'
#' @return A list with `tree`, `taxonomy`, `map`, `nonmono`, and `skipped`.
#' @export
split_nonmono_genera_for_tact <- function(tree,
                                          tax,
                                          nonmono = c("split", "skip", "error"),
                                          nonmono_allocation = c("proportional", "equal", "random"),
                                          exclude_tips = character(0),
                                          seed = 1,
                                          logger = NULL) {
  nonmono <- match.arg(nonmono)
  nonmono_allocation <- match.arg(nonmono_allocation)
  set.seed(seed)

  nonmono_g <- detect_nonmono_genera_for_tact(tree, exclude_tips = character(0))
  excluded_genera <- unique(.tact_get_genus(exclude_tips))
  nonmono_g <- sort(unique(c(
    nonmono_g,
    intersect(excluded_genera, unique(.tact_get_genus(tree$tip.label)))
  )))

  if (!length(nonmono_g)) {
    return(list(tree = tree, taxonomy = tax, map = data.frame(), nonmono = data.frame(), skipped = data.frame()))
  }

  .tact_msg(
    "Non-monophyletic genera detected for TACT handling: ",
    paste(nonmono_g, collapse = ", "),
    logger = logger
  )

  if (nonmono == "error") {
    stop("Non-monophyletic genera present: ", paste(nonmono_g, collapse = ", "), call. = FALSE)
  }

  map <- data.frame(
    original_genus = character(), temp_genus = character(), old_label = character(),
    new_label = character(), cluster = integer(), protected = logical(),
    stringsAsFactors = FALSE
  )
  skipped <- data.frame(
    genus = character(), species_underscore = character(), reason = character(),
    stringsAsFactors = FALSE
  )
  labels <- tree$tip.label

  for (g in nonmono_g) {
    species_in_tax <- tax$species_underscore[tax$genus == g]
    species_in_tree <- grep(paste0("^", g, "_"), labels, value = TRUE)
    species_missing <- setdiff(species_in_tax, species_in_tree)
    protected <- intersect(species_in_tree, exclude_tips)

    if (nonmono == "skip") {
      if (length(species_missing)) {
        skipped <- rbind(skipped, data.frame(
          genus = g,
          species_underscore = species_missing,
          reason = "nonmonophyletic_genus_skipped",
          stringsAsFactors = FALSE
        ))
      }
      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      next
    }

    groups <- .partition_genus_by_mrca_children(tree, g, exclude_tips = exclude_tips)
    if (!length(groups)) {
      if (length(species_missing)) {
        skipped <- rbind(skipped, data.frame(
          genus = g,
          species_underscore = species_missing,
          reason = "only_excluded_branches_available",
          stringsAsFactors = FALSE
        ))
      }
      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      next
    }

    temp_genera <- paste0("TACTTMP", g, sprintf("%02d", seq_along(groups)))

    for (k in seq_along(groups)) {
      old <- groups[[k]]
      new <- sub(paste0("^", g, "_"), paste0(temp_genera[k], "_"), old)
      labels[match(old, labels)] <- new
      map <- rbind(map, data.frame(
        original_genus = g,
        temp_genus = temp_genera[k],
        old_label = old,
        new_label = new,
        cluster = k,
        protected = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    if (length(protected)) {
      temp_ex <- paste0("TACTEXCL", g, sprintf("%02d", seq_along(protected)))
      for (j in seq_along(protected)) {
        old <- protected[j]
        new <- sub(paste0("^", g, "_"), paste0(temp_ex[j], "_"), old)
        labels[match(old, labels)] <- new
        map <- rbind(map, data.frame(
          original_genus = g,
          temp_genus = temp_ex[j],
          old_label = old,
          new_label = new,
          cluster = NA_integer_,
          protected = TRUE,
          stringsAsFactors = FALSE
        ))
      }
    }

    # Remove original rows for all existing backbone tips of this genus, then add
    # rows for the temporary labels now present in the backbone.
    tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_in_tree), , drop = FALSE]
    g_map <- map[map$original_genus == g, , drop = FALSE]
    if (nrow(g_map)) {
      add_existing <- data.frame(
        Family = tax$Family[1],
        genus = g_map$temp_genus,
        species = .tact_to_space(g_map$new_label),
        species_underscore = g_map$new_label,
        stringsAsFactors = FALSE
      )
      tax <- rbind(tax, add_existing)
    }

    if (length(species_missing)) {
      weights <- vapply(groups, length, integer(1))
      if (nonmono_allocation == "equal") weights <- rep(1, length(groups))

      if (nonmono_allocation == "random") {
        assign <- sample(seq_along(groups), length(species_missing), replace = TRUE)
      } else {
        assign <- sample(seq_along(groups), length(species_missing), replace = TRUE, prob = weights)
      }

      tax <- tax[!(tax$genus == g & tax$species_underscore %in% species_missing), , drop = FALSE]
      missing_epithet <- sub(paste0("^", g, "_"), "", species_missing)
      missing_temp_labels <- paste0(temp_genera[assign], "_", missing_epithet)

      add_missing <- data.frame(
        Family = tax$Family[1],
        genus = temp_genera[assign],
        species = .tact_to_space(missing_temp_labels),
        species_underscore = missing_temp_labels,
        stringsAsFactors = FALSE
      )
      tax <- rbind(tax, add_missing)
    }
  }

  tree$tip.label <- labels
  tax <- tax[!duplicated(tax$species_underscore), , drop = FALSE]
  nonmono_log <- data.frame(genus = nonmono_g, action = nonmono, stringsAsFactors = FALSE)

  list(tree = tree, taxonomy = tax, map = map, nonmono = nonmono_log, skipped = skipped)
}

#' Restore temporary TACT pseudo-genus names
#'
#' Converts temporary labels such as `TACTTMPNeoponera01_apicalis` and
#' `TACTEXCLNeoponera01_bucki` back to `Neoponera_apicalis` and
#' `Neoponera_bucki`.
#'
#' @param tree An object of class `phylo` or `multiPhylo`.
#' @param map Mapping data frame returned by [split_nonmono_genera_for_tact()].
#'
#' @return A tree object of the same class with temporary prefixes removed.
#' @export
restore_tact_temp_names <- function(tree, map) {
  if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, restore_tact_temp_names, map = map)
    class(out) <- "multiPhylo"
    return(out)
  }

  if (is.null(map) || !nrow(map)) return(tree)
  labels <- tree$tip.label
  for (i in seq_len(nrow(map))) {
    tg <- map$temp_genus[i]
    og <- map$original_genus[i]
    labels <- sub(paste0("^", tg, "_"), paste0(og, "_"), labels)
  }
  tree$tip.label <- labels
  tree
}

#' Normalize TACT output labels from spaces to underscores
#'
#' TACT/DendroPy may emit labels with spaces because underscores in unquoted
#' Newick labels are interpreted as spaces. This helper supports both `phylo`
#' and `multiPhylo` outputs.
#'
#' @param tree A `phylo` or `multiPhylo` object.
#'
#' @return A tree object of the same class with spaces converted to underscores
#'   in tip labels.
#' @keywords internal
.tact_normalize_output_labels <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, .tact_normalize_output_labels)
    class(out) <- "multiPhylo"
    return(out)
  }
  tree$tip.label <- .tact_to_underscore(tree$tip.label)
  tree
}

#' Remove invalid node labels before writing TACT output
#'
#' Some TACT Nexus outputs contain node labels that do not match the number of
#' internal nodes after reading with `ape`. `ape::write.tree()` warns in that
#' case. These labels are not used by the MacroPhyloMaker TACT wrapper, so this
#' helper drops them before writing cleaned Newick trees. Supports both `phylo`
#' and `multiPhylo` objects.
#'
#' @param tree A `phylo` or `multiPhylo` object.
#'
#' @return A tree object of the same class with invalid `node.label` removed.
#' @keywords internal
.tact_drop_invalid_node_labels <- function(tree) {
  if (inherits(tree, "multiPhylo")) {
    out <- lapply(tree, .tact_drop_invalid_node_labels)
    class(out) <- "multiPhylo"
    return(out)
  }

  if (!is.null(tree$node.label) && length(tree$node.label) != ape::Nnode(tree)) {
    tree$node.label <- NULL
  }
  tree
}

#' Read a TACT output tree
#'
#' Reads Newick outputs with `ape::read.tree()` and Nexus outputs with
#' `ape::read.nexus()`. TACT can write more than one tree, so this may return a
#' `phylo` or `multiPhylo` object.
#'
#' @param path Character path to a TACT output tree.
#'
#' @return A `phylo` or `multiPhylo` object.
#' @keywords internal
.tact_read_output_tree <- function(path) {
  if (grepl("\\.(nex|nexus)(\\.tre)?$", path, ignore.case = TRUE) ||
      grepl("\\.nexus\\.tre$", path, ignore.case = TRUE)) {
    return(ape::read.nexus(path))
  }
  ape::read.tree(path)
}

# -----------------------------------------------------------------------------
# TACT execution
# -----------------------------------------------------------------------------

#' Run `tact_build_taxonomic_tree`
#'
#' @param work_dir Working directory mounted into Docker or used locally.
#' @param taxonomy_csv TACT rank CSV path.
#' @param taxonomy_tree Output taxonomy tree path.
#' @param tact_runner One of `"docker"` or `"system"`.
#' @param docker_image Docker image name.
#' @param tact_build_bin Local executable for taxonomy-tree construction.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list with stdout, stderr, status, and taxonomy path.
#' @export
run_tact_build_taxonomy_external <- function(work_dir,
                                             taxonomy_csv,
                                             taxonomy_tree,
                                             tact_runner = c("docker", "system"),
                                             docker_image = "jonchang/tact",
                                             tact_build_bin = "tact_build_taxonomic_tree",
                                             logger = NULL) {
  tact_runner <- match.arg(tact_runner)

  stdout_file <- file.path(work_dir, "tact_build_taxonomy_stdout.log")
  stderr_file <- file.path(work_dir, "tact_build_taxonomy_stderr.log")

  if (tact_runner == "docker") {
    args <- c(
      "run", "--rm",
      "-v", paste0(normalizePath(work_dir), ":/workdir"),
      "-w", "/workdir",
      docker_image,
      "tact_build_taxonomic_tree",
      basename(taxonomy_csv),
      "--output", basename(taxonomy_tree)
    )

    .tact_msg("TACT taxonomy-build Docker command: docker ", paste(args, collapse = " "), logger = logger)
    status <- system2("docker", args = args, stdout = stdout_file, stderr = stderr_file)
  } else {
    args <- c(taxonomy_csv, "--output", taxonomy_tree)
    .tact_msg("TACT taxonomy-build command: ", tact_build_bin, " ", paste(args, collapse = " "), logger = logger)
    status <- system2(tact_build_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  }

  if (!identical(status, 0L)) {
    err <- if (file.exists(stderr_file)) {
      paste(tail(readLines(stderr_file, warn = FALSE), 80), collapse = "\n")
    } else {
      "<no stderr>"
    }
    stop(
      "TACT taxonomy-tree construction failed with exit status ", status,
      ".\nStderr tail:\n", err,
      "\nFull stderr: ", stderr_file,
      call. = FALSE
    )
  }

  invisible(list(stdout = stdout_file, stderr = stderr_file, status = status, taxonomy_tree = taxonomy_tree))
}

#' Run `tact_add_taxa`
#'
#' Calls TACT using either Docker or a system-installed executable.
#'
#' @param work_dir Working directory mounted into Docker or used locally.
#' @param backbone_file TACT-ready backbone tree.
#' @param taxonomy_file TACT taxonomy tree produced by `tact_build_taxonomic_tree`.
#' @param output_prefix TACT output prefix.
#' @param outgroups Optional outgroup labels to pass to TACT.
#' @param tact_runner One of `"docker"` or `"system"`.
#' @param docker_image Docker image name.
#' @param tact_bin Local `tact_add_taxa` executable.
#' @param extra_args Additional command-line arguments for TACT.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list with stdout, stderr, and status.
#' @export
run_tact_external <- function(work_dir,
                              backbone_file,
                              taxonomy_file,
                              output_prefix,
                              outgroups = NULL,
                              tact_runner = c("docker", "system"),
                              docker_image = "jonchang/tact",
                              tact_bin = "tact_add_taxa",
                              extra_args = character(0),
                              logger = NULL) {
  tact_runner <- match.arg(tact_runner)
  stdout_file <- file.path(work_dir, "tact_stdout.log")
  stderr_file <- file.path(work_dir, "tact_stderr.log")

  if (tact_runner == "docker") {
    args <- c(
      "run", "--rm",
      "-v", paste0(normalizePath(work_dir), ":/workdir"),
      "-w", "/workdir",
      docker_image,
      "tact_add_taxa",
      "--backbone", basename(backbone_file),
      "--taxonomy", basename(taxonomy_file),
      "--output", basename(output_prefix)
    )
    if (!is.null(outgroups) && length(outgroups) && all(nzchar(outgroups))) {
      args <- c(args, "--outgroups", paste(outgroups, collapse = ","))
    }
    args <- c(args, extra_args)
    .tact_msg("TACT Docker command: docker ", paste(args, collapse = " "), logger = logger)
    status <- system2("docker", args = args, stdout = stdout_file, stderr = stderr_file)
  } else {
    args <- c(
      "--backbone", backbone_file,
      "--taxonomy", taxonomy_file,
      "--output", output_prefix
    )
    if (!is.null(outgroups) && length(outgroups) && all(nzchar(outgroups))) {
      args <- c(args, "--outgroups", paste(outgroups, collapse = ","))
    }
    args <- c(args, extra_args)
    .tact_msg("TACT command: ", tact_bin, " ", paste(args, collapse = " "), logger = logger)
    status <- system2(tact_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  }

  if (!identical(status, 0L)) {
    err <- if (file.exists(stderr_file)) {
      paste(tail(readLines(stderr_file, warn = FALSE), 80), collapse = "\n")
    } else {
      "<no stderr>"
    }
    stop(
      "TACT failed with exit status ", status,
      ".\nStderr tail:\n", err,
      "\nFull stderr: ", stderr_file,
      call. = FALSE
    )
  }

  invisible(list(stdout = stdout_file, stderr = stderr_file, status = status))
}

#' Find a TACT output tree in a work directory
#'
#' TACT commonly writes `{output}.newick.tre` and `{output}.nexus.tre`.
#'
#' @param work_dir Character work directory.
#' @param output_prefix Character output prefix passed to TACT.
#'
#' @return Character path or `NA_character_`.
#' @keywords internal
.find_tact_output_tree <- function(work_dir, output_prefix) {
  stem <- basename(output_prefix)

  candidates <- c(
    file.path(work_dir, paste0(stem, ".newick.tre")),
    file.path(work_dir, paste0(stem, ".nexus.tre")),
    file.path(work_dir, paste0(stem, ".tre")),
    file.path(work_dir, paste0(stem, ".nwk")),
    file.path(work_dir, paste0(stem, ".newick"))
  )

  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) return(NA_character_)
  candidates[which.max(file.info(candidates)$mtime)]
}

# -----------------------------------------------------------------------------
# Main wrapper
# -----------------------------------------------------------------------------

#' Run TACT taxonomy-based grafting for a MacroPhyloMaker backbone
#'
#' Prepares a backbone tree and taxonomy for TACT, optionally handles
#' non-monophyletic genera by splitting them into temporary pseudo-genera,
#' protects excluded taxa or clades from receiving grafted species, runs TACT via
#' Docker or a system executable, and restores temporary names in the final tree.
#'
#' TACT itself is not distributed with MacroPhyloMaker. This function only
#' prepares inputs and calls external TACT commands.
#'
#' @param backbone_tree Character path to a Newick backbone tree, or a `phylo`
#'   object.
#' @param taxonomy Character path to a taxonomy file, or a data frame.
#' @param out_prefix Character output prefix. The wrapper writes audit files and
#'   a TACT work directory using this prefix.
#' @param family Character family/root label for taxonomy rows.
#' @param taxonomy_format One of `"antwiki"`, `"simple"`, or `"tact_csv"`.
#' @param taxonomy_output Currently retained for compatibility. The recommended
#'   path always writes CSV and lets TACT build the taxonomy tree.
#' @param tact_runner One of `"docker"`, `"system"`, or `"none"`. Use `"none"`
#'   to prepare files without running TACT.
#' @param docker_image Docker image name for TACT.
#' @param tact_bin Local executable name or path for `tact_add_taxa`.
#' @param outgroups Optional character vector of outgroup taxa to pass to TACT.
#' @param seed Integer random seed for stochastic preprocessing decisions.
#' @param nonmono One of `"split"`, `"skip"`, or `"error"`, controlling treatment
#'   of non-monophyletic genera.
#' @param nonmono_allocation One of `"proportional"`, `"equal"`, or `"random"`,
#'   controlling allocation of missing species among pseudo-genera.
#' @param genus_only One of `"replace_random_species"`, `"temporary_species"`, or
#'   `"keep"`, controlling genus-only backbone tips.
#' @param species_code One of `"keep"`, `"temporary_species"`, or `"drop"`,
#'   controlling code-like labels such as `Eburopone_MG04`.
#' @param exclude_taxa Optional character vector of exact taxa to protect from
#'   grafting. If present in the backbone, the corresponding branch is renamed to
#'   a protected `TACTEXCL...` pseudo-genus. If absent from the backbone but
#'   present in taxonomy, it is removed from TACT taxonomy.
#' @param exclude_mrca Optional clade exclusions. Supply a character vector of
#'   tip labels or a list of such vectors; each vector is resolved to an MRCA and
#'   all descendant tips are protected.
#' @param exclude_missing One of `"warn"`, `"error"`, or `"ignore"`, controlling
#'   behavior when `exclude_taxa` are not found in the backbone.
#' @param sample_with_replacement Logical. Whether genus-only labels can sample
#'   an already-used species if no unused species remain.
#' @param fallback_prefix Character prefix for temporary placeholder species.
#' @param extra_tact_args Character vector of additional arguments passed to
#'   `tact_add_taxa`.
#' @param keep_temp Logical. Whether to keep the TACT work directory after a
#'   successful run.
#' @param logger Optional logger object.
#'
#' @return Invisibly returns a list of class `tact_grafting_result` containing
#'   the final tree, prepared tree and taxonomy, audit tables, paths, and
#'   parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- run_tact_grafting(
#'   backbone_tree = "project/results/grafted/final_tree.tre",
#'   taxonomy = "project/tables/antwiki-valid-species-8Mar2026.txt",
#'   out_prefix = "project/results/tact/Formicidae_complete_tact",
#'   taxonomy_format = "antwiki",
#'   tact_runner = "docker",
#'   docker_image = "jonchang/tact",
#'   outgroups = NULL,
#'   nonmono = "split",
#'   exclude_taxa = "Neoponera_bucki"
#' )
#' }
run_tact_grafting <- function(backbone_tree,
                              taxonomy,
                              out_prefix,
                              family = "Formicidae",
                              taxonomy_format = c("antwiki", "simple", "tact_csv"),
                              taxonomy_output = c("newick", "csv"),
                              tact_runner = c("docker", "system", "none"),
                              docker_image = "jonchang/tact",
                              tact_bin = "tact_add_taxa",
                              outgroups = NULL,
                              seed = 42,
                              nonmono = c("split", "skip", "error"),
                              nonmono_allocation = c("proportional", "equal", "random"),
                              genus_only = c("replace_random_species", "temporary_species", "keep"),
                              species_code = c("keep", "temporary_species", "drop"),
                              exclude_taxa = NULL,
                              exclude_mrca = NULL,
                              exclude_missing = c("warn", "error", "ignore"),
                              sample_with_replacement = TRUE,
                              fallback_prefix = "sp",
                              extra_tact_args = character(0),
                              keep_temp = TRUE,
                              logger = NULL) {
  .tact_require(c("ape", "stringr", "phytools"))

  taxonomy_format <- match.arg(taxonomy_format)
  taxonomy_output <- match.arg(taxonomy_output)
  tact_runner <- match.arg(tact_runner)
  nonmono <- match.arg(nonmono)
  nonmono_allocation <- match.arg(nonmono_allocation)
  genus_only <- match.arg(genus_only)
  species_code <- match.arg(species_code)
  exclude_missing <- match.arg(exclude_missing)

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  if (is.null(logger)) logger <- .tact_make_logger(out_prefix)

  .tact_section("TACT grafting: input", logger = logger)
  .tact_msg("Backbone tree: ", if (inherits(backbone_tree, "phylo")) "<phylo object>" else backbone_tree, logger = logger)
  .tact_msg("Taxonomy format: ", taxonomy_format, logger = logger)
  .tact_msg("TACT runner: ", tact_runner, logger = logger)

  tree <- if (inherits(backbone_tree, "phylo")) backbone_tree else ape::read.tree(backbone_tree)
  tax <- read_tact_taxonomy(taxonomy, taxonomy_format = taxonomy_format, family = family)

  .tact_msg("Backbone tips: ", ape::Ntip(tree), logger = logger)
  .tact_msg("Taxonomy species rows: ", nrow(tax), logger = logger)

  .tact_section("TACT grafting: label preparation", logger = logger)
  prep <- prepare_tact_backbone_labels(
    tree = tree,
    tax = tax,
    genus_only = genus_only,
    species_code = species_code,
    sample_with_replacement = sample_with_replacement,
    fallback_prefix = fallback_prefix,
    seed = seed,
    logger = logger
  )
  tree <- prep$tree
  tax <- prep$taxonomy

  exclude_tips <- resolve_tact_exclusions(
    tree,
    exclude_taxa = exclude_taxa,
    exclude_mrca = exclude_mrca,
    exclude_missing = exclude_missing
  )
  if (length(exclude_tips)) {
    .tact_msg("TACT-excluded backbone tips: ", length(exclude_tips), logger = logger)
  }

  tax_excl <- .remove_excluded_taxa_not_in_backbone(
    tax = tax,
    exclude_taxa = exclude_taxa,
    backbone_tips = tree$tip.label,
    logger = logger
  )
  tax <- tax_excl$taxonomy
  excluded_taxonomy <- tax_excl$removed

  .tact_section("TACT grafting: non-monophyletic genera", logger = logger)
  spl <- split_nonmono_genera_for_tact(
    tree = tree,
    tax = tax,
    nonmono = nonmono,
    nonmono_allocation = nonmono_allocation,
    exclude_tips = exclude_tips,
    seed = seed,
    logger = logger
  )
  tree_tact <- spl$tree
  tax_tact <- spl$taxonomy

  outgroups_tact <- outgroups
  if (!is.null(outgroups_tact) && length(outgroups_tact)) {
    all_map <- rbind(
      if (nrow(prep$label_map)) {
        data.frame(old_label = prep$label_map$old, new_label = prep$label_map$new)
      } else {
        data.frame(old_label = character(), new_label = character())
      },
      if (nrow(spl$map)) {
        spl$map[, c("old_label", "new_label")]
      } else {
        data.frame(old_label = character(), new_label = character())
      }
    )
    if (nrow(all_map)) {
      idx <- match(outgroups_tact, all_map$old_label)
      outgroups_tact[!is.na(idx)] <- all_map$new_label[idx[!is.na(idx)]]
    }
  }

  work_dir <- paste0(out_prefix, "_tact_work")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  backbone_file <- file.path(work_dir, "backbone_tact_ready.tre")
  taxonomy_csv <- file.path(work_dir, "taxonomy_tact_ready.csv")
  taxonomy_file <- file.path(work_dir, "taxonomy_tact_ready.tre")
  raw_output_prefix <- file.path(work_dir, "tact_output")

  ape::write.tree(tree_tact, backbone_file)
  write_tact_taxonomy_csv(tax_tact, taxonomy_csv, family = family)

  if (tact_runner != "none") {
    run_tact_build_taxonomy_external(
      work_dir = work_dir,
      taxonomy_csv = taxonomy_csv,
      taxonomy_tree = taxonomy_file,
      tact_runner = tact_runner,
      docker_image = docker_image,
      logger = logger
    )
  } else {
    .tact_msg(
      "TACT runner set to 'none'; taxonomy CSV written but taxonomy tree not built by TACT.",
      logger = logger
    )
  }

  utils::write.table(
    prep$label_map,
    paste0(out_prefix, "_label_replacements.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$map,
    paste0(out_prefix, "_nonmono_temp_name_map.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$nonmono,
    paste0(out_prefix, "_nonmono_genera.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    spl$skipped,
    paste0(out_prefix, "_skipped_taxa.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  utils::write.table(
    excluded_taxonomy,
    paste0(out_prefix, "_excluded_taxonomy.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  writeLines(
    if (length(exclude_tips)) exclude_tips else character(0),
    paste0(out_prefix, "_excluded_tips.txt")
  )

  .tact_section("TACT grafting: run external TACT", logger = logger)
  .tact_msg("TACT-ready backbone: ", backbone_file, logger = logger)
  .tact_msg("TACT-ready taxonomy CSV: ", taxonomy_csv, logger = logger)
  .tact_msg("TACT-ready taxonomy tree: ", taxonomy_file, logger = logger)

  if (tact_runner != "none") {
    run_tact_external(
      work_dir = work_dir,
      backbone_file = backbone_file,
      taxonomy_file = taxonomy_file,
      output_prefix = raw_output_prefix,
      outgroups = outgroups_tact,
      tact_runner = tact_runner,
      docker_image = docker_image,
      tact_bin = tact_bin,
      extra_args = extra_tact_args,
      logger = logger
    )
  } else {
    .tact_msg("TACT runner set to 'none'; prepared files only.", logger = logger)
  }

  raw_tree_path <- NA_character_
  final_tree <- NA
  final_path <- NA_character_

  if (tact_runner != "none") {
    raw_tree_path <- .find_tact_output_tree(work_dir, raw_output_prefix)
    if (!is.na(raw_tree_path) && file.exists(raw_tree_path)) {
      raw_tree <- .tact_read_output_tree(raw_tree_path)
      raw_tree <- .tact_normalize_output_labels(raw_tree)
      final_tree <- restore_tact_temp_names(raw_tree, spl$map)
      final_tree <- .tact_drop_invalid_node_labels(final_tree)
      final_path <- paste0(out_prefix, "_tacted_cleaned.tre")
      ape::write.tree(
        phy = final_tree,
        file = final_path
      )
      .tact_msg(
        "Cleaned TACT tip count(s): ",
        paste(.tact_n_tips(final_tree), collapse = ", "),
        logger = logger
      )

      .tact_msg(
        "Remaining temporary TACT labels: ",
        paste(.tact_n_temp_labels(final_tree), collapse = ", "),
        logger = logger
      )
      .tact_msg("Raw TACT output tree detected: ", raw_tree_path, logger = logger)
      .tact_msg("Cleaned TACT tree written: ", final_path, logger = logger)
    } else {
      .tact_msg(
        "TACT completed, but no output tree was detected. Check work directory: ",
        work_dir,
        logger = logger
      )
    }
  } else {
    .tact_msg(
      "TACT runner set to 'none'; no TACT output tree will be read or cleaned.",
      logger = logger
    )
  }

  if (!keep_temp && tact_runner != "none") {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  }

  paths <- list(
    tact_ready_backbone = backbone_file,
    tact_ready_taxonomy_csv = taxonomy_csv,
    tact_ready_taxonomy = taxonomy_file,
    raw_tact_output_tree = raw_tree_path,
    cleaned_tree = final_path,
    label_replacements = paste0(out_prefix, "_label_replacements.tsv"),
    nonmono_map = paste0(out_prefix, "_nonmono_temp_name_map.tsv"),
    nonmono_genera = paste0(out_prefix, "_nonmono_genera.tsv"),
    skipped_taxa = paste0(out_prefix, "_skipped_taxa.tsv"),
    excluded_taxonomy = paste0(out_prefix, "_excluded_taxonomy.tsv"),
    excluded_tips = paste0(out_prefix, "_excluded_tips.txt"),
    log = logger$file,
    work_dir = work_dir
  )

  res <- list(
    tree = final_tree,
    tact_ready_tree = tree_tact,
    tact_ready_taxonomy = tax_tact,
    label_map = prep$label_map,
    nonmono_map = spl$map,
    nonmono_genera = spl$nonmono,
    skipped_taxa = spl$skipped,
    excluded_taxonomy = excluded_taxonomy,
    excluded_tips = exclude_tips,
    paths = paths,
    parameters = list(
      seed = seed,
      nonmono = nonmono,
      nonmono_allocation = nonmono_allocation,
      genus_only = genus_only,
      species_code = species_code,
      tact_runner = tact_runner,
      docker_image = docker_image,
      taxonomy_output = taxonomy_output
    )
  )
  class(res) <- c("tact_grafting_result", class(res))
  invisible(res)
}
