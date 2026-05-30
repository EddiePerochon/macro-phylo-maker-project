#' Extract a clade (genus or MRCA-defined) and optionally append a sister outgroup tip
#'
#' Extracts a focal clade from a larger tree and, if requested, appends \emph{one}
#' outgroup tip chosen from the sister clade by minimal patristic distance to the
#' ingroup. Before extraction, it can (i) resolve polytomies at random,
#' (ii) enforce strictly positive branch lengths, and (iii) canonicalize labels to
#' \code{Genus_species} and collapse intraspecific duplicates. All steps emit
#' messages through AntPhyloMaker's lightweight logger.
#'
#' @details
#' \strong{Workflow}
#' \itemize{
#'   \item \emph{Preprocess (optional)}: resolve polytomies and/or nudge non-positive
#'         edge lengths to small positive values.
#'   \item \emph{Define the ingroup}:
#'         \itemize{
#'           \item If \code{genus} is given, the ingroup is all tips whose
#'                 canonicalized labels start with \code{paste0(genus, "_")}.
#'                 Optionally collapse intraspecific duplicates inside the genus.
#'           \item If \code{mrca_tips} is given, the ingroup is the clade defined by
#'                 the MRCA of those anchors; optionally collapse duplicates by binomial.
#'         }
#'   \item \emph{Non-monophyly (genus mode)}: if genus tips are not monophyletic,
#'         either error out (\code{nonmono="error"}) or prune away non-genus
#'         descendants within the MRCA (\code{nonmono="prune_extras"}, default).
#'   \item \emph{Outgroup (see below)}: add a single outgroup tip if requested.
#'   \item \emph{Write outputs}: the extracted subtree (and optional TSVs with
#'         renames/dropped tips) are written with filenames aligned to the logger
#'         file stem (e.g., \code{<stem>_with_outgroup.tre}).
#' }
#'
#' \strong{Outgroup selection}
#'
#' \describe{
#'   \item{\code{outgroup = "sister_one"}}{Finds the ingroup MRCA's \emph{sister}
#'   clade in the input tree. Among its tips, chooses the label with \emph{minimal}
#'   patristic distance to any ingroup tip (based on the cophenetic matrix).
#'   Ties are broken deterministically by lexicographic order of candidate labels.
#'   The chosen single tip is appended to the output (so the result contains the
#'   ingroup plus that one outgroup tip). If no unique sister exists, the MRCA is
#'   at the root, or no eligible sister tip remains after pruning, \code{"sister_one"}
#'   silently reduces to \code{"none"} and no outgroup is added.}
#'
#'   \item{\code{outgroup = "none"}}{No outgroup is added; the output contains the
#'   ingroup only.}
#' }
#'
#' \strong{Files and naming}
#'
#' If a real smart logger is in use, outputs are named from the logger's file stem,
#' for example:
#' \itemize{
#'   \item \code{<stem>_with_outgroup.tre} or \code{<stem>_no_outgroup.tre}
#'   \item \code{<stem>_renamed.tsv} (canonicalization map, when available)
#'   \item \code{<stem>_dropped.tsv} (intraspecific duplicates removed)
#' }
#' When \code{logger = NULL}, \code{TRUE}, or a directory path (see \emph{logger}
#' below), a logger is created for this call and its stem is used for all sidecars.
#'
#' @param tree A \code{phylo} object.
#' @param genus Optional character scalar. If supplied, the ingroup is defined as
#'   tips whose canonicalized label begins with \code{paste0(genus, "_")}.
#' @param mrca_tips Optional character vector of tip labels whose MRCA defines
#'   the clade. Exactly one of \code{genus} or \code{mrca_tips} must be supplied.
#' @param outgroup One of \code{"sister_one"} (default; see above) or \code{"none"}.
#' @param clean One of \code{"genus_species"} (default in genus mode) or \code{"none"}.
#'   With \code{genus}, canonicalize and collapse intraspecific duplicates; with
#'   \code{mrca_tips}, collapse duplicates by binomial if \code{"genus_species"}.
#' @param nonmono For non-monophyly in \code{genus} mode: \code{"prune_extras"} (default)
#'   keeps only genus tips within the MRCA; \code{"error"} aborts.
#' @param resolve_polytomies Logical; if \code{TRUE} (default) resolves polytomies randomly.
#' @param force_positive_lengths Logical; if \code{TRUE} (default) enforces strictly positive edges.
#' @param seed Integer seed for randomized steps (polytomy resolution; rare tie situations).
#' @param write_tree Logical; if \code{TRUE} (default) write the extracted subtree \code{.tre}.
#' @param tree_path Optional explicit output path for the Newick; auto-named from the
#'   logger stem if \code{NULL}.
#' @param write_renames Logical; if \code{TRUE} write TSV of raw→canonical renames (ingroup).
#' @param renames_path Optional TSV path for renames; auto-named if \code{NULL}.
#' @param write_drops Logical; if \code{TRUE} write TSV of dropped duplicate tips (ingroup).
#' @param drops_path Optional TSV path for dropped duplicates; auto-named if \code{NULL}.
#' @param out_dir Base directory for auto-named files (default \code{"."}).
#' @param logger Logging control.  Optional AntPhyloMaker smart logger (from \code{".make_logger()"}").
#'   Either NULL (create a per-run logger in `out_dir`) or a smart logger object created by \code{".make_logger()"}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{tree}: extracted \code{phylo} (with/without outgroup),
#'   \item \code{outgroup}: character (chosen outgroup) or \code{NA},
#'   \item \code{ingroup_tips}: character vector of ingroup tips in the output,
#'   \item \code{paths}: named list of written file paths (\code{tree}, \code{renames}, \code{drops}, when present).
#' }
#'
#'
#' @examples
#' \dontrun{
#' ## Example 1: genus mode with outgroup
#' tr <- ape::read.tree(text = "((Camponotus_a:1,Camponotus_b:1):1,(Formica_c:1,Lasius_d:1):1);")
#' res <- extract_clade_with_outgroup(
#'   tr,
#'   genus = "Camponotus", outgroup = "sister_one",
#'   out_dir = tempdir()
#' )
#' res$paths$tree # "<stem>_with_outgroup.tre"
#'
#' ## Example 2: MRCA-defined clade without outgroup, custom output path
#' nwk <- "((A:1,B:1):1,(C:1,D:1):1);"
#' tr2 <- ape::read.tree(text = nwk)
#' res2 <- extract_clade_with_outgroup(
#'   tr2,
#'   mrca_tips = c("A", "B"), outgroup = "none",
#'   write_tree = TRUE, tree_path = file.path(tempdir(), "AB_subtree.tre")
#' )
#' }
#'
#' @seealso \code{prepare_clade_template()}, \code{graft_template_onto_backbone()},
#'   \code{graft_many_clades()}, \code{.make_logger()}, \code{log_msg()}, \code{log_section()}.
#' @md
#' @export

extract_clade_with_outgroup <- function(
  tree,
  genus = NULL,
  mrca_tips = NULL,
  outgroup = c("sister_one", "none"),
  clean = c("genus_species", "none"),
  nonmono = c("prune_extras", "error"),
  resolve_polytomies = TRUE,
  force_positive_lengths = TRUE,
  seed = 42L,
  write_tree = TRUE,
  tree_path = NULL,
  write_renames = TRUE,
  renames_path = NULL,
  write_drops = TRUE,
  drops_path = NULL,
  out_dir = ".",
  logger = NULL
) {
  stopifnot(inherits(tree, "phylo"))
  outgroup <- match.arg(outgroup)
  nonmono <- match.arg(nonmono)

  if (is.null(genus) && is.null(mrca_tips)) stop("Provide either 'genus' or 'mrca_tips'.")
  if (!is.null(genus) && !is.null(mrca_tips)) stop("Provide only one of 'genus' or 'mrca_tips'.")

  clean <- match.arg(clean)
  if (is.null(genus) && identical(clean, "genus_species")) clean <- "none"

  # ---- logger setup ----------------------------------------------------------
  # Build a stable prefix for file names & logs
  prefix <- if (!is.null(genus)) {
    paste0("extract_", genus)
  } else {
    anchors <- paste0(utils::head(mrca_tips, 3), collapse = "-")
    if (length(mrca_tips) > 3) anchors <- paste0(anchors, "-etc")
    paste0("extract_mrca_", anchors)
  }

  created_logger <- FALSE
  if (is.null(logger)) {
    # Create a per-call logger (append console + file)
    logger <- .make_logger(
      outdir = out_dir,
      genus = if (!is.null(genus)) genus else "clade",
      mode = "extract",
      file_prefix = prefix
    )
    created_logger <- TRUE
    on.exit(.close_logger(logger), add = TRUE)
  }

  # Helper: base stem for sidecar files (.tre / .tsv) aligned with the log file
  log_file <- if (inherits(logger, "smart_logger")) logger$file else file.path(out_dir, paste0(prefix, ".log"))
  stem <- sub("\\.log$", "", log_file)

  # ---- sections --------------------------------------------------------------
  log_section(logger, "Settings")
  log_msg(logger, "resolve_polytomies: ", resolve_polytomies)
  log_msg(logger, "force_positive_lengths: ", force_positive_lengths)
  log_msg(logger, "clean: ", clean)
  log_msg(logger, "nonmono: ", nonmono)
  log_msg(logger, "outgroup: ", outgroup)
  log_msg(logger, "seed: ", as.integer(seed))

  tr0 <- tree
  log_section(logger, "Input")
  log_msg(logger, "tips_in: ", length(tr0$tip.label))

  # ---- 0) resolve polytomies / positive lengths -----------------------------
  tr <- tr0
  if (isTRUE(resolve_polytomies) || isTRUE(force_positive_lengths)) {
    tr <- .resolve_polytomies_positive_lengths(tr,
      seed = seed,
      do_resolve = resolve_polytomies,
      do_positive = force_positive_lengths
    )
    log_msg(
      logger, "Preprocess: polytomies=", resolve_polytomies,
      ", positive_lengths=", force_positive_lengths
    )
  }

  # ---- 1) define ingroup & optionally clean/collapse -------------------------
  tips_all <- tr$tip.label
  renames_df <- data.frame(raw = character(0), clean = character(0), stringsAsFactors = FALSE)
  dropped_vec <- character(0)

  if (!is.null(genus)) {
    if (identical(clean, "genus_species")) {
      ci <- .clean_ingroup_labels(tr, genus = genus) # list(tree, renamed, dropped)
      tr <- ci$tree
      if (nrow(ci$renamed)) renames_df <- rbind(renames_df, ci$renamed)
      if (length(ci$dropped)) dropped_vec <- unique(c(dropped_vec, ci$dropped))

      cd <- .collapse_species_duplicates_ingroup(tr, genus = genus) # list(tree, dropped)
      tr <- cd$tree
      if (length(cd$dropped)) dropped_vec <- unique(c(dropped_vec, cd$dropped))
      tips_all <- tr$tip.label
    }
    ingroup_tips <- tips_all[vapply(tips_all, function(x) {
      parts <- strsplit(x, "_", fixed = TRUE)[[1]]
      length(parts) >= 2 && parts[1] == genus
    }, logical(1))]
    if (!length(ingroup_tips)) stop("No tips found for genus '", genus, "' after cleaning.")
    log_msg(logger, "Ingroup (genus=", genus, ") identified: ", length(ingroup_tips), " tips")

    is_mono <- ape::is.monophyletic(tr, ingroup_tips)
    log_msg(logger, "Genus monophyly: ", is_mono)
    if (!is_mono && identical(nonmono, "error")) {
      stop("Genus '", genus, "' is not monophyletic. Use nonmono='prune_extras'.")
    }
    mrca_node <- ape::getMRCA(tr, ingroup_tips)
    if (is.na(mrca_node)) stop("Could not identify MRCA for genus '", genus, "'.")
    sub <- ape::extract.clade(tr, mrca_node)
    if (!is_mono && identical(nonmono, "prune_extras")) {
      keep <- intersect(sub$tip.label, ingroup_tips)
      sub <- if (length(keep)) ape::keep.tip(sub, keep) else sub
      log_msg(logger, "Pruned non-genus tips inside MRCA (prune_extras).")
    }
    ingroup_tips <- sub$tip.label
  } else {
    if (!all(mrca_tips %in% tips_all)) {
      miss <- setdiff(mrca_tips, tips_all)
      stop("mrca_tips not in tree: ", paste(miss, collapse = ", "))
    }
    mrca_node <- ape::getMRCA(tr, mrca_tips)
    if (is.na(mrca_node)) stop("Could not compute MRCA for supplied anchors.")
    sub <- ape::extract.clade(tr, mrca_node)
    ingroup_tips <- sub$tip.label
    log_msg(logger, "Ingroup (MRCA of anchors) size: ", length(ingroup_tips))
    if (identical(clean, "genus_species")) {
      cd2 <- .collapse_species_duplicates_by_binomial(sub) # list(tree, dropped)
      sub <- cd2$tree
      if (length(cd2$dropped)) dropped_vec <- unique(c(dropped_vec, cd2$dropped))
      ingroup_tips <- sub$tip.label
      log_msg(logger, "Collapsed intraspecific duplicates by binomial (MRCA mode).")
    }
  }

  # ---- 2) optional sister outgroup ------------------------------------------
  chosen_out <- NA_character_
  min_d <- NA_real_
  if (!identical(outgroup, "none")) {
    parent <- .parent_of(tr)
    parent_node <- parent[mrca_node]
    if (!is.na(parent_node) && parent_node > 0) {
      kids <- tr$edge[tr$edge[, 1] == parent_node, 2]
      sis <- setdiff(kids, mrca_node)
      if (length(sis) == 1) {
        sis_desc <- phytools::getDescendants(tr, sis)
        sis_tips <- sis_desc[sis_desc <= length(tr$tip.label)]
        sis_labels <- setdiff(tr$tip.label[sis_tips], ingroup_tips)
        if (length(sis_labels) >= 1) {
          D <- ape::cophenetic.phylo(tr)
          mins <- sapply(sis_labels, function(x) min(D[x, ingroup_tips], na.rm = TRUE))
          min_d <- min(mins, na.rm = TRUE)
          candidates <- sort(names(which(mins == min_d)))
          set.seed(as.integer(seed))
          chosen_out <- candidates[1]
          log_msg(logger, "Outgroup (sister_one): ", chosen_out, " (min distance=", sprintf("%.6f", min_d), ")")
        } else {
          log_msg(logger, "No eligible sister tips for outgroup.")
        }
      } else {
        log_msg(logger, "No unique sister clade for MRCA; outgroup skipped.")
      }
    } else {
      log_msg(logger, "MRCA at root or no parent; outgroup skipped.")
    }
  }

  out_tree <- if (is.na(chosen_out)) sub else ape::keep.tip(tr, c(ingroup_tips, chosen_out))
  log_section(logger, "Output")
  log_msg(logger, "ingroup_out: ", length(ingroup_tips))
  log_msg(logger, "tips_out (final tree): ", length(out_tree$tip.label))

  # ---- 3) write files (tree + renames/drops) --------------------------------
  paths <- list()
  # Tree path
  if (isTRUE(write_tree)) {
    if (is.null(tree_path)) {
      suffix <- if (is.na(chosen_out)) "_no_outgroup" else "_with_outgroup"
      tree_path <- paste0(stem, suffix, ".tre")
    }
    dir.create(dirname(normalizePath(tree_path, mustWork = FALSE)), recursive = TRUE, showWarnings = FALSE)
    ape::write.tree(out_tree, file = tree_path)
    log_msg(logger, "Wrote tree: ", tree_path)
    paths$tree <- tree_path
  }

  # Renames TSV (ingroup)
  if (isTRUE(write_renames) && nrow(renames_df)) {
    if (is.null(renames_path)) renames_path <- paste0(stem, "_renamed.tsv")
    utils::write.table(renames_df, file = renames_path, sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg(logger, "Wrote renames TSV: ", renames_path)
    paths$renames <- renames_path
  }

  # Drops TSV (ingroup duplicates)
  if (isTRUE(write_drops) && length(dropped_vec)) {
    if (is.null(drops_path)) drops_path <- paste0(stem, "_dropped.tsv")
    utils::write.table(data.frame(dropped = sort(unique(dropped_vec))),
      file = drops_path, sep = "\t", row.names = FALSE, quote = FALSE
    )
    log_msg(logger, "Wrote drops TSV: ", drops_path)
    paths$drops <- drops_path
  }

  invisible(list(
    tree = out_tree,
    outgroup = chosen_out,
    ingroup_tips = ingroup_tips,
    paths = paths
  ))
}

# ======================= internal helpers =======================

.resolve_polytomies_positive_lengths <- function(tr, seed = 42L,
                                                 do_resolve = TRUE,
                                                 do_positive = TRUE) {
  stopifnot(inherits(tr, "phylo"))
  out <- tr
  if (isTRUE(do_resolve)) {
    set.seed(as.integer(seed))
    out <- ape::multi2di(out, random = TRUE)
  }
  if (isTRUE(do_positive)) {
    if (is.null(out$edge.length)) {
      out$edge.length <- rep(1, nrow(out$edge))
    } else {
      el <- out$edge.length
      pos <- el[is.finite(el) & el > 0]
      eps <- if (length(pos)) min(pos) * 1e-6 else 1e-6
      el[!is.finite(el) | el <= 0] <- eps
      out$edge.length <- el
    }
  }
  out
}

.parent_of <- function(tr) {
  n <- max(tr$edge)
  parent <- rep(NA_integer_, n)
  parent[tr$edge[, 2]] <- tr$edge[, 1]
  parent
}

# list(tree, renamed=data.frame, dropped=character)
.clean_ingroup_labels <- function(tree, genus) {
  labs <- tree$tip.label
  idx <- which(vapply(labs, function(lbl) {
    bin <- .clean_label_to_binomial(lbl, genus_hint = NULL)
    if (is.na(bin)) {
      return(FALSE)
    }
    strsplit(bin, "_", fixed = TRUE)[[1]][1] == genus
  }, logical(1)))
  if (!length(idx)) {
    return(list(
      tree = tree,
      renamed = data.frame(raw = character(0), clean = character(0), stringsAsFactors = FALSE),
      dropped = character(0)
    ))
  }
  cleaned <- vapply(labs[idx], .clean_label_to_binomial, character(1), genus_hint = genus)
  bad <- which(is.na(cleaned) | !nzchar(cleaned))
  keepers <- setdiff(seq_along(idx), bad)
  ren_map <- data.frame(raw = character(0), clean = character(0), stringsAsFactors = FALSE)
  if (length(keepers)) {
    ren_map <- data.frame(raw = labs[idx[keepers]], clean = cleaned[keepers], stringsAsFactors = FALSE)
    labs[idx[keepers]] <- cleaned[keepers]
  }
  dropped <- character(0)
  if (length(bad)) {
    to_drop <- labs[idx[bad]]
    tree <- ape::drop.tip(tree, to_drop)
    dropped <- to_drop
    labs <- tree$tip.label
  }
  tree$tip.label <- labs
  list(tree = tree, renamed = ren_map, dropped = dropped)
}

# list(tree, dropped=character)
.collapse_species_duplicates_ingroup <- function(tr, genus) {
  labs <- tr$tip.label
  idx <- which(vapply(labs, function(x) {
    parts <- strsplit(x, "_", fixed = TRUE)[[1]]
    length(parts) >= 2 && parts[1] == genus
  }, logical(1)))
  if (!length(idx)) {
    return(list(tree = tr, dropped = character(0)))
  }
  ing <- labs[idx]
  keys <- sub("^([A-Za-z]+_[a-z]+).*$", "\\1", ing, perl = TRUE)
  tab <- table(keys)
  if (!any(tab > 1)) {
    return(list(tree = tr, dropped = character(0)))
  }
  D <- ape::cophenetic.phylo(tr)
  to_drop <- character(0)
  for (sp in names(tab[tab > 1])) {
    sp_tips <- ing[keys == sp]
    if (length(sp_tips) < 2) next
    if (ape::is.monophyletic(tr, sp_tips)) {
      mnode <- ape::getMRCA(tr, sp_tips)
      H <- phytools::nodeHeights(tr)
      root_age <- max(H[, 2])
      nh_mrca <- phytools::nodeheight(tr, mnode)
      tip_nodes <- match(sp_tips, labs)
      tip_heights <- setNames(sapply(tip_nodes, function(n) phytools::nodeheight(tr, n)), sp_tips)
      crown_depths <- (root_age - tip_heights) - (root_age - nh_mrca)
      keep <- names(which.min(crown_depths))
      to_drop <- c(to_drop, setdiff(sp_tips, keep))
    } else {
      sums <- sapply(sp_tips, function(x) sum(D[x, setdiff(sp_tips, x)], na.rm = TRUE))
      keep <- names(which.min(sums))
      to_drop <- c(to_drop, setdiff(sp_tips, keep))
    }
  }
  if (length(to_drop)) tr <- ape::drop.tip(tr, unique(to_drop))
  list(tree = tr, dropped = unique(to_drop))
}

# list(tree, dropped=character)
.collapse_species_duplicates_by_binomial <- function(tr) {
  labs <- tr$tip.label
  has_bin <- grepl("^[A-Z][a-z]+_[a-z]+", labs)
  if (!any(has_bin)) {
    return(list(tree = tr, dropped = character(0)))
  }
  keys <- sub("^([A-Za-z]+_[a-z]+).*$", "\\1", labs[has_bin], perl = TRUE)
  tab <- table(keys)
  if (!any(tab > 1)) {
    return(list(tree = tr, dropped = character(0)))
  }
  D <- ape::cophenetic.phylo(tr)
  to_drop <- character(0)
  for (sp in names(tab[tab > 1])) {
    sp_tips <- labs[has_bin][keys == sp]
    if (length(sp_tips) < 2) next
    sums <- sapply(sp_tips, function(x) sum(D[x, setdiff(sp_tips, x)], na.rm = TRUE))
    keep <- names(which.min(sums))
    to_drop <- c(to_drop, setdiff(sp_tips, keep))
  }
  if (length(to_drop)) tr <- ape::drop.tip(tr, unique(to_drop))
  list(tree = tr, dropped = unique(to_drop))
}

.clean_label_to_binomial <- function(label, genus_hint = NULL) {
  if (is.na(label) || !nzchar(label)) {
    return(NA_character_)
  }
  x <- label
  x <- gsub("\\s+", "_", x)
  x <- gsub("_(nr|cf|aff)_", "_", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("_([a-z]+)_(nr|cf|aff)", "_\\1_", x, perl = TRUE)
  x <- sub("(_EX\\d+.*)$", "", x, perl = TRUE)
  x <- sub("(_CASENT\\d+.*)$", "", x, perl = TRUE)
  x <- sub("(_MAMI\\d+.*)$", "", x, perl = TRUE)
  x <- sub("(_D\\d+.*)$", "", x, perl = TRUE)
  m <- regexpr("[A-Z][a-z]+_[a-z]+", x, perl = TRUE)
  if (m[1] == -1) {
    return(NA_character_)
  }
  bin <- substring(x, m[1], m[1] + attr(m, "match.length") - 1)
  parts <- strsplit(bin, "_", fixed = TRUE)[[1]]
  Genus <- parts[1]
  species <- tolower(parts[2])
  if (!is.null(genus_hint) && nzchar(genus_hint) && Genus != genus_hint) Genus <- genus_hint
  paste(Genus, species, sep = "_")
}
