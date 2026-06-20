library(testthat)
library(ape)
library(phytools)

# -------------------------------------------------------------------
# Namespace helpers
# -------------------------------------------------------------------

pkg_ns <- asNamespace("MacroPhyloMaker")
ns_fun <- function(name) get(name, envir = pkg_ns, inherits = FALSE)

normalize_stem_mode <- ns_fun(".normalize_stem_mode")
infer_single_outgroup <- ns_fun(".infer_single_outgroup")
fallback_single_outgroup <- ns_fun(".fallback_single_outgroup")
compute_crown_metrics <- ns_fun(".compute_crown_metrics")
rescale_to_height <- ns_fun(".rescale_to_height")
safe_force_ultrametric <- ns_fun(".safe_force_ultrametric")
collapse_mrca_to_existing_tip <- ns_fun(".collapse_mrca_to_existing_tip")
ensure_positive_edges <- ns_fun(".ensure_positive_edges")

# Optional logger hooks (guarded below)
has_logger_utils <- all(vapply(
  c(".make_logger", ".close_logger", "log_msg", "log_section"),
  exists, logical(1), envir = pkg_ns, inherits = FALSE
))

# -------------------------------------------------------------------
# Tree fixtures
# -------------------------------------------------------------------

make_backbone <- function() {
  # Ultrametric 4-tip backbone
  read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
}

make_backbone_missing_A <- function() {
  read.tree(text = "((X:1,Y:1):1,(C:1,D:1):1);")
}

make_donor_with_outgroup <- function() {
  # Root singleton outgroup; ultrametric donor
  read.tree(text = "((Genus_alpha:1,Genus_beta:1):1,Outgroup:2);")
}

make_donor_single_species_samples <- function() {
  # Intended for optional regression test only
  read.tree(text = "((Genus_alpha_1:1,Genus_alpha_2:1):1,Outgroup:2);")
}

make_template <- function() {
  list(
    ok = TRUE,
    donor_unit = read.tree(text = "(Genus_alpha:1,Genus_beta:1);"),
    stem_fraction = 0.25,
    stem_len = 0.5,
    label_map = NULL,
    template_log = data.frame(
      n_ingroup_kept = 2L,
      outgroup_used = "Outgroup",
      stem_len_inferred = 0.5,
      stem_fraction = 0.25,
      crown_fraction = 0.75,
      stringsAsFactors = FALSE
    ),
    logfile = NA_character_
  )
}

# -------------------------------------------------------------------
# Topology helpers
# -------------------------------------------------------------------

tip_desc_labels <- function(tree, node) {
  desc <- phytools::getDescendants(tree, node)
  tip_idx <- desc[desc <= ape::Ntip(tree)]
  tree$tip.label[tip_idx]
}

expect_tip_set <- function(tree, expected) {
  expect_setequal(unname(tree$tip.label), unname(expected))
}

expect_sister_pair <- function(tree, tip1, tip2) {
  mrca <- ape::getMRCA(tree, c(tip1, tip2))
  expect_false(is.na(mrca))
  labels <- tip_desc_labels(tree, mrca)
  expect_equal(sort(labels), sort(c(tip1, tip2)))
}

expect_tip_absent <- function(tree, tip) {
  expect_false(tip %in% tree$tip.label)
}

expect_edge_lengths_positive <- function(tree) {
  expect_true(!is.null(tree$edge.length))
  expect_true(all(is.finite(tree$edge.length)))
  expect_true(all(tree$edge.length > 0))
}

# -------------------------------------------------------------------
# resolve_authority_binomials()
# -------------------------------------------------------------------

test_that("resolve_authority_binomials handles NULL, vectors, data frames, and files", {
  expect_null(resolve_authority_binomials(NULL))

  expect_equal(
    sort(resolve_authority_binomials(c("Genus_alpha", "Genus_alpha", "Genus_beta"))),
    c("Genus_alpha", "Genus_beta")
  )

  df_gs <- data.frame(
    Genus = c("Genus", "Genus"),
    Species = c("alpha", "beta"),
    stringsAsFactors = FALSE
  )
  expect_equal(
    sort(resolve_authority_binomials(df_gs)),
    c("Genus_alpha", "Genus_beta")
  )

  df_bin <- data.frame(
    Binomial = c("Genus_alpha", "Genus_beta"),
    stringsAsFactors = FALSE
  )
  expect_equal(
    sort(resolve_authority_binomials(df_bin)),
    c("Genus_alpha", "Genus_beta")
  )

  tf <- tempfile(fileext = ".tsv")
  write.table(df_bin, tf, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_equal(
    sort(resolve_authority_binomials(tf)),
    c("Genus_alpha", "Genus_beta")
  )
})

test_that("resolve_authority_binomials rejects malformed data frames", {
  bad <- data.frame(foo = "bar", stringsAsFactors = FALSE)

  expect_error(
    resolve_authority_binomials(bad),
    "Authority data.frame"
  )
})

# -------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------

test_that(".normalize_stem_mode normalizes aliases and defaults", {
  expect_equal(normalize_stem_mode("outgroup"), "outgroup")
  expect_equal(normalize_stem_mode("stem"), "outgroup")
  expect_equal(normalize_stem_mode("with_outgroup"), "outgroup")

  expect_equal(normalize_stem_mode("crown"), "crown")
  expect_equal(normalize_stem_mode("crown_only"), "crown")
  expect_equal(normalize_stem_mode("ingroup_only"), "crown")

  expect_equal(normalize_stem_mode(NA_character_), "outgroup")
  expect_equal(normalize_stem_mode(""), "outgroup")
})

test_that("outgroup inference detects singleton root outgroup and fallback returns a tip", {
  donor <- make_donor_with_outgroup()

  expect_equal(infer_single_outgroup(donor), "Outgroup")

  fb <- fallback_single_outgroup(donor)
  expect_true(is.character(fb))
  expect_equal(length(fb), 1L)
  expect_true(fb %in% donor$tip.label)
})

test_that(".compute_crown_metrics returns finite stem/crown metrics", {
  tr <- make_donor_with_outgroup()
  ig <- c("Genus_alpha", "Genus_beta")

  met <- compute_crown_metrics(tr, ig)

  expect_true(isTRUE(met$ok))
  expect_true(is.finite(met$stem_len))
  expect_true(is.finite(met$crown_depth))
  expect_true(is.finite(met$total))
  expect_gt(met$r, 0)
  expect_lt(met$r, 1)
})

test_that(".rescale_to_height and .ensure_positive_edges behave as expected", {
  tr <- read.tree(text = "(A:2,B:2);")
  scaled <- rescale_to_height(tr, 1)

  h <- max(phytools::nodeHeights(scaled))
  expect_equal(h, 1, tolerance = 1e-8)

  bad <- tr
  bad$edge.length[1] <- 0
  fixed <- ensure_positive_edges(bad)
  expect_edge_lengths_positive(fixed)
})

test_that(".collapse_mrca_to_existing_tip keeps exactly one deterministic descendant tip", {
  bb <- make_backbone()
  mrca <- ape::getMRCA(bb, c("A", "B"))

  out <- collapse_mrca_to_existing_tip(bb, mrca)

  expect_s3_class(out$tree, "phylo")
  expect_equal(length(out$tree$tip.label), 3L)
  expect_equal(out$tip, "A")  # deterministic sort() choice
  expect_true("A" %in% out$tree$tip.label)
  expect_false("B" %in% out$tree$tip.label)
})

# -------------------------------------------------------------------
# prepare_clade_template()
# -------------------------------------------------------------------

test_that("prepare_clade_template infers singleton outgroup and returns unit donor crown", {
  donor <- make_donor_with_outgroup()

  tpl <- prepare_clade_template(
    donor_tree = donor,
    ingroup_tips = c("Genus_alpha", "Genus_beta"),
    chronos_select = "off",
    stem_mode = "outgroup",
    resolve_polytomies = "none",
    logger = NULL,
    expand_ingroup_to_full_donor = TRUE
  )

  expect_true(isTRUE(tpl$ok))
  expect_s3_class(tpl$donor_unit, "phylo")
  expect_tip_set(tpl$donor_unit, c("Genus_alpha", "Genus_beta"))

  h <- max(phytools::nodeHeights(tpl$donor_unit))
  expect_equal(h, 1, tolerance = 1e-8)

  expect_true(is.finite(tpl$stem_fraction))
  expect_gt(tpl$stem_fraction, 0)
  expect_lt(tpl$stem_fraction, 1)

  expect_equal(tpl$template_log$n_ingroup_kept, 2L)
  expect_equal(as.character(tpl$template_log$outgroup_used), "Outgroup")
})

test_that("prepare_clade_template respects crown mode and authority filtering", {
  donor <- make_donor_with_outgroup()

  tpl <- prepare_clade_template(
    donor_tree = donor,
    ingroup_tips = c("Genus_alpha", "Genus_beta"),
    authority = c("Genus_alpha", "Genus_beta"),
    chronos_select = "off",
    stem_mode = "crown",
    resolve_polytomies = "none",
    logger = NULL,
    expand_ingroup_to_full_donor = TRUE
  )

  expect_true(isTRUE(tpl$ok))
  expect_s3_class(tpl$donor_unit, "phylo")
  expect_tip_set(tpl$donor_unit, c("Genus_alpha", "Genus_beta"))
  expect_true(is.na(tpl$stem_fraction))
})

# -------------------------------------------------------------------
# graft_template_onto_backbone()
# -------------------------------------------------------------------

test_that("tip graft replaces the target tip with the donor crown", {
  bb <- make_backbone()
  tpl <- make_template()

  res <- graft_template_onto_backbone(
    backbone_ultra = bb,
    placement = list(type = "tip", tip = "A", stem_mode = "outgroup"),
    template = tpl
  )

  expect_s3_class(res$tree, "phylo")
  expect_tip_set(res$tree, c("B", "C", "D", "Genus_alpha", "Genus_beta"))
  expect_tip_absent(res$tree, "A")
  expect_sister_pair(res$tree, "Genus_alpha", "Genus_beta")
  expect_edge_lengths_positive(res$tree)

  expect_true(is.data.frame(res$log))
  expect_equal(as.character(res$log$target_type), "tip")
  expect_equal(as.character(res$log$target_label), "A")
})

test_that("clade graft replaces the target MRCA clade with the donor crown (outgroup mode)", {
  bb <- make_backbone()
  tpl <- make_template()

  res <- graft_template_onto_backbone(
    backbone_ultra = bb,
    placement = list(type = "clade", anchors = c("A", "B"), stem_mode = "outgroup"),
    template = tpl
  )

  expect_s3_class(res$tree, "phylo")
  expect_tip_set(res$tree, c("C", "D", "Genus_alpha", "Genus_beta"))
  expect_false(any(c("A", "B") %in% res$tree$tip.label))
  expect_sister_pair(res$tree, "Genus_alpha", "Genus_beta")
  expect_edge_lengths_positive(res$tree)

  expect_true(is.data.frame(res$log))
  expect_equal(as.character(res$log$target_type), "clade")
  expect_equal(as.character(res$log$target_label), "A,B")
})

test_that("clade graft in crown mode also replaces the target clade", {
  bb <- make_backbone()
  tpl <- make_template()

  res <- graft_template_onto_backbone(
    backbone_ultra = bb,
    placement = list(type = "clade", anchors = c("A", "B"), stem_mode = "crown"),
    template = tpl
  )

  expect_s3_class(res$tree, "phylo")
  expect_tip_set(res$tree, c("C", "D", "Genus_alpha", "Genus_beta"))
  expect_false(any(c("A", "B") %in% res$tree$tip.label))
  expect_sister_pair(res$tree, "Genus_alpha", "Genus_beta")
  expect_edge_lengths_positive(res$tree)
})

test_that("graft_template_onto_backbone validates anchors and placement type", {
  bb <- make_backbone()
  tpl <- make_template()

  expect_error(
    graft_template_onto_backbone(
      bb,
      placement = list(type = "clade", anchors = c("A", "Z"), stem_mode = "outgroup"),
      template = tpl
    ),
    "Anchors not found"
  )

  expect_error(
    graft_template_onto_backbone(
      bb,
      placement = list(type = "weird"),
      template = tpl
    ),
    "placement\\$type"
  )
})

# -------------------------------------------------------------------
# graft_many_clades()
# -------------------------------------------------------------------

test_that("graft_many_clades processes a pre-built template and writes outputs", {
  bb <- make_backbone()
  tpl <- make_template()

  out_tree <- tempfile(fileext = ".tre")
  out_log  <- tempfile(fileext = ".tsv")

  res <- graft_many_clades(
    backbones = bb,
    donors = list(d1 = tpl),
    placements = list(d1 = list(type = "tip", tip = "A", stem_mode = "outgroup")),
    ultrametric_final = "none",
    out_tree = out_tree,
    out_log = out_log,
    logger = NULL
  )

  expect_s3_class(res$tree, "phylo")
  expect_tip_set(res$tree, c("B", "C", "D", "Genus_alpha", "Genus_beta"))
  expect_true(file.exists(out_tree))
  expect_true(file.exists(out_log))
  expect_true(nrow(res$log) == 1L)
})

test_that("graft_many_clades skips backbones missing the target and preserves others", {
  bb1 <- make_backbone()
  bb2 <- make_backbone_missing_A()
  backbones <- list(bb1, bb2)
  class(backbones) <- "multiPhylo"
  attr(backbones, "order") <- NULL

  tpl <- make_template()

  out_tree <- tempfile(fileext = ".tre")
  out_log  <- tempfile(fileext = ".tsv")

  res <- graft_many_clades(
    backbones = backbones,
    donors = list(d1 = tpl),
    placements = list(d1 = list(type = "tip", tip = "A", stem_mode = "outgroup")),
    ultrametric_final = "none",
    out_tree = out_tree,
    out_log = out_log,
    logger = NULL
  )

  expect_s3_class(res$tree, "multiPhylo")

  # First backbone gets grafted
  expect_true(all(c("Genus_alpha", "Genus_beta") %in% res$tree[[1]]$tip.label))

  # Second backbone lacks A, so it should be skipped as coded
  expect_false(any(c("Genus_alpha", "Genus_beta") %in% res$tree[[2]]$tip.label))

  # Only one successful graft row logged
  expect_true(nrow(res$log) == 1L)
})

# -------------------------------------------------------------------
# run_clade_grafting() end-to-end
# -------------------------------------------------------------------

test_that("run_clade_grafting completes an end-to-end MRCA graft on a simple plan", {
  # This wrapper depends on utilities defined across the package.
  # Guard rather than fail for missing unrelated package modules.
  needed <- c("read_trees_any", "set_global_seed", "safe_drop_tips",
              "extract_ingroup_by_anchors")
  skip_if(any(!vapply(needed, exists, logical(1), envir = pkg_ns, inherits = FALSE)),
          "run_clade_grafting dependencies are defined outside this module")

  bb <- make_backbone()
  donor <- make_donor_with_outgroup()

  bb_path <- tempfile(fileext = ".tre")
  donor_path <- tempfile(fileext = ".tre")
  ape::write.tree(bb, bb_path)
  ape::write.tree(donor, donor_path)

  plan <- data.frame(
    MRCA = "A,B",
    Phylogeny_file_path = donor_path,
    Input_tree_type = "chronogram",
    Stem_mode = "outgroup",
    stringsAsFactors = FALSE
  )
  plan_path <- tempfile(fileext = ".tsv")
  write.table(plan, plan_path, sep = "\t", row.names = FALSE, quote = FALSE)

  out_prefix <- tempfile("clade_graft_")

  res <- run_clade_grafting(
    backbone_path = bb_path,
    plan_path = plan_path,
    out_prefix = out_prefix,
    seed_mode = 42,
    ensure_ultrametric = "none",
    chronos_select = "off",
    ultrametric_final = "none",
    plot_pdf = FALSE,
    logger = NULL
  )

  expect_true(file.exists(paste0(out_prefix, ".tre")))
  expect_true(file.exists(paste0(out_prefix, "_graft_log.tsv")))

  out_tree <- ape::read.tree(paste0(out_prefix, ".tre"))
  expect_tip_set(out_tree, c("C", "D", "Genus_alpha", "Genus_beta"))
})

# -------------------------------------------------------------------
# Optional logger integration test
# -------------------------------------------------------------------

test_that("prepare_clade_template writes expected logger sections when logger utilities exist", {
  skip_if_not(has_logger_utils)

  make_logger  <- ns_fun(".make_logger")
  close_logger <- ns_fun(".close_logger")

  dir <- tempdir()
  logger <- make_logger(
    outdir = dir,
    genus = "Clade",
    mode = "graft",
    file_prefix = "template_test"
  )
  on.exit(try(close_logger(logger), silent = TRUE), add = TRUE)

  donor <- make_donor_with_outgroup()

  tpl <- prepare_clade_template(
    donor_tree = donor,
    ingroup_tips = c("Genus_alpha", "Genus_beta"),
    chronos_select = "off",
    stem_mode = "outgroup",
    resolve_polytomies = "none",
    logger = logger,
    expand_ingroup_to_full_donor = TRUE
  )

  expect_true(isTRUE(tpl$ok))
  expect_true(file.exists(logger$file))

  txt <- paste(readLines(logger$file), collapse = "\n")
  expect_match(txt, "Template: input summary")
  expect_match(txt, "Template: donor ultrametric conversion")
  expect_match(txt, "Template: output")
})

test_that("clade graft logs retained backbone tip in outgroup mode", {
  bb <- make_backbone()
  tpl <- make_template()

  res <- graft_template_onto_backbone(
    backbone_ultra = bb,
    placement = list(type = "clade", anchors = c("A", "B"), stem_mode = "outgroup"),
    template = tpl
  )

  expect_equal(as.character(res$log$retained_tip), "A")
})

test_that("clade graft logs retained backbone tip in crown mode", {
  bb <- make_backbone()
  tpl <- make_template()

  res <- graft_template_onto_backbone(
    backbone_ultra = bb,
    placement = list(type = "clade", anchors = c("A", "B"), stem_mode = "crown"),
    template = tpl
  )

  expect_equal(as.character(res$log$retained_tip), "A")
})