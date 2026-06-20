library(testthat)
library(ape)
library(phytools)

# -------------------------------
# Helpers
# -------------------------------

make_basic_tree <- function() {
  read.tree(text = "((Camponotus_a:1,Camponotus_b:1):1,(Formica_c:1,Lasius_d:1):1);")
}

make_nonmono_tree <- function() {
  read.tree(text = "((Camponotus_a:1,Formica_c:1):1,(Camponotus_b:1,Lasius_d:1):1);")
}

# ---- topology helpers ----

expect_exact_clade <- function(tree, tips) {
  mrca <- ape::getMRCA(tree, tips)

  descendants <- tree$edge[tree$edge[,1] == mrca, 2]

  # collect all descendant tips
  desc <- phytools::getDescendants(tree, mrca)
  tip_idx <- desc[desc <= length(tree$tip.label)]
  labels <- tree$tip.label[tip_idx]

  expect_equal(sort(labels), sort(tips))
}

expect_contains_only_genus <- function(tree, genus) {
  expect_true(all(grepl(paste0("^", genus, "_"), tree$tip.label)))
}

# -------------------------------
# Genus extraction
# -------------------------------

test_that("genus extraction returns correct clade", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    outgroup = "none",
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  expect_true(all(grepl("^Camponotus_", res$ingroup_tips)))
  expect_exact_clade(res$tree, res$ingroup_tips)
})

# -------------------------------
# MRCA extraction
# -------------------------------

test_that("MRCA extraction returns correct clade", {
  tr <- make_basic_tree()

  anchors <- c("Camponotus_a", "Camponotus_b")

  res <- extract_clade_with_outgroup(
    tr,
    mrca_tips = anchors,
    outgroup = "none",
    write_tree = FALSE
  )

  expect_exact_clade(res$tree, anchors)
})

# -------------------------------
# Argument validation
# -------------------------------

test_that("must supply exactly one of genus or mrca_tips", {
  tr <- make_basic_tree()

  expect_error(
    extract_clade_with_outgroup(tr, write_tree = FALSE),
    "Provide either"
  )

  expect_error(
    extract_clade_with_outgroup(tr,
      genus = "Camponotus",
      mrca_tips = c("A")
    ),
    "Provide only one"
  )
})

test_that("invalid MRCA tips error", {
  tr <- make_basic_tree()

  expect_error(
    extract_clade_with_outgroup(
      tr,
      mrca_tips = c("NOT_REAL"),
      write_tree = FALSE
    ),
    "not in tree"
  )
})

# -------------------------------
# Non-monophyly handling
# -------------------------------

test_that("non-monophyletic genus errors when requested", {
  tr <- make_nonmono_tree()

  expect_error(
    extract_clade_with_outgroup(
      tr,
      genus = "Camponotus",
      nonmono = "error",
      write_tree = FALSE
    ),
    "not monophyletic"
  )
})

test_that("non-monophyletic genus prunes extras but preserves genus", {
  tr <- make_nonmono_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    nonmono = "prune_extras",
    write_tree = FALSE
  )

  expect_contains_only_genus(res$tree, "Camponotus")
})

# -------------------------------
# Outgroup logic
# -------------------------------

test_that("sister_one outgroup is valid and outside ingroup", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    outgroup = "sister_one",
    write_tree = FALSE
  )

  expect_false(is.na(res$outgroup))
  expect_true(res$outgroup %in% setdiff(tr$tip.label, res$ingroup_tips))

  # topology: MRCA of ingroup must not include outgroup
  ingroup_mrca <- ape::getMRCA(res$tree, res$ingroup_tips)
  all_tips <- phytools::getDescendants(res$tree, ingroup_mrca)
  labels <- res$tree$tip.label[all_tips]

  expect_false(res$outgroup %in% labels)
})

test_that("outgroup none returns only ingroup", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    outgroup = "none",
    write_tree = FALSE
  )

  expect_true(is.na(res$outgroup))
  expect_exact_clade(res$tree, res$ingroup_tips)
})

# -------------------------------
# Label handling
# -------------------------------

test_that("MRCA matching uses original labels, not cleaned labels", {
  tr <- read.tree(text = "((Camponotus_AAA:1,Camponotus_BBB:1):1,(X:1,Y:1):1);")

  res <- extract_clade_with_outgroup(
    tr,
    mrca_tips = c("Camponotus_AAA", "Camponotus_BBB"),
    clean = "genus_species",
    outgroup = "none",
    write_tree = FALSE
  )

  expect_exact_clade(res$tree, c("Camponotus_AAA", "Camponotus_BBB"))
})

# -------------------------------
# Duplicate collapse
# -------------------------------

test_that("duplicates collapse correctly in genus mode", {
  tr <- read.tree(text = "((A:1,B:1):1,(X:1,Y:1):1);")

  tr$tip.label <- c(
    "Camponotus_a",
    "Camponotus_a_extra",
    "X",
    "Y"
  )

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    clean = "genus_species",
    outgroup = "none",
    write_tree = FALSE
  )

  expect_equal(length(res$ingroup_tips), 1)
  expect_true(grepl("^Camponotus_a", res$ingroup_tips))
})

# -------------------------------
# Logging
# -------------------------------

test_that("logger writes expected structured sections", {
  tr <- make_basic_tree()
  dir <- tempdir()

  extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    out_dir = dir,
    write_tree = TRUE
  )

  log_files <- list.files(dir, pattern = "extract_.*\\.log", full.names = TRUE)
  expect_true(length(log_files) >= 1)

  log_content <- paste(readLines(log_files[1]), collapse = "\n")

  expect_match(log_content, "Settings")
  expect_match(log_content, "Input")
  expect_match(log_content, "Output")
  expect_match(log_content, "tips_in")
  expect_match(log_content, "tips_out")
})