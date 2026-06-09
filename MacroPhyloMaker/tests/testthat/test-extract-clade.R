library(testthat)
library(ape)

make_basic_tree <- function() {
  read.tree(text = "((Camponotus_a:1,Camponotus_b:1):1,(Formica_c:1,Lasius_d:1):1);")
}

make_nonmono_tree <- function() {
  # Camponotus not monophyletic
  read.tree(text = "((Camponotus_a:1,Formica_c:1):1,(Camponotus_b:1,Lasius_d:1):1);")
}

test_that("genus extraction returns correct tips", {
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
  expect_equal(length(res$ingroup_tips), 2)
  expect_true(all(res$ingroup_tips %in% res$tree$tip.label))
})

test_that("MRCA extraction returns correct clade", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    mrca_tips = c("Camponotus_a", "Camponotus_b"),
    outgroup = "none",
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  expect_equal(sort(res$ingroup_tips),
               sort(c("Camponotus_a", "Camponotus_b")))
})

test_that("must supply exactly one of genus or mrca_tips", {
  tr <- make_basic_tree()

  expect_error(
    extract_clade_with_outgroup(tr, write_tree = FALSE),
    "Provide either"
  )

  expect_error(
    extract_clade_with_outgroup(tr, genus = "Camponotus",
                               mrca_tips = c("A")),
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

test_that("non-monophyletic genus prunes extras", {
  tr <- make_nonmono_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    nonmono = "prune_extras",
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  expect_true(all(grepl("^Camponotus_", res$tree$tip.label)))
})

test_that("sister_one outgroup is added correctly", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    outgroup = "sister_one",
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  expect_false(is.na(res$outgroup))
  expect_true(res$outgroup %in% c("Formica_c", "Lasius_d"))
  expect_equal(length(res$tree$tip.label), length(res$ingroup_tips) + 1)
})

test_that("outgroup none returns only ingroup", {
  tr <- make_basic_tree()

  res <- extract_clade_with_outgroup(
    tr,
    genus = "Camponotus",
    outgroup = "none",
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  expect_true(is.na(res$outgroup))
  expect_equal(length(res$tree$tip.label), length(res$ingroup_tips))
})

test_that("MRCA matching uses original labels, not cleaned labels", {
  tr <- read.tree(text = "((Camponotus_AAA:1,Camponotus_BBB:1):1,(X:1,Y:1):1);")

  # These would collapse to Camponotus_aaa / bbb if cleaned BEFORE matching
  res <- extract_clade_with_outgroup(
    tr,
    mrca_tips = c("Camponotus_AAA", "Camponotus_BBB"),
    clean = "genus_species",
    outgroup = "none",
    write_tree = FALSE
  )

  expect_equal(length(res$ingroup_tips), 2)
})

test_that("duplicates collapse in genus mode", {
  tr <- read.tree(text = "((A:1,B:1):1,(X:1,Y:1):1);")

  # Force duplicate binomials
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
    write_tree = FALSE,
    write_renames = FALSE,
    write_drops = FALSE
  )

  # Only one species should remain after collapsing
  expect_equal(length(res$ingroup_tips), 1)
  expect_equal(length(grep("^Camponotus_a", res$ingroup_tips)), 1)
})

test_that("logger writes expected sections", {
  tr <- make_basic_tree()
  dir <- tempdir()

  res <- extract_clade_with_outgroup(
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