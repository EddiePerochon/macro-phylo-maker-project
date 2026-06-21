library(testthat)
library(ape)

# -------------------------------
# Helpers
# -------------------------------

make_tree <- function() {
  read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
}

make_plan <- function(fun, tip, sister = NA_character_) {
  data.frame(
    Function = fun,
    GraftedTip = tip,
    Sister = sister,
    stringsAsFactors = FALSE
  )
}

# ---- topological helpers ----

expect_sister_tip <- function(tree, tip1, tip2) {
  mrca <- ape::getMRCA(tree, c(tip1, tip2))
  children <- tree$edge[tree$edge[, 1] == mrca, 2]

  child_tips <- children[children <= length(tree$tip.label)]
  labels <- tree$tip.label[child_tips]

  expect_equal(sort(labels), sort(c(tip1, tip2)))
}

expect_sister_clade <- function(tree, new_tip, clade_tips) {
  clade_mrca <- ape::getMRCA(tree, clade_tips)
  full_mrca <- ape::getMRCA(tree, c(new_tip, clade_tips))

  expect_false(clade_mrca == full_mrca)

  children <- tree$edge[tree$edge[, 1] == full_mrca, 2]

  child_tips <- children[children <= length(tree$tip.label)]
  child_labels <- tree$tip.label[child_tips]

  expect_true(new_tip %in% child_labels)
  expect_true(any(children == clade_mrca))
}

expect_within_clade <- function(tree, tip, clade_tips) {
  clade_mrca <- ape::getMRCA(tree, clade_tips)
  tip_mrca <- ape::getMRCA(tree, c(tip, clade_tips[1]))

  # tip MRCA must be inside or equal to clade MRCA
  expect_true(tip_mrca >= clade_mrca)
}

# -------------------------------
# Error handling
# -------------------------------

test_that("unknown graft function errors", {
  tr <- make_tree()
  plan <- make_plan("not_a_function", "X")

  expect_error(
    apply_grafts_from_table(tr, plan),
    "Unrecognized Function value"
  )
})

# -------------------------------
# Sister-to-tip
# -------------------------------

test_that("graft sister to tip adds tip and places correctly", {
  tr <- make_tree()

  plan <- make_plan("graft_sister_to_tip", "X", "A")

  res <- apply_grafts_from_table(tr, plan)
  tree <- res$tree

  # existence
  expect_true("X" %in% tree$tip.label)

  # count increases by 1
  expect_equal(length(tree$tip.label), length(tr$tip.label) + 1)

  # topology
  expect_sister_tip(tree, "A", "X")
})

# -------------------------------
# Sister-to-clade
# -------------------------------

test_that("graft sister to clade places tip correctly", {
  tr <- make_tree()

  plan <- make_plan("graft_sister_to_clade", "X", "A;B")

  res <- apply_grafts_from_table(tr, plan)
  tree <- res$tree

  expect_true("X" %in% tree$tip.label)

  expect_sister_clade(tree, "X", c("A", "B"))
})

# -------------------------------
# Within-clade random
# -------------------------------

test_that("within-clade graft places tip inside clade", {
  tr <- make_tree()

  plan <- make_plan("graft_within_clade_random", "X", "A;B")

  res <- apply_grafts_from_table(tr, plan, seed = 1)
  tree <- res$tree

  expect_true("X" %in% tree$tip.label)

  expect_within_clade(tree, "X", c("A", "B"))
})

# -------------------------------
# Parameters
# -------------------------------

test_that("branch position parameters accepted silently", {
  tr <- make_tree()

  plan <- make_plan("graft_sister_to_tip", "X", "A")

  expect_silent(
    apply_grafts_from_table(
      tr, plan,
      default_shape1 = 2,
      default_shape2 = 5,
      default_min_frac = 0.1,
      default_max_frac = 0.9
    )
  )
})

# -------------------------------
# Ingroup extraction
# -------------------------------

test_that("extract ingroup by anchors works", {
  tr <- make_tree()

  res <- extract_ingroup_by_anchors(tr, c("A", "B"))

  expect_equal(sort(res$tree$tip.label), c("A", "B"))
})

test_that("ingroup extraction handles missing anchors", {
  tr <- make_tree()

  res <- extract_ingroup_by_anchors(tr, c("A", "NOT"))

  expect_equal(res$tree$tip.label, tr$tip.label)
})

test_that("ingroup extraction requires >=2 anchors", {
  tr <- make_tree()

  res <- extract_ingroup_by_anchors(tr, "A")

  expect_equal(res$tree$tip.label, tr$tip.label)
})

# -------------------------------
# Utilities
# -------------------------------

test_that("parse_sister handles separators", {
  x <- "A;B, C | D"

  res <- parse_sister(x)

  expect_equal(sort(res), c("A", "B", "C", "D"))
})

test_that("normalize_fun cleans names", {
  expect_equal(normalize_fun("Graft_Sister_To_Tip"), "graftsistertotip")
})

test_that("get_num_or returns default correctly", {
  df <- data.frame(x = "bad")

  expect_equal(get_num_or(df, "x", 5), 5)
})

# -------------------------------
# Logging
# -------------------------------

test_that("apply_grafts returns log when requested", {
  tr <- make_tree()
  plan <- make_plan("graft_sister_to_tip", "X", "A")

  res <- apply_grafts_from_table(tr, plan, return_log = TRUE)

  expect_true(is.list(res))
  expect_true("log" %in% names(res))
  expect_true(nrow(res$log) == 1)
})

# -------------------------------
# Output writing
# -------------------------------

test_that("pipeline_write_outputs creates files", {
  tr <- make_tree()
  log <- data.frame(step = 1)

  dir <- tempdir()
  prefix <- file.path(dir, "test")

  paths <- pipeline_write_outputs(tr, log, prefix, plot_pdf = FALSE)

  expect_true(file.exists(paths$tree))
  expect_true(file.exists(paths$log))
})
