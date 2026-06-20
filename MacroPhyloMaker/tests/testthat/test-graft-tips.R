library(testthat)
library(ape)

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

test_that("unknown graft function errors", {
  tr <- make_tree()
  plan <- make_plan("not_a_function", "X")

  expect_error(
    apply_grafts_from_table(tr, plan),
    "Unrecognized Function value"
  )
})

test_that("graft next to tip adds a new taxon", {
  tr <- make_tree()

  plan <- data.frame(
    Function = "graft_sister_to_tip",
    GraftedTip = "X",
    Sister = "A",
    stringsAsFactors = FALSE
  )

  res <- apply_grafts_from_table(tr, plan)

  expect_true("X" %in% res$tree$tip.label)
  expect_equal(length(res$tree$tip.label), 5)
})

test_that("within clade graft places tip inside clade", {
  tr <- make_tree()

  plan <- data.frame(
    Function = "graft_within_clade_random",
    GraftedTip = "X",
    Sister = "A;B",
    stringsAsFactors = FALSE
  )

  res <- apply_grafts_from_table(tr, plan, seed = 1)

  expect_true("X" %in% res$tree$tip.label)
})

test_that("branch position parameters accepted", {
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

test_that("apply_grafts returns log when requested", {
  tr <- make_tree()
  plan <- make_plan("graft_sister_to_tip", "X", "A")

  res <- apply_grafts_from_table(tr, plan, return_log = TRUE)

  expect_true(is.list(res))
  expect_true("log" %in% names(res))
})

test_that("pipeline_write_outputs creates files", {
  tr <- make_tree()
  log <- data.frame(step = 1)

  dir <- tempdir()
  prefix <- file.path(dir, "test")

  paths <- pipeline_write_outputs(tr, log, prefix, plot_pdf = FALSE)

  expect_true(file.exists(paths$tree))
  expect_true(file.exists(paths$log))
})