library(testthat)
library(ape)

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

make_tree_small <- function() {
  read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
}

make_tree_large <- function() {
  read.tree(text = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
}

make_tree_long_labels <- function() {
  tr <- make_tree_small()
  tr$tip.label <- c(
    "Short",
    "Medium_label",
    "A_very_very_long_tip_label_that_should_increase_width",
    "Another_long_label"
  )
  tr
}

# -------------------------------------------------------------------
# suggest_pdf_size()
# -------------------------------------------------------------------

test_that("suggest_pdf_size returns bounded numeric dimensions", {
  tr <- make_tree_small()

  dims <- suggest_pdf_size(tr)

  expect_type(dims$width, "double")
  expect_type(dims$height, "double")
  expect_true(is.finite(dims$width))
  expect_true(is.finite(dims$height))
  expect_gte(dims$width, 5)
  expect_lte(dims$width, 20)
  expect_gte(dims$height, 6)
  expect_lte(dims$height, 200)
})

test_that("suggest_pdf_size increases height with more tips before clamping", {
  tr_small <- make_tree_small()
  tr_large <- make_tree_large()

  dims_small <- suggest_pdf_size(tr_small, min_h = 0)
  dims_large <- suggest_pdf_size(tr_large, min_h = 0)

  expect_gt(dims_large$height, dims_small$height)
})

test_that("suggest_pdf_size increases width for longer labels when auto_width is TRUE", {
  tr_short <- make_tree_small()
  tr_long  <- make_tree_long_labels()

  dims_short <- suggest_pdf_size(tr_short, auto_width = TRUE)
  dims_long  <- suggest_pdf_size(tr_long, auto_width = TRUE)

  expect_gt(dims_long$width, dims_short$width)
})

test_that("suggest_pdf_size respects fixed width when auto_width is FALSE", {
  tr <- make_tree_long_labels()

  dims <- suggest_pdf_size(
    tr,
    pdf_width = 9,
    auto_width = FALSE
  )

  expect_equal(dims$width, 9)
})

test_that("suggest_pdf_size errors for invalid cex", {
  tr <- make_tree_small()

  expect_error(
    suggest_pdf_size(tr, cex = 0),
    class = "simpleError"
  )
})

# -------------------------------------------------------------------
# plot_tree_autosize()
# -------------------------------------------------------------------

test_that("plot_tree_autosize writes a PDF and returns final dimensions", {
  tr <- make_tree_small()
  pdf_path <- tempfile(fileext = ".pdf")

  dims <- plot_tree_autosize(
    tr,
    pdf_path = pdf_path,
    pdf_auto = TRUE
  )

  expect_true(file.exists(pdf_path))
  expect_true(file.info(pdf_path)$size > 0)
  expect_true(is.list(dims))
  expect_true(all(c("width", "height") %in% names(dims)))
  expect_true(is.finite(dims$width))
  expect_true(is.finite(dims$height))
})

test_that("plot_tree_autosize uses plot_fn when supplied", {
  tr <- make_tree_small()
  pdf_path <- tempfile(fileext = ".pdf")
  seen_n_tip <- NULL

  custom_plot <- function(x) {
    seen_n_tip <<- ape::Ntip(x)
    graphics::plot(x, cex = 0.5)
  }

  plot_tree_autosize(
    tr,
    pdf_path = pdf_path,
    plot_fn = custom_plot
  )

  expect_equal(seen_n_tip, ape::Ntip(tr))
  expect_true(file.exists(pdf_path))
})

test_that("plot_tree_autosize works in tight mode", {
  tr <- make_tree_small()
  pdf_path <- tempfile(fileext = ".pdf")

  expect_no_error(
    plot_tree_autosize(
      tr,
      pdf_path = pdf_path,
      tight = TRUE
    )
  )

  expect_true(file.exists(pdf_path))
})

test_that("plot_tree_autosize respects explicit dimensions when pdf_auto is FALSE", {
  tr <- make_tree_small()
  pdf_path <- tempfile(fileext = ".pdf")

  dims <- plot_tree_autosize(
    tr,
    pdf_path = pdf_path,
    pdf_width = 7,
    pdf_height = 11,
    pdf_auto = FALSE
  )

  expect_equal(dims$width, 7)
  expect_equal(dims$height, 11)
})

# -------------------------------------------------------------------
# plot_tree_paged()
# -------------------------------------------------------------------

test_that("plot_tree_paged writes the expected number of page PDFs", {
  tr <- make_tree_large()
  prefix <- tempfile("paged_tree_")

  n_pages <- plot_tree_paged(
    tr,
    out_prefix = prefix,
    tips_per_page = 2
  )

  expect_equal(n_pages, 3)

  expected_files <- sprintf("%s_page%02d.pdf", prefix, 1:3)
  expect_true(all(file.exists(expected_files)))
  expect_true(all(file.info(expected_files)$size > 0))
})

test_that("plot_tree_paged works when pdf_auto is FALSE", {
  tr <- make_tree_large()
  prefix <- tempfile("paged_tree_fixed_")

  n_pages <- plot_tree_paged(
    tr,
    out_prefix = prefix,
    tips_per_page = 3,
    pdf_auto = FALSE
  )

  expect_equal(n_pages, 2)
  expected_files <- sprintf("%s_page%02d.pdf", prefix, 1:2)
  expect_true(all(file.exists(expected_files)))
})

test_that("plot_tree_paged errors for invalid tips_per_page", {
  tr <- make_tree_small()

  expect_error(
    plot_tree_paged(tr, out_prefix = tempfile("bad_"), tips_per_page = 0),
    class = "simpleError"
  )
})
