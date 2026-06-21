library(testthat)
library(ape)
library(phytools)

# -------------------------------------------------------------------
# Fixtures
# -------------------------------------------------------------------

make_tree <- function() {
  read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
}

# Root of a rooted phylo
root_node <- function(tr) {
  setdiff(unique(tr$edge[, 1]), tr$edge[, 2])[1]
}

tip_desc_labels <- function(tree, node) {
  desc <- phytools::getDescendants(tree, node)
  tip_idx <- desc[desc <= ape::Ntip(tree)]
  tree$tip.label[tip_idx]
}

expect_sister_pair <- function(tree, tip1, tip2) {
  mrca <- ape::getMRCA(tree, c(tip1, tip2))
  expect_false(is.na(mrca))
  expect_equal(
    sort(tip_desc_labels(tree, mrca)),
    sort(c(tip1, tip2))
  )
}

expect_sister_to_clade <- function(tree, new_tip, clade_tips) {
  clade_mrca <- ape::getMRCA(tree, clade_tips)
  full_mrca <- ape::getMRCA(tree, c(new_tip, clade_tips))

  expect_false(is.na(clade_mrca))
  expect_false(is.na(full_mrca))
  expect_false(clade_mrca == full_mrca)

  kids <- tree$edge[tree$edge[, 1] == full_mrca, 2]
  kid_tips <- kids[kids <= ape::Ntip(tree)]
  kid_labels <- tree$tip.label[kid_tips]

  expect_true(new_tip %in% kid_labels)
  expect_true(any(kids == clade_mrca))
}

expect_tip_inside_clade <- function(tree, new_tip, clade_tips) {
  cl_mrca <- ape::getMRCA(tree, clade_tips)
  desc <- tip_desc_labels(tree, cl_mrca)
  expect_true(new_tip %in% desc)
}

# -------------------------------------------------------------------
# Internal helpers via namespace
# -------------------------------------------------------------------

ns <- asNamespace("MacroPhyloMaker")
incoming_edge_index <- get("incoming_edge_index", envir = ns)
draw_depth_fraction <- get("draw_depth_fraction", envir = ns)
tree_depth_helper <- get("tree_depth", envir = ns)
bind_tip_at <- get("bind_tip_at", envir = ns)

# -------------------------------------------------------------------
# Internal helper tests
# -------------------------------------------------------------------

test_that("incoming_edge_index finds the correct edge for a tip and an internal node", {
  tr <- make_tree()

  idx_tip_A <- match("A", tr$tip.label)
  ei_tip <- incoming_edge_index(tr, idx_tip_A)
  expect_equal(length(ei_tip), 1L)
  expect_equal(tr$edge[ei_tip, 2], idx_tip_A)

  mrca_AB <- ape::getMRCA(tr, c("A", "B"))
  ei_node <- incoming_edge_index(tr, mrca_AB)
  expect_equal(length(ei_node), 1L)
  expect_equal(tr$edge[ei_node, 2], mrca_AB)
})

test_that("draw_depth_fraction respects bounds and falls back safely", {
  set.seed(1)
  x <- replicate(50, draw_depth_fraction(
    shape1 = 2, shape2 = 5,
    min_frac = 0.2, max_frac = 0.4
  ))
  expect_true(all(x >= 0.2))
  expect_true(all(x <= 0.4))

  set.seed(1)
  y <- replicate(50, draw_depth_fraction(
    shape1 = -1, shape2 = 0,
    min_frac = 0.3, max_frac = 0.6
  ))
  expect_true(all(y >= 0.3))
  expect_true(all(y <= 0.6))

  # invalid interval should reset to [0,1]
  set.seed(1)
  z <- replicate(50, draw_depth_fraction(
    shape1 = 1, shape2 = 1,
    min_frac = 0.9, max_frac = 0.1
  ))
  expect_true(all(z >= 0))
  expect_true(all(z <= 1))
})

test_that("tree_depth helper equals max node height", {
  tr <- make_tree()
  expect_equal(tree_depth_helper(tr), max(phytools::nodeHeights(tr)))
})

test_that("bind_tip_at adds the tip and preserves requested attachment geometry", {
  tr <- make_tree()
  Htot <- tree_depth_helper(tr)

  # bind next to A at half its incident edge
  node_A <- match("A", tr$tip.label)
  tr2 <- bind_tip_at(tr, "X", where_node = node_A, position = 0.5, attach_age = 1.5)

  expect_true("X" %in% tr2$tip.label)
  expect_equal(length(tr2$tip.label), length(tr$tip.label) + 1)

  # New tip length should be Htot - attach_age = 0.5
  # Find the incoming edge for X after binding
  node_X <- match("X", tr2$tip.label)
  ei_X <- incoming_edge_index(tr2, node_X)
  expect_equal(length(ei_X), 1L)
  expect_equal(tr2$edge.length[ei_X], max(0, Htot - 1.5), tolerance = 1e-8)
})

# -------------------------------------------------------------------
# graft_sister_to_tip()
# -------------------------------------------------------------------

test_that("graft_sister_to_tip errors on missing target tip", {
  tr <- make_tree()

  expect_error(
    graft_sister_to_tip(tr, new_label = "X", sister_to = "ZZZ"),
    "not found"
  )
})

test_that("graft_sister_to_tip adds a tip as true sister and captures placement log", {
  tr <- make_tree()

  out <- graft_sister_to_tip(
    tr,
    new_label = "X",
    sister_to = "A",
    shape1 = 2,
    shape2 = 5,
    min_frac = 0.25,
    max_frac = 0.75,
    .capture = TRUE
  )

  expect_s3_class(out$tree, "phylo")
  expect_true("X" %in% out$tree$tip.label)
  expect_equal(length(out$tree$tip.label), 5L)
  expect_sister_pair(out$tree, "A", "X")

  lg <- out$log
  expect_true(is.list(lg))
  expect_true(is.finite(lg$edge_length))
  expect_true(is.finite(lg$chosen_frac))
  expect_true(lg$chosen_frac >= 0.25)
  expect_true(lg$chosen_frac <= 0.75)
  expect_true(is.finite(lg$position_from_child))
  expect_true(is.finite(lg$attach_age))
  expect_true(is.finite(lg$new_tip_length))
  expect_true(is.finite(lg$tree_depth))
})

# -------------------------------------------------------------------
# graft_sister_to_clade()
# -------------------------------------------------------------------

test_that("graft_sister_to_clade errors when fewer than two valid clade tips remain", {
  tr <- make_tree()

  expect_error(
    graft_sister_to_clade(tr, new_label = "X", clade_tips = c("A", "ZZZ")),
    "at least two existing tips"
  )
})

test_that("graft_sister_to_clade adds a tip as sister to MRCA-defined clade", {
  tr <- make_tree()

  out <- graft_sister_to_clade(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    shape1 = 2,
    shape2 = 2,
    min_frac = 0.1,
    max_frac = 0.9,
    .capture = TRUE
  )

  expect_s3_class(out$tree, "phylo")
  expect_true("X" %in% out$tree$tip.label)
  expect_equal(length(out$tree$tip.label), 5L)
  expect_sister_to_clade(out$tree, "X", c("A", "B"))

  lg <- out$log
  expect_true(is.list(lg))
  expect_true(is.finite(lg$chosen_frac))
  expect_true(lg$chosen_frac >= 0.1)
  expect_true(lg$chosen_frac <= 0.9)
})

test_that("graft_sister_to_clade handles root MRCA by using the root polytomy path", {
  tr <- make_tree()

  # Entire tree -> MRCA is root
  out <- graft_sister_to_clade(
    tr,
    new_label = "X",
    clade_tips = c("A", "B", "C", "D"),
    .capture = TRUE
  )

  expect_true("X" %in% out$tree$tip.label)
  expect_equal(length(out$tree$tip.label), 5L)
  expect_true(is.na(out$log$edge_index))
  expect_equal(out$log$where_node, root_node(tr))
})

# -------------------------------------------------------------------
# graft_within_clade_random()
# -------------------------------------------------------------------

test_that("graft_within_clade_random validates labels and duplicate new labels", {
  tr <- make_tree()

  expect_error(
    graft_within_clade_random(tr, new_label = "A", clade_tips = c("A", "B")),
    "already exists"
  )

  expect_error(
    graft_within_clade_random(tr, new_label = "X", clade_tips = c("A", "ZZZ")),
    "at least two existing tips"
  )
})

test_that("graft_within_clade_random places the new tip inside the target clade", {
  tr <- make_tree()

  out <- graft_within_clade_random(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    include_terminal = TRUE,
    include_stem = FALSE,
    shape1 = 2,
    shape2 = 5,
    min_frac = 0.2,
    max_frac = 0.8,
    return_log = TRUE
  )

  expect_s3_class(out$tree, "phylo")
  expect_true("X" %in% out$tree$tip.label)
  expect_equal(length(out$tree$tip.label), 5L)
  expect_tip_inside_clade(out$tree, "X", c("A", "B"))

  lg <- out$log
  expect_true(lg$chosen_frac >= 0.2)
  expect_true(lg$chosen_frac <= 0.8)
  expect_true(lg$edge_type %in% c("terminal", "internal", "stem"))
})

test_that("graft_within_clade_random can be forced to use only terminal edges", {
  tr <- make_tree()

  out <- graft_within_clade_random(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    include_terminal = TRUE,
    include_stem = FALSE,
    return_log = TRUE
  )

  # For cherry (A,B), inside-clade edges are terminal only
  expect_equal(out$log$edge_type, "terminal")
})

test_that("graft_within_clade_random can be forced to use only the stem edge", {
  tr <- make_tree()

  out <- graft_within_clade_random(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    include_terminal = FALSE,
    include_stem = TRUE,
    return_log = TRUE
  )

  expect_equal(out$log$edge_type, "stem")
})

test_that("graft_within_clade_random errors when no candidate edges remain", {
  tr <- make_tree()

  expect_error(
    graft_within_clade_random(
      tr,
      new_label = "X",
      clade_tips = c("A", "B"),
      include_terminal = FALSE,
      include_stem = FALSE
    ),
    "No internal edges|No candidate edges"
  )
})

test_that("graft_within_clade_random is reproducible under global seed", {
  tr <- make_tree()

  set.seed(42)
  out1 <- graft_within_clade_random(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    include_terminal = TRUE,
    include_stem = FALSE,
    return_log = TRUE
  )

  set.seed(42)
  out2 <- graft_within_clade_random(
    tr,
    new_label = "X",
    clade_tips = c("A", "B"),
    include_terminal = TRUE,
    include_stem = FALSE,
    return_log = TRUE
  )

  expect_equal(out1$log$chosen_edge_index, out2$log$chosen_edge_index)
  expect_equal(out1$log$chosen_frac, out2$log$chosen_frac)
  expect_equal(sort(out1$tree$tip.label), sort(out2$tree$tip.label))
})
