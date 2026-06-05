#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(yourpkg)  # <-- replace with your package name
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

to_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("1","true","t","yes","y")
}

to_num <- function(x, default) {
  if (is.null(x)) return(default)
  as.numeric(x)
}

# ---- parse args ----
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    if (grepl("^--", args[i])) {
      key <- sub("^--", "", args[i])
      if (i < length(args) && !grepl("^--", args[i+1])) {
        out[[key]] <- args[i+1]
        i <- i + 1
      } else {
        out[[key]] <- TRUE
      }
    }
    i <- i + 1
  }
  out
}

a <- parse_args(args)

if (is.null(a$tree)) {
  stop("Missing --tree")
}

tr <- ape::read.tree(a$tree)

mode <- a$mode %||% "single"

if (mode == "single") {

  plot_tree_autosize(
    tree = tr,
    pdf_path = a$out %||% paste0(a$tree, ".pdf"),
    cex = to_num(a$cex, 0.4),
    tight = to_bool(a$tight, FALSE)
  )

} else if (mode == "paged") {

  plot_tree_paged(
    tree = tr,
    out_prefix = a$out_prefix %||% "tree_paged",
    tips_per_page = as.integer(to_num(a$tips_per_page, 2500))
  )

} else {
  stop("Unknown mode")
}