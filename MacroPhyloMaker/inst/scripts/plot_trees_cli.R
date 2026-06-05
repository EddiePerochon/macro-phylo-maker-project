#!/usr/bin/env Rscript
# ============================================================
# Tree plotting utilities (stand-alone, CLI)
# - suggest_pdf_size(): same heuristic as pipeline
# - plot_tree_autosize(): if tight=TRUE, use pipeline's single-page layout
# - plot_tree_paged(): paginated plotting
# Requires: ape
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required but not installed. Install with: install.packages('ape')", call. = FALSE)
  }
  library(ape)
})

# --- 1) Same autosize helper as in species-grafting-realms.sh ---
suggest_pdf_size <- function(tree,
                             pdf_width  = 6,
                             cex        = 0.40,
                             target_tpi = 15,  # tips per inch
                             min_h      = 6,
                             max_h      = 200,
                             buffer_in  = 2,
                             auto_width = TRUE,
                             char_per_in_at_cex = 10,
                             min_w = 5,
                             max_w = 20) {
  stopifnot(inherits(tree, "phylo"))
  ntip <- ape::Ntip(tree)
  # Effective tips per inch scales with cex (pipeline's original rule):
  tpi_eff <- target_tpi * (0.40 / cex)^0.5
  core_h  <- ntip / max(1, tpi_eff)
  h_in    <- core_h + buffer_in
  h_in    <- max(min_h, min(max_h, h_in))

  if (isTRUE(auto_width)) {
    labs       <- as.character(tree$tip.label)
    max_chars  <- if (length(labs)) max(nchar(labs), na.rm = TRUE) else 12
    cpi_eff    <- char_per_in_at_cex * (0.40 / cex)  # smaller cex -> more chars/in
    label_w_in <- max_chars / max(1, cpi_eff)
    w_in       <- label_w_in + 2
    pdf_width  <- max(min_w, min(max_w, w_in))
  }

  list(width = pdf_width, height = h_in)
}

# --- 2) Single-page plotting with "pipeline-tight" fallback ---
plot_tree_autosize <- function(tree,
                               pdf_path,
                               cex          = 0.40,
                               pdf_width    = NULL,
                               pdf_height   = NULL,
                               pdf_auto     = TRUE,
                               target_tpi   = 15,
                               min_h        = 6,
                               max_h        = 200,
                               buffer_in    = 2,
                               auto_width   = TRUE,
                               char_per_in_at_cex = 10,
                               min_w        = 5,
                               max_w        = 20,
                               plot_fn      = NULL,   # optional function(tr) to plot custom
                               ladderize_it = TRUE,
                               # when TRUE, use the exact pipeline layout
                               tight        = FALSE,
                               # when tight=FALSE, use these (margin-safe) defaults:
                               mar_safe     = c(4, 3.5, 3, 1) + 0.1,
                               oma_safe     = c(0, 0, 0, 0),
                               ...) {
  stopifnot(inherits(tree, "phylo"))
  tr <- if (isTRUE(ladderize_it)) ladderize(tree) else tree

  # Default size if user supplies nothing and pdf_auto is FALSE: 6x20
  w <- if (is.null(pdf_width)) 6 else as.numeric(pdf_width)
  h <- if (is.null(pdf_height)) 20 else as.numeric(pdf_height)

  # Autosize (same as pipeline)
  if (isTRUE(pdf_auto) && (is.null(pdf_width) || is.null(pdf_height))) {
    dims <- suggest_pdf_size(tr,
                             pdf_width  = w,
                             cex        = cex,
                             target_tpi = target_tpi,
                             min_h      = min_h,
                             max_h      = max_h,
                             buffer_in  = buffer_in,
                             auto_width = auto_width,
                             char_per_in_at_cex = char_per_in_at_cex,
                             min_w      = min_w,
                             max_w      = max_w)
    if (is.null(pdf_height)) h <- dims$height
    if (is.null(pdf_width))  w <- dims$width
  }

  pdf(pdf_path, width = w, height = h)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  if (isTRUE(tight)) {
    # --- PIPELINE MODE (match species script) ---
    op <- par(mar = c(1, 1, 1, 1)); on.exit(par(op), add = TRUE)
    if (is.function(plot_fn)) {
      plot_fn(tr)
    } else {
      plot(tr, cex = cex)  # same call as in pipeline
    }
  } else {
    # --- MARGIN-SAFE MODE (standard margins, no clipping) ---
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = mar_safe, oma = oma_safe)
    if (is.function(plot_fn)) {
      plot_fn(tr)
    } else {
      plot(tr, cex = cex, no.margin = FALSE, ...)  # let ape manage margins
    }
  }

  invisible(list(width = w, height = h))
}

# --- 3) Multi-page plotting ---
plot_tree_paged <- function(tree,
                            out_prefix    = "tree_paged",
                            tips_per_page = 2500,
                            cex           = 0.35,
                            pdf_auto      = TRUE,
                            target_tpi    = 16,
                            min_h         = 6,
                            max_h         = 200,
                            buffer_in     = 1.5,
                            auto_width    = TRUE,
                            char_per_in_at_cex = 10,
                            min_w         = 5,
                            max_w         = 20,
                            ladderize_it  = TRUE,
                            ...) {
  stopifnot(inherits(tree, "phylo"))
  tr   <- if (isTRUE(ladderize_it)) ladderize(tree) else tree
  n    <- ape::Ntip(tr)
  idx  <- seq_len(n)
  pages <- split(idx, ceiling(idx / tips_per_page))

  for (k in seq_along(pages)) {
    block  <- pages[[k]]
    keep   <- tr$tip.label[block]
    subtr  <- ape::keep.tip(tr, keep)

    dims <- if (pdf_auto) {
      suggest_pdf_size(subtr,
                       pdf_width  = 6,
                       cex        = cex,
                       target_tpi = target_tpi,
                       min_h      = min_h,
                       max_h      = max_h,
                       buffer_in  = buffer_in,
                       auto_width = auto_width,
                       char_per_in_at_cex = char_per_in_at_cex,
                       min_w      = min_w,
                       max_w      = max_w)
    } else list(width = 6, height = 20)

    pdf_path <- sprintf("%s_page%02d.pdf", out_prefix, k)
    pdf(pdf_path, width = dims$width, height = dims$height)

    op <- par(mar = c(1, 1, 1, 1)); on.exit(par(op), add = TRUE)
    plot(subtr, cex = cex, ...)
    dev.off()
  }
  invisible(length(pages))
}

# -----------------------------
# CLI utilities
# -----------------------------
to_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("1","true","t","yes","y","on")) return(TRUE)
  if (x %in% c("0","false","f","no","n","off")) return(FALSE)
  # bare flag presence like --tight with no value:
  if (identical(x, "")) return(TRUE)
  default
}

to_num <- function(x, default = NA_real_) {
  if (is.null(x)) return(default)
  as.numeric(x)
}

to_vec_num <- function(x, n = NULL, default = NULL) {
  if (is.null(x) || x == "") return(default)
  v <- as.numeric(strsplit(x, "[, ]+")[[1]])
  if (!is.null(n) && length(v) != n) {
    stop(sprintf("Expected %d numeric values, got: %s", n, paste(v, collapse = ",")))
  }
  v
}

basename_noext <- function(path) {
  b <- basename(path)
  sub("\\.[^.]*$", "", b)
}

print_help <- function() {
  cat("
plot_trees_cli.R — plot Newick trees to PDF (single page or paged)

USAGE:
  Rscript plot_trees_cli.R --tree <path> [--mode single|paged] [options]

REQUIRED:
  --tree <path>             Path to Newick/Nexus tree file.

COMMON OPTIONS:
  --mode single|paged       Plot mode (default: single)
  --cex <num>               Tip label cex (default: 0.40 single, 0.35 paged)
  --ladderize <bool>        Ladderize tree (default: TRUE)
  --pdf_auto <bool>         Auto-size height/width (default: TRUE)
  --auto_width <bool>       Auto-width based on label lengths (default: TRUE)
  --min_h <num>             Minimum PDF height in inches (default: 6)
  --max_h <num>             Maximum PDF height in inches (default: 200)
  --buffer_in <num>         Extra inches beyond computed core height
                            (default: 2 for single, 1.5 for paged)
  --char_per_in_at_cex <n>  Characters per inch at cex=0.40 (default: 10)
  --min_w <num>             Min PDF width when auto_width=TRUE (default: 5)
  --max_w <num>             Max PDF width when auto_width=TRUE (default: 20)

SINGLE-PAGE OPTIONS:
  --out <path>              Output PDF path (default: <tree>.pdf)
  --tight <bool>            Pipeline-tight layout (par(mar=1,1,1,1)) (default: TRUE)
  --width <num>             PDF width in inches (overrides autosizing if set)
  --height <num>            PDF height in inches (overrides autosizing if set)
  --mar \"a,b,c,d\"         Margins when --tight=FALSE (default: 4,3.5,3,1)
  --oma \"a,b,c,d\"         Outer margins when --tight=FALSE (default: 0,0,0,0)

PAGED OPTIONS:
  --out_prefix <str>        Prefix for per-page PDFs (default: <tree_basename>_pages)
  --tips_per_page <int>     Tips per page (default: 2500)

EXAMPLES:
  # Single-page, pipeline-tight look, autosized
  Rscript plot_trees_cli.R --tree ../data/Formicidae.tacted.newick.tre --mode single --cex 0.10 --tight true

  # Single-page, margin-safe mode with manual size
  Rscript plot_trees_cli.R --tree ants.tre --mode single --tight false --width 8 --height 30 --cex 0.2

  # Paged, 3000 tips per page
  Rscript plot_trees_cli.R --tree ants.tre --mode paged --tips_per_page 3000 --cex 0.10

NOTES:
  * Install 'ape' if missing: install.packages('ape')
  * Boolean values accept: true/false, 1/0, yes/no, on/off
  * Defaults mirror your pipeline's heuristics.
\n")
}

# Minimal argument parser (supports --key=value or --key value, and bare --flag)
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
    print_help()
    quit(save = "no", status = 0)
  }
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--", a)) {
      kv <- sub("^--", "", a)
      if (grepl("=", kv)) {
        key <- sub("=.*$", "", kv)
        val <- sub("^[^=]*=", "", kv)
        out[[key]] <- val
      } else {
        # bare flag or a separate value follows
        if (i < length(args) && !grepl("^--", args[[i+1]])) {
          out[[kv]] <- args[[i+1]]
          i <- i + 1
        } else {
          out[[kv]] <- ""  # presence-only flag (treated as TRUE)
        }
      }
    } else {
      warning("Ignoring positional argument: ", a)
    }
    i <- i + 1
  }
  out
}

safe_read_tree <- function(path) {
  if (!file.exists(path)) stop("Tree file not found: ", path, call. = FALSE)
  # Try Newick first, then Nexus
  tr <- tryCatch(ape::read.tree(path), error = function(e) NULL)
  if (inherits(tr, "phylo")) return(tr)
  tr <- tryCatch(ape::read.nexus(path), error = function(e) NULL)
  if (inherits(tr, "phylo")) return(tr)
  stop("Could not read tree. Ensure it's Newick (.tre/.tree/.nwk) or Nexus.", call. = FALSE)
}

# -----------------------------
# Main
# -----------------------------
main <- function() {
  a <- parse_args()

  # Required
  tree_path <- a[["tree"]]
  if (is.null(tree_path)) {
    stop("Missing required --tree <path>. Use --help for usage.", call. = FALSE)
  }
  mode <- tolower(a[["mode"]] %||% "single")
  if (!mode %in% c("single","paged")) stop("--mode must be 'single' or 'paged'")

  tr <- safe_read_tree(tree_path)

  # Common options
  ladderize_it <- to_bool(a[["ladderize"]], TRUE)
  pdf_auto     <- to_bool(a[["pdf_auto"]], TRUE)
  auto_width   <- to_bool(a[["auto_width"]], TRUE)
  min_h        <- to_num(a[["min_h"]], 6)
  max_h        <- to_num(a[["max_h"]], 200)
  char_per_in  <- to_num(a[["char_per_in_at_cex"]], 10)
  min_w        <- to_num(a[["min_w"]], 5)
  max_w        <- to_num(a[["max_w"]], 20)

  if (mode == "single") {
    cex         <- to_num(a[["cex"]], 0.40)
    tight       <- to_bool(a[["tight"]], TRUE)
    width       <- to_num(a[["width"]], NA_real_)
    height      <- to_num(a[["height"]], NA_real_)
    target_tpi  <- to_num(a[["target_tpi"]], 15)
    buffer_in   <- to_num(a[["buffer_in"]], 2)
    mar_safe    <- to_vec_num(a[["mar"]], n = 4, default = c(4, 3.5, 3, 1) + 0.1)
    oma_safe    <- to_vec_num(a[["oma"]], n = 4, default = c(0, 0, 0, 0))
    out_path    <- a[["out"]]
    if (is.null(out_path) || out_path == "") {
      out_path <- paste0(tree_path, ".pdf")
    }

    message("Single-page plotting → ", out_path)
    dims <- plot_tree_autosize(
      tree        = tr,
      pdf_path    = out_path,
      cex         = cex,
      pdf_width   = if (is.na(width)) NULL else width,
      pdf_height  = if (is.na(height)) NULL else height,
      pdf_auto    = pdf_auto,
      target_tpi  = target_tpi,
      min_h       = min_h,
      max_h       = max_h,
      buffer_in   = buffer_in,
      auto_width  = auto_width,
      char_per_in_at_cex = char_per_in,
      min_w       = min_w,
      max_w       = max_w,
      ladderize_it = ladderize_it,
      tight       = tight,
      mar_safe    = mar_safe,
      oma_safe    = oma_safe
    )
    message(sprintf("Done. PDF size: width=%.2f in, height=%.2f in", dims$width, dims$height))

  } else if (mode == "paged") {
    cex            <- to_num(a[["cex"]], 0.35)
    target_tpi     <- to_num(a[["target_tpi"]], 16)
    buffer_in      <- to_num(a[["buffer_in"]], 1.5)
    tips_per_page  <- as.integer(to_num(a[["tips_per_page"]], 2500))
    out_prefix     <- a[["out_prefix"]]
    if (is.null(out_prefix) || out_prefix == "") {
      out_prefix <- paste0(basename_noext(tree_path), "_pages")
    }

    message(sprintf("Paged plotting → prefix '%s', tips/page=%d", out_prefix, tips_per_page))
    npages <- plot_tree_paged(
      tree        = tr,
      out_prefix  = out_prefix,
      tips_per_page = tips_per_page,
      cex         = cex,
      pdf_auto    = pdf_auto,
      target_tpi  = target_tpi,
      min_h       = min_h,
      max_h       = max_h,
      buffer_in   = buffer_in,
      auto_width  = auto_width,
      char_per_in_at_cex = char_per_in,
      min_w       = min_w,
      max_w       = max_w,
      ladderize_it = ladderize_it
    )
    message(sprintf("Done. Wrote %d page(s).", npages))
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Run
tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(save = "no", status = 1)
})