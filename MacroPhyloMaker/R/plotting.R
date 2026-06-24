#' Suggest a PDF size for plotting a tree
#'
#' Estimates a reasonable PDF size for plotting a phylogeny, using:
#' \itemize{
#'   \item the number of tips (to estimate required height), and
#'   \item the maximum tip-label length (to estimate required width, if enabled).
#' }
#'
#' The heuristic is designed for practical plotting of large trees rather than
#' exact device-space calculation. Height is estimated from an effective
#' tips-per-inch target that scales with \code{cex}; width is optionally
#' estimated from the longest label and an approximate characters-per-inch rule.
#'
#' @param tree A rooted or unrooted \code{phylo} object.
#' @param pdf_width Numeric base width in inches. Used directly when
#'   \code{auto_width = FALSE}; otherwise treated as a fallback value.
#' @param cex Numeric label size passed to plotting.
#' @param target_tpi Target tips per inch used to derive plot height.
#' @param min_h,max_h Minimum and maximum PDF height (in inches).
#' @param buffer_in Extra vertical buffer (in inches) added on top of the
#'   estimated core tree height.
#' @param auto_width Logical; if \code{TRUE}, estimate width from the longest
#'   tip label. If \code{FALSE}, use \code{pdf_width} as supplied.
#' @param char_per_in_at_cex Approximate characters-per-inch at \code{cex = 0.40},
#'   used to estimate width when \code{auto_width = TRUE}.
#' @param min_w,max_w Minimum and maximum PDF width (in inches).
#'
#' @details
#' The size suggestion is heuristic, not device-exact. In particular:
#' \itemize{
#'   \item Height depends only on tip count and \code{cex}, not branch lengths.
#'   \item Width depends only on the longest tip label and \code{cex}, not on
#'         actual rendered string width.
#'   \item This function is intended to provide practical defaults for plotting
#'         large trees, especially in automated workflows.
#' }
#'
#' @return A named list with numeric elements:
#' \itemize{
#'   \item \code{width}: suggested PDF width in inches
#'   \item \code{height}: suggested PDF height in inches
#' }
#'
#' @examples
#' \dontrun{
#' tr <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
#' suggest_pdf_size(tr)
#' }
#'
#' @export
suggest_pdf_size <- function(
    tree,
    pdf_width = 6,
    cex = 0.40,
    target_tpi = 15,
    min_h = 6,
    max_h = 200,
    buffer_in = 2,
    auto_width = TRUE,
    char_per_in_at_cex = 10,
    min_w = 5,
    max_w = 20
) {
  stopifnot(inherits(tree, "phylo"))
  stopifnot(is.numeric(cex), length(cex) == 1L, is.finite(cex), cex > 0)

  ntip <- ape::Ntip(tree)

  # Effective tips-per-inch scales with cex
  tpi_eff <- target_tpi * (0.40 / cex)^0.5
  core_h <- ntip / max(1, tpi_eff)
  h_in <- core_h + buffer_in
  h_in <- max(min_h, min(max_h, h_in))

  if (isTRUE(auto_width)) {
    labs <- as.character(tree$tip.label)
    max_chars <- if (length(labs)) max(nchar(labs), na.rm = TRUE) else 12
    cpi_eff <- char_per_in_at_cex * (0.40 / cex)
    label_w <- max_chars / max(1, cpi_eff)
    pdf_width <- max(min_w, min(max_w, label_w + 2))
  }

  list(width = as.numeric(pdf_width), height = as.numeric(h_in))
}

#' Plot a tree to a single PDF page with optional autosizing
#'
#' Writes a single-page PDF containing a phylogeny, optionally using
#' \code{suggest_pdf_size()} to determine width and/or height.
#'
#' @inheritParams suggest_pdf_size
#' @param pdf_path Output PDF file path.
#' @param pdf_height Numeric PDF height in inches. If \code{NULL}, may be
#'   computed automatically when \code{pdf_auto = TRUE}.
#' @param pdf_auto Logical; if \code{TRUE}, compute width and/or height whenever
#'   either is missing.
#' @param plot_fn Optional function taking a single \code{phylo} argument; if
#'   supplied, this function is used instead of the default
#'   \code{graphics::plot()} call.
#' @param ladderize_it Logical; if \code{TRUE}, ladderize the tree before plotting.
#' @param tight Logical; if \code{TRUE}, use minimal margins for compact plotting.
#' @param mar_safe Numeric margin vector passed to \code{graphics::par(mar = ...)}
#'   when \code{tight = FALSE}.
#' @param oma_safe Numeric outer-margin vector passed to
#'   \code{graphics::par(oma = ...)} when \code{tight = FALSE}.
#' @param ... Additional arguments passed to \code{graphics::plot()} when
#'   \code{plot_fn} is not supplied.
#'
#' @details
#' If \code{pdf_auto = TRUE} and either \code{pdf_width} or \code{pdf_height}
#' is missing, \code{suggest_pdf_size()} is used to fill in the missing
#' dimension(s). If both dimensions are provided, they are used as-is.
#'
#' When \code{plot_fn} is supplied, it is responsible for plotting the tree and
#' receives the (optionally ladderized) \code{phylo} object as its only
#' argument. Additional \code{...} arguments are only used by the default plot path.
#'
#' @return Invisibly returns a named list with:
#' \itemize{
#'   \item \code{width}: final PDF width used
#'   \item \code{height}: final PDF height used
#' }
#'
#' @examples
#' \dontrun{
#' tr <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
#' plot_tree_autosize(tr, tempfile(fileext = ".pdf"))
#' }
#'
#' @export
plot_tree_autosize <- function(
    tree,
    pdf_path,
    cex = 0.40,
    pdf_width = NULL,
    pdf_height = NULL,
    pdf_auto = TRUE,
    target_tpi = 15,
    min_h = 6,
    max_h = 200,
    buffer_in = 2,
    auto_width = TRUE,
    char_per_in_at_cex = 10,
    min_w = 5,
    max_w = 20,
    plot_fn = NULL,
    ladderize_it = TRUE,
    tight = FALSE,
    mar_safe = c(4, 3.5, 3, 1) + 0.1,
    oma_safe = c(0, 0, 0, 0),
    ...
) {
  stopifnot(inherits(tree, "phylo"))
  stopifnot(is.character(pdf_path), length(pdf_path) == 1L, nzchar(pdf_path))
  stopifnot(is.numeric(cex), length(cex) == 1L, is.finite(cex), cex > 0)

  tr <- if (isTRUE(ladderize_it)) ape::ladderize(tree) else tree

  w <- if (is.null(pdf_width)) 6 else as.numeric(pdf_width)
  h <- if (is.null(pdf_height)) 20 else as.numeric(pdf_height)

  if (isTRUE(pdf_auto) && (is.null(pdf_width) || is.null(pdf_height))) {
    dims <- suggest_pdf_size(
      tr,
      pdf_width = w,
      cex = cex,
      target_tpi = target_tpi,
      min_h = min_h,
      max_h = max_h,
      buffer_in = buffer_in,
      auto_width = auto_width,
      char_per_in_at_cex = char_per_in_at_cex,
      min_w = min_w,
      max_w = max_w
    )
    if (is.null(pdf_height)) h <- dims$height
    if (is.null(pdf_width)) w <- dims$width
  }

  grDevices::pdf(pdf_path, width = w, height = h)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  if (isTRUE(tight)) {
    op <- graphics::par(mar = c(1, 1, 1, 1))
    on.exit(graphics::par(op), add = TRUE)

    if (is.function(plot_fn)) {
      plot_fn(tr)
    } else {
      graphics::plot(tr, cex = cex, ...)
    }
  } else {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)

    graphics::par(mar = mar_safe, oma = oma_safe)

    if (is.function(plot_fn)) {
      plot_fn(tr)
    } else {
      graphics::plot(tr, cex = cex, no.margin = FALSE, ...)
    }
  }

  invisible(list(width = w, height = h))
}

#' Plot a tree across multiple PDF pages
#'
#' Splits the tip labels of a tree into consecutive blocks of
#' \code{tips_per_page} and writes one PDF per block.
#'
#' @param tree A \code{phylo} object.
#' @param out_prefix Output prefix; page files are named
#'   \code{<out_prefix>_pageNN.pdf}.
#' @param tips_per_page Positive integer giving the number of tips per page.
#' @param cex Label size passed to plotting.
#' @param pdf_auto Logical; if \code{TRUE}, compute page dimensions for each
#'   subtree using \code{suggest_pdf_size()}.
#' @param target_tpi,min_h,max_h,buffer_in,auto_width,char_per_in_at_cex,min_w,max_w
#'   Parameters forwarded to \code{suggest_pdf_size()}.
#' @param ladderize_it Logical; if \code{TRUE}, ladderize the full tree before
#'   splitting it into pages.
#' @param ... Additional arguments passed to \code{graphics::plot()}.
#'
#' @details
#' Paging is based on the current tip order of the (optionally ladderized) tree,
#' not on clade boundaries. This means paged output is intended for readable
#' display of large trees, not for biologically meaningful partitioning.
#'
#' Each page is generated from a subtree produced by \code{ape::keep.tip()} on
#' the corresponding block of tip labels.
#'
#' @return Invisibly returns the number of pages written.
#'
#' @examples
#' \dontrun{
#' tr <- ape::read.tree(text = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")
#' plot_tree_paged(tr, out_prefix = tempfile("tree"), tips_per_page = 2)
#' }
#'
#' @export
plot_tree_paged <- function(
    tree,
    out_prefix = "tree_paged",
    tips_per_page = 2500,
    cex = 0.35,
    pdf_auto = TRUE,
    target_tpi = 16,
    min_h = 6,
    max_h = 200,
    buffer_in = 1.5,
    auto_width = TRUE,
    char_per_in_at_cex = 10,
    min_w = 5,
    max_w = 20,
    ladderize_it = TRUE,
    ...
) {
  stopifnot(inherits(tree, "phylo"))
  stopifnot(is.character(out_prefix), length(out_prefix) == 1L, nzchar(out_prefix))
  stopifnot(is.numeric(tips_per_page), length(tips_per_page) == 1L,
            is.finite(tips_per_page), tips_per_page >= 1)

  tr <- if (isTRUE(ladderize_it)) ape::ladderize(tree) else tree

  n <- ape::Ntip(tr)
  idx <- seq_len(n)
  pages <- split(idx, ceiling(idx / as.integer(tips_per_page)))

  for (k in seq_along(pages)) {
    block <- pages[[k]]
    keep <- tr$tip.label[block]
    subtr <- ape::keep.tip(tr, keep)

    dims <- if (isTRUE(pdf_auto)) {
      suggest_pdf_size(
        subtr,
        pdf_width = 6,
        cex = cex,
        target_tpi = target_tpi,
        min_h = min_h,
        max_h = max_h,
        buffer_in = buffer_in,
        auto_width = auto_width,
        char_per_in_at_cex = char_per_in_at_cex,
        min_w = min_w,
        max_w = max_w
      )
    } else {
      list(width = 6, height = 20)
    }

    pdf_path <- sprintf("%s_page%02d.pdf", out_prefix, k)
    grDevices::pdf(pdf_path, width = dims$width, height = dims$height)

    op <- graphics::par(mar = c(1, 1, 1, 1))
    tryCatch(
      {
        if (ape::Ntip(subtr) < 2L) {
          graphics::plot.new()
          graphics::text(
            x = 0.5, y = 0.5,
            labels = subtr$tip.label,
            cex = cex
          )
        } else {
          graphics::plot(subtr, cex = cex, ...)
        }
      },
      finally = {
        graphics::par(op)
        grDevices::dev.off()
      }
    )
  }

  invisible(length(pages))
}