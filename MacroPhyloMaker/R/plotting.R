#' Suggest a PDF size for plotting a large tree
#'
#' @param tree A `phylo`.
#' @param pdf_width Base width in inches (or `NULL` to auto).
#' @param cex Label size used for plotting.
#' @param target_tpi Target tips-per-inch (scaled with `cex`).
#' @return A list with `width` and `height` in inches.
#' @export
suggest_pdf_size <- function(
  tree,
  pdf_width = 6,
  cex = 0.40,
  target_tpi = 15, # tips per inch
  min_h = 6,
  max_h = 200,
  buffer_in = 2,
  auto_width = TRUE,
  char_per_in_at_cex = 10,
  min_w = 5,
  max_w = 20
) {
  stopifnot(inherits(tree, "phylo"))
  ntip <- ape::Ntip(tree)
  # Effective tips per inch scales with cex (pipeline's original rule)
  tpi_eff <- target_tpi * (0.40 / cex)^0.5
  core_h <- ntip / max(1, tpi_eff)
  h_in <- core_h + buffer_in
  h_in <- max(min_h, min(max_h, h_in))

  if (isTRUE(auto_width)) {
    labs <- as.character(tree$tip.label)
    max_chars <- if (length(labs)) max(nchar(labs), na.rm = TRUE) else 12
    cpi_eff <- char_per_in_at_cex * (0.40 / cex) # smaller cex => more chars/in
    label_w <- max_chars / max(1, cpi_eff)
    pdf_width <- max(min_w, min(max_w, label_w + 2))
  }
  list(width = pdf_width, height = h_in)
}

#' Plot a tree on a single PDF page with optional autosizing
#' @inheritParams suggest_pdf_size
#' @param pdf_path Output PDF path.
#' @param pdf_auto If `TRUE`, auto-compute width/height unless both supplied.
#' @param plot_fn Optional function(tr) to plot custom content.
#' @param ladderize_it Ladderize prior to plotting.
#' @param tight If `TRUE`, minimal margins (pipeline-like).
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
  tr <- if (isTRUE(ladderize_it)) ape::ladderize(tree) else tree

  # Default size if nothing supplied and pdf_auto = FALSE
  w <- if (is.null(pdf_width)) 6 else as.numeric(pdf_width)
  h <- if (is.null(pdf_height)) 20 else as.numeric(pdf_height)

  # Autosize
  if (isTRUE(pdf_auto) && (is.null(pdf_width) || is.null(pdf_height))) {
    dims <- suggest_pdf_size(
      tr,
      pdf_width = w, cex = cex,
      target_tpi = target_tpi, min_h = min_h, max_h = max_h,
      buffer_in = buffer_in, auto_width = auto_width,
      char_per_in_at_cex = char_per_in_at_cex, min_w = min_w, max_w = max_w
    )
    if (is.null(pdf_height)) h <- dims$height
    if (is.null(pdf_width)) w <- dims$width
  }

  grDevices::pdf(pdf_path, width = w, height = h)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

  if (isTRUE(tight)) {
    op <- graphics::par(mar = c(1, 1, 1, 1))
    on.exit(graphics::par(op), add = TRUE)
    if (is.function(plot_fn)) plot_fn(tr) else graphics::plot(tr, cex = cex, ...)
  } else {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::par(mar = mar_safe, oma = oma_safe)
    if (is.function(plot_fn)) plot_fn(tr) else graphics::plot(tr, cex = cex, no.margin = FALSE, ...)
  }
  invisible(list(width = w, height = h))
}

#' Plot a tree across multiple PDF pages
#' @param out_prefix Output prefix (each page becomes `<prefix>_pageNN.pdf`).
#' @param tips_per_page Number of tips per page.
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
  tr <- if (isTRUE(ladderize_it)) ape::ladderize(tree) else tree

  n <- ape::Ntip(tr)
  idx <- seq_len(n)
  pages <- split(idx, ceiling(idx / tips_per_page))

  for (k in seq_along(pages)) {
    block <- pages[[k]]
    keep <- tr$tip.label[block]
    subtr <- ape::keep.tip(tr, keep)

    dims <- if (pdf_auto) {
      suggest_pdf_size(
        subtr,
        pdf_width = 6, cex = cex, target_tpi = target_tpi,
        min_h = min_h, max_h = max_h, buffer_in = buffer_in,
        auto_width = auto_width, char_per_in_at_cex = char_per_in_at_cex,
        min_w = min_w, max_w = max_w
      )
    } else {
      list(width = 6, height = 20)
    }

    pdf_path <- sprintf("%s_page%02d.pdf", out_prefix, k)
    grDevices::pdf(pdf_path, width = dims$width, height = dims$height)
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    op <- graphics::par(mar = c(1, 1, 1, 1))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot(subtr, cex = cex, ...)
    grDevices::dev.off()
  }
  invisible(length(pages))
}
