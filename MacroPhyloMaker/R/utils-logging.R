#' Lightweight logging utilities (internal)
#'
#' A very small logger used across AntPhyloMaker to write a single
#' session header and then mirror messages to **both** the console
#' and a per-run log file.
#'
#' - Call `.make_logger()` once at the start of a task to open a log file.
#' - Use `log_msg()` for single-line messages.
#' - Use `log_section()` to print a titled box section.
#' - Always call `.close_logger()` (typically via `on.exit()`).
#'
#' The logger writes UTF-8 box characters by default. Set the option
#' `options(AntPhyloMaker.ascii = TRUE)` to force ASCII characters.
#'
#' @keywords internal
#' @noRd

# ---- box-drawing character set (with ASCII fallback) -------------------------
.box_chars <- function(ascii = getOption("AntPhyloMaker.ascii", FALSE)) {
  if (isTRUE(ascii)) {
    list(h = "-", tl = "+", tr = "+", bl = "+", br = "+", v = "|")
  } else {
    list(h = "─", tl = "┌", tr = "┐", bl = "└", br = "┘", v = "│")
  }
}

# ---- create/open a logger ----------------------------------------------------
#' @keywords internal
#' @noRd
.make_logger <- function(outdir,
                         genus,
                         enable = TRUE,
                         mode = c("generic", "extract", "graft"),
                         file_prefix = NULL) {
  mode <- match.arg(mode)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d-%H%M%S")

  # Choose filename stem
  stem <- if (!is.null(file_prefix) && nzchar(file_prefix)) {
    file_prefix
  } else {
    switch(mode,
      extract = sprintf("%s_extract", genus),
      graft   = sprintf("%s_graft", genus),
      generic = sprintf("%s_log", genus)
    )
  }
  logfile <- file.path(outdir, sprintf("%s_%s.log", stem, ts))
  con <- file(logfile, open = "a", encoding = "UTF-8")

  # Header once
  title <- switch(mode,
    extract = "Genus extraction",
    graft   = "Genus grafting",
    generic = "AntPhyloMaker"
  )
  header <- sprintf(
    "=== %s log: %s | %s ===",
    title, genus, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  writeLines(header, con)
  cat(header, "\n")

  structure(
    list(
      enabled = isTRUE(enable),
      file    = logfile,
      con     = con,
      mode    = mode,
      name    = genus,
      header  = header
    ),
    class = "smart_logger"
  )
}

# ---- close a logger (safe) ---------------------------------------------------
#' @keywords internal
#' @noRd
.close_logger <- function(logger) {
  if (inherits(logger, "smart_logger")) {
    try(close(logger$con), silent = TRUE)
  }
  invisible(NULL)
}

# ---- write a message to console and log file ---------------------------------
#' @keywords internal
#' @noRd
log_msg <- function(logger, ..., .timestamp = FALSE) {
  line <- paste0(..., collapse = "")
  if (isTRUE(.timestamp)) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), line)
  }

  # Console
  cat(line, "\n")

  # File (if enabled)
  if (inherits(logger, "smart_logger") && isTRUE(logger$enabled)) {
    writeLines(line, logger$con)
    flush(logger$con)
  }
  invisible(NULL)
}

# ---- draw a titled section box -----------------------------------------------
#' @keywords internal
#' @noRd
log_section <- function(logger, title) {
  chars <- .box_chars()
  bar <- paste(rep(chars$h, max(10, nchar(title) + 2L)), collapse = "")
  log_msg(logger, "")
  log_msg(logger, chars$tl, bar, chars$tr)
  log_msg(logger, chars$v, " ", title, " ", chars$v)
  log_msg(logger, chars$bl, bar, chars$br)
  invisible(NULL)
}
