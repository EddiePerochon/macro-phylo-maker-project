# Encoding-robust AntWiki readers (no external deps)
# Place this file in R/ (e.g., R/antwiki_readers.R) and run devtools::document()

# ---- Internal helpers ------------------------------------------------------

#' Heuristically detect text encoding from raw bytes (UTF-16 BOM/byte-pattern, UTF-8 vs CP1252)
#' @keywords internal
#' @noRd
.detect_text_encoding <- function(path) {
  rb <- readBin(path, what = "raw", n = file.info(path)$size)
  n <- length(rb)
  if (!n) {
    return("UTF-8")
  }

  # BOM detection first
  if (n >= 2) {
    b2 <- as.integer(rb[1:2])
    if (identical(b2, c(0xFF, 0xFE))) {
      return("UTF-16LE")
    }
    if (identical(b2, c(0xFE, 0xFF))) {
      return("UTF-16BE")
    }
  }
  if (n >= 3) {
    b3 <- as.integer(rb[1:3])
    if (identical(b3, c(0xEF, 0xBB, 0xBF))) {
      return("UTF-8")
    }
  }

  # UTF-16 without BOM? Look for lots of NUL bytes on one parity
  zeros <- which(rb == as.raw(0x00))
  p_zero <- length(zeros) / n
  if (p_zero > 0.20) { # many NULs -> likely UTF-16
    # R indices are 1-based. LE tends to have NUL at even positions for ASCII range.
    even_zeros <- sum(zeros %% 2 == 0)
    odd_zeros <- sum(zeros %% 2 == 1)
    if (even_zeros >= odd_zeros) {
      return("UTF-16LE")
    } else {
      return("UTF-16BE")
    }
  }

  # Try reading as UTF-8; if invalid or replacement chars are observed, fallback to CP1252
  read_with <- function(enc) {
    con <- file(path, open = "r", encoding = enc)
    on.exit(close(con), add = TRUE)
    had_warn <- FALSE
    warn_msg <- NULL
    txt <- withCallingHandlers(
      readLines(con, warn = FALSE),
      warning = function(w) {
        had_warn <<- TRUE
        warn_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    list(txt = txt, had_warn = had_warn, warn_msg = warn_msg)
  }

  res <- read_with("UTF-8")
  has_repl <- length(res$txt) && any(grepl("\\ufffd", res$txt))
  if (isTRUE(res$had_warn) || isTRUE(has_repl) || !length(res$txt)) {
    return("windows-1252")
  }
  "UTF-8"
}

#' Clean a header/name: strip BOM/replacement/NBSP and trim
#' @keywords internal
#' @noRd
.clean_name <- function(x) {
  x <- sub("^\\ufeff", "", x) # BOM
  x <- sub("^\\ufffd+", "", x) # leading replacement chars
  x <- gsub("\\u00A0", " ", x, useBytes = FALSE) # NBSP
  trimws(x)
}

# ---- Public API ------------------------------------------------------------

#' Read a delimited text file after normalizing to UTF-8 (robust, no guessing packages)
#'
#' @title Read a delimited text file (normalized to UTF-8)
#' @description
#' Deterministic+heuristic decoding:
#'
#' * Honors **UTF-16 BOM** (LE/BE) if present.
#' * Detects **UTF-16 without BOM** via the NUL-byte parity pattern.
#' * Otherwise tries **UTF-8 first**, and falls back to **Windows-1252 (CP1252)** when
#'   invalid input or replacement characters are encountered.
#'
#' After decoding, the function strips BOM, normalizes non-breaking spaces, and parses
#' with [utils::read.table()]. All character columns are returned in UTF-8.
#'
#' @param path Character path to the file.
#' @param sep Field separator (default "\t" for TSV).
#' @param header Logical; file has a header line? Default `TRUE`.
#' @param quote Quoting characters. For TSVs from spreadsheets, `""` is often safest.
#' @param dec Decimal mark. Default ".".
#' @param fix_html_entities Logical; if `TRUE`, replaces "&" entities (e.g., "&amp;") with
#'   plain `&` in character columns.
#' @param ... Passed to [utils::read.table()].
#' @return A `data.frame` whose character columns are encoded in UTF-8.
#' @md
#' @export
read_table_utf8 <- function(path,
                            sep = "\t",
                            header = TRUE,
                            quote = "",
                            dec = ".",
                            fix_html_entities = TRUE,
                            ...) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)

  enc <- .detect_text_encoding(path)

  # Read entire file as text using detected encoding
  con <- file(path, open = "r", encoding = enc)
  on.exit(close(con), add = TRUE)
  txt <- suppressWarnings(readLines(con, warn = FALSE))

  # If empty and we guessed UTF-8, try CP1252 as last resort
  if (!length(txt) && identical(enc, "UTF-8")) {
    con2 <- file(path, open = "r", encoding = "windows-1252")
    on.exit(close(con2), add = TRUE)
    txt <- suppressWarnings(readLines(con2, warn = FALSE))
  }

  if (!length(txt)) stop("No readable text lines found in file: ", path, call. = FALSE)

  # Strip an optional BOM in first line; normalize to UTF-8 and NBSP -> space
  txt[1] <- sub("^\\ufeff", "", txt[1])
  txt <- enc2utf8(txt)
  txt <- gsub("\\u00A0", " ", txt, useBytes = FALSE)

  # Parse via in-memory text connection
  tc <- textConnection(txt)
  on.exit(close(tc), add = TRUE)
  df <- utils::read.table(tc,
    sep = sep,
    header = header,
    quote = quote,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    dec = dec,
    comment.char = "",
    fill = TRUE,
    ...
  )

  # Clean column names
  names(df) <- vapply(names(df), .clean_name, character(1))

  # Ensure character columns are UTF-8 and optionally decode & entities
  df[] <- lapply(df, function(x) {
    if (is.character(x)) {
      x <- enc2utf8(x)
      if (isTRUE(fix_html_entities)) x <- gsub("&amp;", "&", x, fixed = TRUE)
    }
    x
  })

  df
}

#' Read an AntWiki valid species table (always normalized to UTF-8)
#'
#' @title Read an AntWiki "Valid Species" table (UTF-8 normalized)
#' @description
#' Reads the AntWiki tab-delimited table, converts it to UTF-8 using the robust
#' reader [AntPhyloMaker::read_table_utf8()], then adds a convenience column
#' `binom = paste(Genus, Species, sep = "_")` after trimming whitespace,
#' lowercasing the species epithet, and dropping empty genus/species rows.
#' Duplicate `binom` entries are removed.
#'
#' @param path Character path to the AntWiki plain-text table (tab-delimited).
#'   The table must contain columns `TaxonName`, `Genus`, and `Species`.
#' @return A `data.frame` with at least columns `TaxonName`, `Genus`, `Species`, and `binom`,
#'   with all character data in UTF-8.
#' @seealso [AntPhyloMaker::read_table_utf8()] for the encoding-normalizing reader.
#' @md
#' @export
read_antwiki_table <- function(path) {
  df <- read_table_utf8(path, sep = "\t", header = TRUE, quote = "")

  # Validate required columns with a lenient cleanup pass if needed
  required_cols <- c("TaxonName", "Genus", "Species")
  miss <- setdiff(required_cols, names(df))
  if (length(miss)) {
    names(df) <- gsub("^[^A-Za-z0-9_]+", "", names(df))
    names(df) <- gsub("\r$", "", names(df)) # stray CR at end, just in case
    miss <- setdiff(required_cols, names(df))
  }
  if (length(miss)) {
    stop("AntWiki table must include columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  df$Genus <- trimws(df$Genus)
  df$Species <- tolower(trimws(df$Species))

  keep <- nzchar(df$Genus) & nzchar(df$Species)
  if (any(!keep)) {
    warning(sprintf("Dropped %d row(s) with empty Genus/Species.", sum(!keep)), call. = FALSE)
  }
  df <- df[keep, , drop = FALSE]

  df$binom <- paste(df$Genus, df$Species, sep = "_")

  if (any(duplicated(df$binom))) {
    ndup <- sum(duplicated(df$binom))
    warning(sprintf("Removed %d duplicate binom row(s).", ndup), call. = FALSE)
    df <- df[!duplicated(df$binom), , drop = FALSE]
  }

  df
}

#' List valid species for a genus from an AntWiki table
#'
#' @title List valid species by genus from an AntWiki table
#' @param ant_df A `data.frame` as returned by [AntPhyloMaker::read_antwiki_table()].
#' @param genus Focal genus (case-sensitive).
#' @return Character vector of `"Genus_species"` for that genus.
#' @examples
#' \dontrun{
#' ant <- read_antwiki_table("AntWiki_Valid_Species_23Feb2026.txt")
#' species_for_genus_from_antwiki(ant, "Camponotus")[1:10]
#' }
#' @md
#' @export
species_for_genus_from_antwiki <- function(ant_df, genus) {
  if (!all(c("Genus", "binom") %in% names(ant_df))) {
    stop("ant_df must include columns 'Genus' and 'binom' (see read_antwiki_table()).", call. = FALSE)
  }
  unique(ant_df$binom[ant_df$Genus == genus])
}

# ---- Tiny diagnostics you can keep or remove -------------------------------

#' Print first bytes as hex and the parsed column names (for debugging)
#' @keywords internal
#' @noRd
antwiki_debug_header <- function(path) {
  rb <- readBin(path, what = "raw", n = min(256, file.info(path)$size))
  cat("First 64 bytes (hex):\n",
    paste(sprintf("%02X", as.integer(rb[seq_len(min(64, length(rb)))])), collapse = " "),
    "\n\n",
    sep = ""
  )
  enc <- .detect_text_encoding(path)
  cat("Detected encoding:", enc, "\n")
  df <- try(read_table_utf8(path), silent = TRUE)
  if (inherits(df, "try-error")) {
    cat("read_table_utf8 failed:\n", as.character(df), "\n")
  } else {
    cat("Column names:\n")
    print(names(df))
  }
  invisible(enc)
}
