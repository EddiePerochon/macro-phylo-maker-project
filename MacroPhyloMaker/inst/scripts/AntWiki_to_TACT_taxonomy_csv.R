#!/usr/bin/env Rscript
# ------------------------------------------------------------
# Convert AntWiki valid species table (tab-delimited) into a
# TACT-style taxonomy CSV with columns:
#   Family, genus, genus.species
#
# Species-rank only = TaxonName has exactly two tokens: "Genus species".
# Deduplicates by genus.species, sorts by genus then genus.species.
#
# Usage:
#   Rscript AntWiki_to_TACT_taxonomy_csv_cli.R \
#     --ant AntWiki_valid_species_16Oct2025.txt \ # 14532 valid species
#     --out Formicidae.csv \
#     --family Formicidae \
#     --lowercase-family FALSE \
#     --preview 10
# 
# Data from https://www.antwiki.org/wiki/Downloadable_Data
# https://www.antwiki.org/w/images/9/9e/AntWiki_Valid_Species.txt
#
# ------------------------------------------------------------

# ---------- simple CLI parsing (base R) ----------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()

if (length(args) > 0) {
  i <- 1L
  while (i <= length(args)) {
    token <- args[i]
    if (startsWith(token, "--")) {
      token <- substring(token, 3L)
      if (grepl("=", token, fixed = TRUE)) {
        parts <- strsplit(token, "=", fixed = TRUE)[[1]]
        kv[[parts[1]]] <- parts[2]
        i <- i + 1L
      } else {
        key <- token
        val <- if (i + 1L <= length(args) && !startsWith(args[i + 1L], "--")) {
          args[i + 1L]
        } else {
          "TRUE"  # allow boolean flags like --lowercase-family
        }
        kv[[key]] <- val
        i <- i + ifelse(val == "TRUE" && (i + 1L > length(args) || startsWith(args[i + 1L], "--")), 1L, 2L)
      }
    } else {
      i <- i + 1L
    }
  }
}

get_opt <- function(name, default = NULL, required = FALSE) {
  if (!is.null(kv[[name]])) return(kv[[name]])
  if (required) {
    stop(sprintf("Missing required option --%s. Use --help for usage.", name), call. = FALSE)
  }
  default
}

to_logical <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  x <- toupper(as.character(x))
  x %in% c("1", "TRUE", "T", "YES", "Y")
}

# ---------- help ----------
if (!is.null(kv[["help"]]) || !is.null(kv[["h"]])) {
  cat("
AntWiki_to_TACT_taxonomy_csv_cli.R

Required:
  --ant <path>          AntWiki table (tab-delimited .txt)
  --out <path>          Output CSV file

Optional:
  --family <name>       Family literal for first column (default: Formicidae)
  --lowercase-family <TRUE|FALSE>
                        Lowercase the family literal (default: FALSE)
  --preview <n>         Print first n rows to console (default: 10)

Example:
  Rscript AntWiki_to_TACT_taxonomy_csv_cli.R \\
    --ant AntWiki_valid_species_16Oct2025.txt \\
    --out Formicidae.csv \\
    --family Formicidae \\
    --lowercase-family FALSE \\
    --preview 10

")
  quit(status = 0)
}

# ---------- options ----------
ant_file        <- get_opt("ant", required = TRUE)
out_file        <- get_opt("out", required = TRUE)
family_literal  <- get_opt("family", "Formicidae")
lowercase_fam   <- to_logical(get_opt("lowercase-family", "FALSE"), default = FALSE)
preview_n       <- as.integer(get_opt("preview", "10"))

if (lowercase_fam) family_literal <- tolower(family_literal)

# ---------- read AntWiki table (tab-delimited) ----------
ant <- tryCatch(
  read.table(
    ant_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE
  ),
  error = function(e) {
    stop(sprintf("Failed to read AntWiki table: %s", e$message), call. = FALSE)
  }
)

required_cols <- c("TaxonName", "Genus", "Species")
missing_cols <- setdiff(required_cols, names(ant))
if (length(missing_cols) > 0) {
  stop(sprintf("Input file is missing required columns: %s",
               paste(missing_cols, collapse = ", ")), call. = FALSE)
}

# ---------- species-rank only (TaxonName == 'Genus species') ----------
is_species_rank <- sapply(strsplit(trimws(ant$TaxonName), "\\s+"), length) == 2

sp <- ant[is_species_rank, c("Genus", "Species")]
sp <- sp[nchar(sp$Genus) > 0 & nchar(sp$Species) > 0, , drop = FALSE]

# ---------- build output ----------
# NOTE: header stays 'genus.species' but the cell values use a SPACE.
out <- data.frame(
  Family = family_literal,
  genus = sp$Genus,
  `genus.species` = paste(sp$Genus, sp$Species),
  stringsAsFactors = FALSE
)

# deduplicate & sort
out <- out[!duplicated(out$`genus.species`), ]
ord <- order(out$genus, out$`genus.species`)
out <- out[ord, , drop = FALSE]

# ---------- write CSV without quotes ----------
# write.csv() would add quotes; use write.table() with quote = FALSE
write.table(
  out,
  file = out_file,
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  fileEncoding = "UTF-8"
)

# ---------- report ----------
cat(sprintf("Parsed %d valid %s species (species-rank only; subspecies excluded).\n",
            nrow(out), family_literal))

# ---------- preview ----------
if (!is.na(preview_n) && preview_n > 0) {
  print(utils::head(out, preview_n))
}
