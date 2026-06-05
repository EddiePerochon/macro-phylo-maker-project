#!/usr/bin/env Rscript
# ------------------------------------------------------------
# Fill genus-only tips in a NEWICK tree with random species
# from an AntWiki tab-delimited table (species-rank only).
#
# Usage example:
#   Rscript fill_genus_tips_with_species_cli.R \
#     --ant AntWiki_valid_species_16Oct2025.txt \
#     --tree backbone-tree.txt \
#     --out-tree backbone-tree_species_filled.tre \
#     --out-map label_replacements.csv \
#     --seed 123 \
#     --sample-with-replacement TRUE \
#     --fallback-prefix sp
#
# Notes:
#   - Species-rank = TaxonName has exactly two tokens (Genus species).
#   - Genus-only tips are those with no underscore in tip label.
#   - Ensures unique tip labels; if a duplicate appears, appends _dup1, _dup2...
#   - If genus has no species in table, uses <Genus>_<fallback-prefix>##.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(stringr)
})

# ---------- CLI parsing (base R) ----------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()

if (length(args) > 0) {
  # Accept flags like --key value OR --flag=value
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
          "TRUE" # allow boolean flags like --sample-with-replacement
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

# Help text
if (!is.null(kv[["help"]]) || !is.null(kv[["h"]])) {
  cat("
fill_genus_tips_with_species_cli.R

Required:
  --ant <path>         AntWiki table (tab-delimited)
  --tree <path>        NEWICK tree file

Optional:
  --out-tree <path>    Output NEWICK (default: backbone-tree_species_filled.tre)
  --out-map <path>     Output CSV mapping (default: label_replacements.csv)
  --seed <int>         RNG seed for reproducibility (default: 123)
  --sample-with-replacement <TRUE|FALSE>
                       If a genus occurs more times than species available,
                       sample with replacement (default: TRUE)
  --fallback-prefix <str>
                       Prefix for fallback species names if genus not found
                       in the table (default: sp)

Examples:
  Rscript fill_genus_tips_with_species_cli.R \\
    --ant AntWiki_valid_species_16Oct2025.txt \\
    --tree backbone-tree.txt \\
    --out-tree backbone-tree_species_filled.tre \\
    --out-map label_replacements.csv \\
    --seed 42 \\
    --sample-with-replacement FALSE \\
    --fallback-prefix cf

")
  quit(status = 0)
}

# ---------- Options ----------
ant_file   <- get_opt("ant",  required = TRUE)
tree_file  <- get_opt("tree", required = TRUE)
out_tree   <- get_opt("out-tree", "backbone-tree_species_filled.tre")
out_map    <- get_opt("out-map",  "label_replacements.csv")
seed_val   <- as.integer(get_opt("seed", "123"))
with_repl  <- to_logical(get_opt("sample-with-replacement", "TRUE"), default = TRUE)
fb_prefix  <- get_opt("fallback-prefix", "sp")

set.seed(seed_val)

# ---------- Read AntWiki + species-rank filter ----------
ant <- read_tsv(
  ant_file,
  col_types = cols(.default = col_character())
)

ant_sp <- ant %>%
  mutate(TaxonName = str_squish(TaxonName)) %>%
  filter(!is.na(TaxonName), str_count(TaxonName, "\\S+") == 2) %>%
  select(Genus, Species, TaxonName) %>%
  filter(!is.na(Genus), !is.na(Species), Genus != "", Species != "") %>%
  mutate(Genus = str_trim(Genus),
         Species = str_trim(Species))

if (nrow(ant_sp) == 0) {
  stop("No species-rank rows detected in AntWiki table (TaxonName must be exactly 'Genus species').", call. = FALSE)
}

sp_by_genus <- split(ant_sp$Species, ant_sp$Genus)

# ---------- Read tree ----------
tr <- read.tree(tree_file)
if (is.null(tr$tip.label)) stop("Tree has no tip labels.", call. = FALSE)

labs <- tr$tip.label
is_genus_only <- !str_detect(labs, "_")
genus_only    <- labs[is_genus_only]

message(sprintf("Tree tips: %d | genus-only: %d | genus+species: %d",
                length(labs), length(genus_only), sum(!is_genus_only)))

# ---------- Helpers ----------
make_unique <- function(x, used) {
  if (!(x %in% used)) return(x)
  i <- 1L
  cand <- paste0(x, "_dup", i)
  while (cand %in% used) {
    i <- i + 1L
    cand <- paste0(x, "_dup", i)
  }
  cand
}

# ---------- Build assignments per genus ----------
need_per_genus <- table(genus_only)
assignments <- vector("list", length(need_per_genus))
names(assignments) <- names(need_per_genus)

for (g in names(need_per_genus)) {
  n <- as.integer(need_per_genus[[g]])
  pool <- sp_by_genus[[g]]

  if (is.null(pool) || length(pool) == 0L) {
    species_draw <- sprintf("%s%02d", fb_prefix, seq_len(n))  # e.g., sp01, sp02
  } else {
    # Try without replacement if possible, else depending on flag:
    if (length(pool) >= n) {
      species_draw <- sample(pool, n, replace = FALSE)
    } else {
      if (with_repl) {
        species_draw <- sample(pool, n, replace = TRUE)
      } else {
        # fill with all available unique species, then fallback for overflow
        species_draw <- c(
          sample(pool, length(pool), replace = FALSE),
          sprintf("%s%02d", fb_prefix, seq_len(n - length(pool)))
        )
      }
    }
  }
  assignments[[g]] <- paste0(g, "_", species_draw)
}

# ---------- Apply replacements ----------
used_labels <- labs[!is_genus_only]
new_labs <- labs
map_df <- tibble(old = character(), new = character())

idx_genus_only <- which(is_genus_only)
leftover <- lapply(assignments, identity)

for (i in idx_genus_only) {
  g <- labs[i]
  proposal <- leftover[[g]][1]
  leftover[[g]] <- leftover[[g]][-1]

  final <- make_unique(proposal, used_labels)
  used_labels <- c(used_labels, final)

  new_labs[i] <- final
  map_df <- bind_rows(map_df, tibble(old = g, new = final))
}

tr$tip.label <- new_labs

# ---------- Write outputs ----------
write.tree(tr, file = out_tree)
write_csv(map_df, out_map)

# ---------- Reporting ----------
n_fallback <- sum(str_detect(map_df$new, paste0("_", fb_prefix, "[0-9]+$")))
if (n_fallback > 0) {
  warn_genera <- unique(str_replace(map_df$new[str_detect(map_df$new, paste0("_", fb_prefix, "[0-9]+$"))],
                                    paste0("_", fb_prefix, "[0-9]+$"), ""))
  message(sprintf("NOTE: %d tips used fallback '%s##' because no (or insufficient) species found for these genera: %s",
                  n_fallback, fb_prefix, paste(unique(warn_genera), collapse = ", ")))
}

message(sprintf("Done.\n  NEWICK: %s\n  Mapping: %s", out_tree, out_map))