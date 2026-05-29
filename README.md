# Macro Phylo Maker Project

This repository contains:

- **`MacroPhyloMaker/`** — an R package for phylogenetic grafting
- **`project/`** — input trees and tables used to generate results

The goal is to allow fully reproducible construction of large phylogenies by:
1. Extracting clades from published trees
2. Preparing donor templates
3. Grafting clades and tips onto a backbone chronogram

---

# Setup Instructions

## Clone the repository

```bash
git clone https://github.com/marekborowiec/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

## Install required R packages
Inside directory `macro-phylo-maker-project` open `R` command line interface and install:
```
install.packages(c(
  "devtools",
  "ape",
  "phytools",
  "here"
))
```

## Load the development package
From the same directory as above.
```
devtools::load_all("MacroPhyloMaker")
```
This will activate all package functions. Now you should be able to run code examples below. I tested clade grafting to make sure it works. Your current mileage for other functions may vary for the moment.

Some explainers first.

# Repository Structure
```
MacroPhyloMaker/        # R package (core functions)
project/
  backbones/                                         # backbone trees
    Borowiec-129to158Ma-independent-genus-only.tre   # Borowiec et al. 2025 tree
    genus.tree                                       # as above but with missing genera grafted
  published/            # input phylogenies
  tables/               # graft tables and taxonomy authority
  results/              # outputs (trees, logs, PDFs)
```

# Internal Functions
These are defined in the package and used throughout pipelines:

`extract_clade_with_outgroup()` — extract focal clades with optional sister outgroup
`run_tip_grafting()` — apply tip-level grafting plans
`run_clade_grafting()` — full clade-grafting pipeline from TSV plan
`graft_many_clades()` — batch grafting engine
`prepare_clade_template()` — converts donor phylogenies into graftable templates

# Important
All files/scripts must use `here::here(...)`, not relative paths. **This is something that I had to adjust to make this portable, but I am still working to fix it throughout by converting from relative paths.** Clade grafting example below, at least, should work now.

Do not use:

absolute paths `(/home/...)`, `setwd()`, `"../"`

Example:

```
tree_path <- here::here(
  "project", "published", "attini",
  "hanisch2022",
  "Pbruchi_MC1_SN_mcmctree_combr1-4_rename.nwk"
)
```

# Workflow

## Extract clades from source trees

Example:
```
solenopsis <- read.tree(
  here::here(
    "project", "published", "formicidae", "nelsen2018",
    "Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre"
  )
)

extract_clade_with_outgroup(
  solenopsis,
  genus = NULL,
  mrca_tips = c("Solenopsis_papuana", "Solenopsis_xyloni"),
  outgroup = "sister_one",
  clean = "genus_species",
  nonmono = "prune_extras",
  resolve_polytomies = TRUE,
  force_positive_lengths = TRUE,
  seed = 42L,
  write_tree = TRUE,
  tree_path = NULL,
  write_renames = TRUE,
  renames_path = NULL,
  write_drops = TRUE,
  drops_path = NULL,

  out_dir = here::here(
    "project", "published", "formicidae", "nelsen2018"
  )
)
```

This:

* extracts an ingroup defined by MRCA
* optionally adds the closest sister outgroup
* writes Newick + log files

The function `extract_clade_with_outgroup()` is used to extract a focal clade from a larger phylogeny and optionally append a single outgroup taxon for downstream grafting.

### Basic usage

You must provide:

- a tree (phylo object), and exactly one of:
- genus → extract all species in a genus or
- mrca_tips → extract a clade defined by the MRCA of specified tips

### Defining the ingroup

Two modes are available:

-**Genus mode (genus = "Genus")**

- selects all tips matching Genus_species
- optionally collapses duplicate species names
- if the genus is non-monophyletic: 

  - "prune_extras" (default): keeps only genus members within the MRCA,
  - "error": aborts

-**MRCA mode (mrca_tips = c("taxon1","taxon2"))**

- extracts the clade subtended by those anchors

### Outgroup selection

- outgroup = "sister_one" (default):

  - finds the sister clade to the ingroup
  - selects one tip with minimum patristic distance to the ingroup
  - appends that tip to the output tree

- outgroup = "none":

 - returns the ingroup only

## Grafting tips

Example:
```
run_tip_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path     = here::here("project", "tables", "grafted_genera.tsv"),
  out_prefix    = here::here("project", "results", "grafted", "genus"),
  seed_mode     = 42,
  ingroup_anchors = c("Martialis", "Camponotus")
)
```
This applies a TSV-defined grafting plan:

* attaches new tips (in this case, genera missing from Borowiec et al. 2025)
* writes resulting tree + logs

## Clade grafting

Example:
```
run_clade_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path = here::here("project", "tables", "clades-to-graft-clean.tsv"),
  authority = here::here("project", "tables", "antwiki-valid-species-8Mar2026.txt"),
  out_prefix = here::here("project", "results", "grafted",
                          "backbone_clade_grafted_new_bby_test"),
  seed_mode = 42,
  chronos_select = "auto",
  ultrametric_final = "none",
  plot_pdf = TRUE,
  pdf_auto = TRUE,
  plot_cex = 0.35
)
```

This:

* reads donor phylogenies
* prepares templates (chronograms if needed)
* grafts onto backbone (tip or MRCA placement)
* ensures ultrametricity of final tree

## Outputs

All outputs are written to:

`project/results/`

Including:

`.tre` — Newick trees

`.pdf` — visualizations

`_graft_log.tsv` — detailed logs

`_tips.txt` — final tip lists

