# Macro Phylo Maker Project

This repository contains:

- **`MacroPhyloMaker/`** — an R package for phylogenetic grafting
- **`project/`** — input data, tables, and scripts used to generate results

The goal is to allow fully reproducible construction of large phylogenies by:
1. Extracting clades from published trees
2. Preparing donor templates
3. Grafting clades and tips onto a backbone chronogram

---

# Setup Instructions

## Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

## Open the R project
**From R interface**, open the file:
```
macro-phylo-maker-project.Rproj
```

## Install required R packages
```
install.packages(c(
  "devtools",
  "ape",
  "phytools",
  "here"
))
```

## Load the development package
```
devtools::load_all("MacroPhyloMaker")
```

# Repository Structure
```
MacroPhyloMaker/        # R package (core functions)
project/
  published/            # input phylogenies
  tables/               # TSV metadata and taxonomy
  scripts/              # analysis scripts
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
All files/scripts must use `here::here(...)`, not relative paths. **This is something that I had to adjust to make this portable, but I am still working to fix it throughout by converting from relative paths.**

Do not use:

absolute paths (/home/...)
setwd()
"../"

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
tree <- read.tree(
  here::here("project", "published", "attini", "hanisch2022",
             "Pbruchi_MC1_SN_mcmctree_combr1-4_rename.nwk")
)

res <- extract_clade_with_outgroup(
  tree,
  mrca_tips = c("Mycocepurus_goeldii_1495", "Pseudoatta_sp_nov_1478"),
  outgroup = "sister_one",
  resolve_polytomies = TRUE,
  force_positive_lengths = TRUE,
  out_dir = here::here("project", "published", "attini", "hanisch2022")
)
```

This:

* extracts an ingroup defined by MRCA
* optionally adds the closest sister outgroup
* writes Newick + log files

## Grafting tips

Example:
```
run_tip_grafting(
  backbone_path = here::here("project", "backbone", "genus.tre"),
  plan_path     = here::here("project", "tables", "grafted_genera.tsv"),
  out_prefix    = here::here("project", "results", "grafted", "genus"),
  seed_mode     = 42,
  ingroup_anchors = c("Martialis", "Camponotus")
)
```
This applies a TSV-defined grafting plan:

* attaches new tips or clades
* writes resulting tree + logs

## Clade grafting

Example:
```
run_clade_grafting(
  backbone_path = here::here("project", "results", "grafted", "genus.tre"),
  plan_path     = here::here("project", "tables", "clades-to-graft.tsv"),
  authority     = here::here("project", "tables", "antwiki-valid-species-8Mar2026.txt"),
  out_prefix    = here::here("project", "results", "grafted", "backbone_clade_grafted"),
  seed_mode     = 42
)
```

This:

* reads donor phylogenies
* prepares templates (chronograms if needed)
* grafts onto backbone (tip or MRCA placement)
* ensures ultrametricity of final tree

## Outputs

All outputs are written to:

project/results/

Including:

`.tre` — Newick trees

`.pdf` — visualizations

`_graft_log.tsv` — detailed logs

`_tips.txt` — final tip lists

