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

### Example:

```
tree_path <- here::here(
  "project", "published", "attini",
  "hanisch2022",
  "Pbruchi_MC1_SN_mcmctree_combr1-4_rename.nwk"
)
```

# Workflow

## Extract clades from source trees

### Example:
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
- `genus` → extract all species in a genus or
- `mrca_tips` → extract a clade defined by the MRCA of specified tips

### Defining the ingroup

Two modes are available:

**Genus mode (genus = "Genus")**

- selects all tips matching `Genus_species`
- optionally collapses duplicate species names
- if the genus is non-monophyletic: 

  - `"prune_extras"` (default): keeps only genus members within the MRCA,
  - `"error"`: aborts

**MRCA mode (mrca_tips = c("taxon1","taxon2"))**

- extracts the clade subtended by those anchors

### Outgroup selection

- `outgroup = "sister_one"` (default):
  - finds the sister clade to the ingroup
  - selects one tip with minimum patristic distance to the ingroup
  - appends that tip to the output tree

- `outgroup = "none"`:
  - returns the ingroup only

### Tree preprocessing

By default, the function prepares trees for downstream use:

- `resolve_polytomies = TRUE` → randomly resolves polytomies
- `force_positive_lengths = TRUE` → removes zero/negative branch lengths
- `clean = "genus_species"` → standardizes tip labels and collapses duplicates

Set seed to make these steps reproducible.

### Outputs

The function returns a list:

- `tree` → extracted subtree
- `outgroup` → chosen outgroup tip (if any)
- `ingroup_tips` → tips retained
- `paths` → output file paths

If `write_tree = TRUE`, the subtree is written automatically:
```
<stem>_with_outgroup.tre
<stem>_renamed.tsv
<stem>_dropped.tsv
```

All files are written to out_dir unless explicitly overridden.

### Typical use in this project

This function is used to:

- Extract clades from published phylogenies
- Standardize labels and remove duplicates
- Add a single outgroup tip
- Produce donor trees for clade grafting

These outputs are then referenced in the clade grafting plan file.

## Grafting tips
The function `run_tip_grafting()` inserts tips into a backbone phylogeny according to a predefined grafting table. This can be the first grafting step, producing a backbone used for downstream clade grafting, or used to attach isolated tips missing from grafted clades.

### Example:
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

### Inputs

- Backbone tree (`backbone_path`)
  - A phylogeny (Newick)
  - Must already contain the taxa used as anchors for grafting
  - Typically produced or curated prior to grafting

- Grafting table (`plan_path`)
  - Tab-delimited file specifying how each genus should be inserted
  - Each row corresponds to one grafting operation

**Grafting table structure**

Minimum required columns: `GraftedTip`, `Function`, `Sister`

Example:
```
GraftedTip   Function                  Sister
Poneracantha graft sister to clade     Gnamptogenys,Typhlomyrmex
Alfaria      graft sister to tip       Poneracantha
```

**Supported grafting operations**

- graft sister to tip
  - Inserts a new genus as sister to an existing tip
- graft sister to clade
  - Inserts a genus as sister to the MRCA of two or more taxa
- graft within clade random
  - Inserts a genus somewhere within a specified clade
  - Placement is randomized (controlled by `seed_mode`)

### Branch attachment along a distribution
For grafting operations, the exact point where a new lineage attaches along a branch is drawn from a continuous distribution between two relative positions on that branch.

- The attachment location is expressed as a fraction of branch length:
  - `0` = immediately at the parent node
  - `1` = at the descendant tip
- The interval for possible placement is controlled by:
  - `min_frac` → lower bound
  - `max_frac` → upper bound
- The function samples a position randomly within this interval, using a distribution defined by:
  - `shape1`, `shape2` (parameters of a Beta distribution)

**Interpretation**

- Uniform placement
  - `shape1` = 1, `shape2` = 1
  - attachment equally likely anywhere between `min_frac` and `max_frac`
- Bias toward the base of the branch
  - `shape1` < `shape2`
  - attachment closer to the parent node
- Bias toward the tip
  - `shape1` > `shape2`
  - attachment closer to the descendant lineage

**Example**
```
min_frac = 0.1
max_frac = 0.9
shape1 = 1
shape2 = 1
```
→ Attachment is drawn uniformly between 10% and 90% along the branch.

### How graft placement works
For each row:

1. Identify the placement location:
  - single taxon → tip
  - multiple taxa → MRCA of those taxa
2. Insert the new genus (`GraftedTip`) relative to that location
3. Assign branch lengths according to the specified model (or defaults)
4. Repeat for all rows in the table

### Reproducibility

`seed_mode` ensures consistent placement for stochastic operations (e.g., random grafts, branch placement drawn from distribution)
Always set a seed for reproducible pipelines

### Outputs
The function writes:
```
<out_prefix>.tre → updated backbone tree
<out_prefix>_graft_log.tsv → detailed record of graft operations
optional PDF plots if enabled
```
All outputs are written to the directory specified by `out_prefix`.

### Typical workflow
Tip grafting is the first major step in the assembly pipeline:

- Start with a genus-level backbone tree
- Add missing genera using `run_tip_grafting()`
- Use the resulting tree as input for clade grafting

## Clade grafting
The function `run_clade_grafting()` integrates full donor phylogenies (clades) into a backbone tree. It takes a TSV plan describing donor trees and their placement, converts donors into standardized templates, and grafts them onto the backbone while preserving branch‑length structure.

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

### Inputs

- Backbone tree (`backbone_path`)
  - An ultrametric genus-level tree (typically output of tip grafting)
  - Must include anchor taxa used for placement

- Clade grafting plan (`plan_path`)
  - TSV file specifying donor phylogenies and where to graft them
  - Each row = one clade graft

- Authority file (`authority`)
  - Optional species list used to standardize and filter donor tips

**Grafting plan structure**

Minimum required columns: `MRCA`, `Phylogeny_file_path`

Example:
```
MRCA                          Phylogeny_file_path
Atta,Acromyrmex               project/published/attini/tree.tre
Camponotus                    project/published/camponotus/tree.tre
```
**Placement modes**
Placement is determined automatically from the MRCA column:

Clade graft (MRCA mode)

Multiple taxa (comma-separated)
Donor is grafted at the MRCA of those taxa

Tip replacement (single label)

Single taxon
Donor replaces that terminal branch

**How clade grafting works**
For each row:

Read donor tree
Prepare template:

clean tip labels
infer ingroup and (optionally) outgroup
convert to chronogram if needed
compute donor stem fraction r

Identify placement in backbone:

MRCA of anchors (clade mode), or
specific tip (tip mode)

Modify backbone:

drop existing taxa at the target clade
retain one representative tip (internally)

Graft donor:

insert scaled donor crown along a branch
placement depth determined by r or a distribution

Repeat for all rows

**Stem vs crown grafting**
Controlled by the optional Stem_mode column in the plan:

outgroup (default)

uses donor stem length (more realistic timing)

crown

ignores stem, grafts using crown-only placement

**Time scaling and chronograms**

If a donor tree is not ultrametric:

converted to a chronogram using ape::chronos()
best-fitting model selected automatically (chronos_select = "auto")

If already ultrametric:

used directly

**Branch placement**

Graft position along a branch is determined by:

donor stem fraction (r)
or a Beta-distributed position (crown mode)

Ensures:

consistent scaling between donor and backbone
no zero-length branches
preserved ultrametric structure

**Outputs**
The function writes:
```
<out_prefix>.tre → final grafted tree
<out_prefix>_graft_log.tsv → all graft operations
<out_prefix>.pdf → tree visualization (optional)
<out_prefix>_tips.txt → final tip list
```


