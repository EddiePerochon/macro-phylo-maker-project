# Macro Phylo Maker Project

This repository contains:

- **`MacroPhyloMaker/`** — an R package for phylogenetic grafting
- **`project/`** — input trees and tables used to generate results

The goal is to allow fully reproducible construction of large phylogenies by:
1. Extracting clades from published trees
2. Preparing donor templates
3. Grafting clades and tips onto a backbone chronogram

This workflow was developed by Marek Borowiec [marek.borowiec@colostate.edu](mailto:marek.borowiec@colostate.edu) and Eddie Pérochon [eddie.perochon@hotmail.com](mailto:eddie.perochon@hotmail.com).

If you use our approach, please cite:
```
Pérochon, E., Bertelsmeier, C., Borowiec, M.L. (2026). Assembling the ant tree of life through scalable phylogenetic synthesis. Journal, volume, pages. doi:XXXX.
```

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
  "here",
  "readr"
))
```
You may need to install additional dependencies for `devtools` and others to run. This will vary a lot by system and environment, but on my Ubuntu 22.04 installing these packages was necessary:
```
sudo apt update

sudo apt install \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev
    libwebp-dev \
    libuv1-dev
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
res <- read.tree(
  here::here(
    "project", "published", "formicidae", "nelsen2018",
    "Dryad_Supplementary_File_7_ML_TREE_treepl_185.tre"
  )
)

extract_clade_with_outgroup(
  res,
  genus = NULL,
  mrca_tips = c("Solenopsis_papuana", "Solenopsis_xyloni"),
  outgroup = "sister_one",
  clean = "none",
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
- note that `clean = "genus_species"` will not work in this mode; grafting functions will still clean 

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

- Clade graft (MRCA mode)
  - Multiple taxa (comma-separated)
  - Donor is grafted at the MRCA of those taxa

- Tip replacement (single label)
  - Single taxon
  - Donor replaces that terminal branch

**How clade grafting works**
For each row:

1. Read donor tree
2. Prepare template:
  - clean tip labels
  - infer ingroup and (optionally) outgroup
  - convert to chronogram if needed
  - compute donor stem fraction r
3. Identify placement in backbone:
  - MRCA of anchors (clade mode), or
  - specific tip (tip mode)
4. Modify backbone:
  - drop existing taxa at the target clade
  - retain one representative tip (internally)
5. Graft donor:
  - insert scaled donor crown along a branch
  - placement depth determined by r or a distribution
6. Repeat for all rows

**Stem vs crown grafting**
Controlled by the optional `Stem_mode` column in the plan:

- outgroup (default)
  - uses donor stem length (more realistic timing)

- crown
  - ignores stem, grafts using crown-only placement

**Time scaling and chronograms**

- If a donor tree is not ultrametric:
  - converted to a chronogram using ape::chronos()
  - best-fitting model selected automatically (chronos_select = "auto")

- If already ultrametric:
  - used directly

**Outputs**
The function writes:
```
<out_prefix>.tre → final grafted tree
<out_prefix>_graft_log.tsv → all graft operations
<out_prefix>.pdf → tree visualization (optional)
<out_prefix>_tips.txt → final tip list
```
---

## Time-informed species grafting with ChronoSTA

After tip grafting and clade grafting, some species may still remain unplaced because they occur in source trees that partially overlap with trees already used for clade grafting. In these cases, choosing a single donor tree for a clade can leave species from alternative, overlapping source trees behind. The function `run_chronosta_grafting()` performs a final time-informed grafting step to recover these species while minimizing disruption to the existing backbone topology.

This step uses ChronoSTA as a controlled gap-filling procedure. The reference tree is the output of the previous MacroPhyloMaker grafting steps. Donor trees are searched for species absent from the reference tree, recalibrated to the reference timescale when needed, split into smaller subtrees, optionally prefused when they overlap, and then merged back into local regions of the reference tree. This allows additional species to be incorporated without rebuilding the entire tree from scratch.

### Note on ChronoSTA compatibility patch

The workflow downloads the upstream `chronosta.py` script from the [Chrono-STA repository](https://github.com/josebarbamontoya/chrono-sta) at runtime. To maintain compatibility with recent NumPy/Pandas versions, MacroPhyloMaker applies a small runtime patch to the temporary copy of `chronosta.py` used for each analysis. The patch replaces an in-place modification of `m.values` with a writable copy before calling `np.fill_diagonal()`.

The original downloaded script is not modified. The patch is applied only to the temporary copy written into the ChronoSTA run directory, and the log records when the patch is applied. This preserves provenance while allowing the workflow to run reproducibly with current Python environments.

Users may alternatively provide their own Chrono-STA script with `chronosta_script = "path/to/chronosta.py"`. If the script already contains the compatibility patch, MacroPhyloMaker will detect this and skip patching.

ChronoSTA is developed independently and distributed by its authors. MacroPhyloMaker does not vendor a modified copy of ChronoSTA. Instead, it downloads or uses a user-provided upstream `chronosta.py` script and applies a documented runtime compatibility patch to the temporary copy used for each analysis. Users should cite ChronoSTA and comply with its license when distributing modified Chrono-STA code.

When using this step of the workflow please cite:
```
Barba-Montoya, J., Craig, J. M., & Kumar, S. (2025). Integrating phylogenies with chronology to assemble the tree of life. Frontiers in Bioinformatics, 5, 1571568.
```

### Example

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-8Jul2026.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 19982018,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc = c(
    "Monomorium",
    "Camponotus",
    "Rogeria",
    "Syllophopsis",
    "Neoponera"
  )
)
```

### Inputs

- `reference_tree`
  - A Newick tree path or a `phylo` object.
  - Usually the output of clade grafting followed by genus reconstitution.
  - This tree should be ultrametric and should contain the species already incorporated by previous grafting steps.

- `donor_tree_dir`
  - Folder containing additional donor trees in `.nwk` or `.tre` format.
  - These trees should already have harmonized tip labels in `Genus_species` format.
  - Donor trees are scanned for species that are absent from the reference tree.

- `out_prefix`
  - Prefix for all output files.
  - Output files, logs, intermediate tables, and final trees are written using this prefix.

- `split_seed`
  - Seed controlling stochastic subtree splitting.
  - Set this value for reproducible runs.

- `ref_weight`
  - Weight assigned to the reference tree during Chrono-STA merging.
  - Internally, this is implemented by duplicating the reference tree when Chrono-STA is run.
  - Larger values make the final merged tree more strongly anchored to the reference topology.

- `paraph_exc`
  - Genera that are allowed to remain non-monophyletic in the final cleanup step.
  - These are cases where non-monophyly is already known or accepted from source phylogenies.

### What the function does

For each donor tree, `run_chronosta_grafting()`:

1. Identifies species present in the donor tree but absent from the reference tree.
2. Identifies shared nodes between the donor tree and reference tree.
3. Optionally recalibrates donor trees to the reference timescale using `chronos()`.
4. Splits donor trees at the genus level to avoid imposing deep donor-tree topology onto the reference.
5. Further splits donor trees into smaller subtrees containing:
   - missing species,
   - representative shared species,
   - and local outgroups where possible.
6. Optionally prefuses overlapping subtrees with ChronoSTA.
7. Merges each subtree locally with the corresponding region of the reference tree.
8. Rescales the ChronoSTA output to the original reference-tree timescale.
9. Grafts the merged subtree back onto the reference.
10. Optionally restores genus-level monophyly by removing newly introduced conflicts while preserving known non-monophyletic genera.

### ChronoSTA setup

The function can automatically download the ChronoSTA Python script if it is not found. Python must be available on the system, and the following Python packages are required:

```bash
python -m pip install biopython pandas numpy scipy matplotlib
```

Because many systems have more than one Python installation, it is safest to choose one Python executable explicitly and use it consistently for setup, dependency checks, and the final grafting call. For example, on a system using Anaconda:

```r
py <- "/home/marek/anaconda3/bin/python3"

setup_chronosta_env(python = py)
check_chronosta_python(py)
```

Then pass the same Python executable to `run_chronosta_grafting()`:

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-8Jul2026.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 19982018,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc = c(
    "Monomorium",
    "Camponotus",
    "Rogeria",
    "Syllophopsis",
    "Neoponera"
  ),
  python = py
)
```

If the dependency check reports missing Python modules even after setup, confirm that R and ChronoSTA are using the same Python:

```r
system2(py, c("-c", "import sys; print(sys.executable)"), stdout = TRUE)
system2(py, c("-m", "pip", "--version"), stdout = TRUE)
check_chronosta_python(py)
```

For Anaconda users, dependencies can also be installed with conda:

```bash
/home/marek/anaconda3/bin/conda install -c conda-forge biopython pandas numpy scipy matplotlib
```

If you already have a local copy of `chronosta.py`, you can provide it explicitly:

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-8Jul2026.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  chronosta_script = here::here(
    "software", "chrono-sta", "code", "chronosta.py"
  )
)
```

### Outputs

The function writes several files:

```text
<out_prefix>_final_tree.nwk
<out_prefix>_final_tree_monophyletic.nwk
<out_prefix>_overlap.tsv
<out_prefix>_splits.tsv
<out_prefix>_prefusion_overlap.tsv
<out_prefix>_prefusion_groups.tsv
<out_prefix>_monophyly_removed.tsv
logs/chronosta_grafting_YYYYMMDD_HHMMSS.log
```

The returned object contains:

```r
res$tree              # final tree used for downstream analyses
res$raw_tree          # tree before optional monophyly cleanup
res$overlap           # donor tree overlap summary
res$split_log         # subtree splitting summary
res$monophyly         # monophyly cleanup details
res$paths             # file paths to outputs
res$parameters        # parameters used in the run
```

### Typical workflow

ChronoSTA grafting is used after the main tip and clade grafting steps:

1. Start with a genus-level backbone.
2. Add missing genera with `run_tip_grafting()`.
3. Add species-level clades with `run_clade_grafting()`.
4. Reconstitute genus tips removed during clade grafting, if needed, using `run_tip_grafting()`.
5. Add remaining species from overlapping donor trees using `run_chronosta_grafting()`.

This final step is intended for species that are present in published source trees but could not be incorporated by direct clade grafting because their source trees partially overlap with, but were not selected over, other donor trees.


### ChronoSTA troubleshooting

#### `setup_chronosta_env()` succeeds, but `run_chronosta_grafting()` says Python packages are missing

This usually means that setup installed packages into one Python environment, while the final grafting call is using a different Python executable. Always set a Python path explicitly and use it in all three places:

```r
py <- "/home/marek/anaconda3/bin/python3"

setup_chronosta_env(python = py)
check_chronosta_python(py)

res <- run_chronosta_grafting(
  ...,
  python = py
)
```

#### `sh: 1: Syntax error: word unexpected (expecting ")")`

This indicates that the Python dependency check is being interpreted by the shell rather than passed safely to Python. Use the patched version of `check_chronosta_python()` that writes a temporary Python script and runs that file, rather than sending an inline Python expression through `python -c`.

#### Avoid unnecessary upgrades in a shared Python environment

ChronoSTA requires `biopython`, `pandas`, `numpy`, `scipy`, and `matplotlib`, but it does not require the newest available versions. If you are using an Anaconda base environment that also supports other analyses, avoid unnecessary `--upgrade` installs unless needed. Prefer:

```bash
/home/marek/anaconda3/bin/python3 -m pip install biopython pandas numpy scipy matplotlib
```

or:

```bash
/home/marek/anaconda3/bin/conda install -c conda-forge biopython pandas numpy scipy matplotlib
```

#### Checking imports manually

A direct check outside MacroPhyloMaker is:

```r
py <- "/home/marek/anaconda3/bin/python3"

system2(
  py,
  c(
    "-c",
    "import Bio, pandas, numpy, scipy, matplotlib; print('ChronoSTA imports OK')"
  ),
  stdout = TRUE,
  stderr = TRUE
)
```

If this command works but `check_chronosta_python(py)` fails, the package helper function should be updated.


---



## Taxonomic completion with TACT

After tip grafting, clade grafting, and optional ChronoSTA gap-filling, the tree may still be incomplete relative to the current species-level taxonomy. `run_tact_grafting()` performs a final taxonomic completion step using TACT. MacroPhyloMaker does not distribute TACT; instead, it prepares TACT-ready inputs, calls an external TACT installation or Docker image, restores temporary labels, removes scaffold taxa, and writes validation reports.

### What this step does

The wrapper:

1. Reads the current backbone tree.
2. Reads the AntWiki taxonomy authority file.
3. Uses `TaxonName` as the authoritative name field for AntWiki input.
4. Keeps only binomial species names, i.e. rows of the form `Genus species`.
5. Omits trinomials, subspecies, and AntWiki rows where the genus parsed from `TaxonName` disagrees with the `Genus` column.
6. Detects non-monophyletic genera and temporarily splits them into TACT-safe pseudo-genera.
7. Optionally assigns species to biogeographic realms using type-locality country and a country-to-realm table.
8. Optionally uses biogeographic realms to assign missing species preferentially to same-realm backbone anchors within each genus.
9. Runs TACT to add missing species and simulate branching times.
10. Restores temporary labels in the final tree.
11. Removes scaffold terminals such as code-like species and placeholder labels, e.g. `Eburopone_CM02`, `Recurvidris_TH01`, and `Uwari_sp`.
12. Writes validation reports comparing the final tree with the processed taxonomy.

### TACT installation

The easiest route is Docker:

```bash
docker run --rm jonchang/tact tact_add_taxa --help
```

If this command succeeds, the wrapper can be run with:

```r
tact_runner = "docker"
docker_image = "jonchang/tact"
```

A local system installation can also be used with:

```r
tact_runner = "system"
tact_bin = "tact_add_taxa"
```

Use `tact_runner = "none"` to prepare TACT input files without running TACT.

### AntWiki taxonomy handling

For `taxonomy_format = "antwiki"`, the wrapper uses the `TaxonName` column rather than constructing names from `Genus` and `Species`. This avoids artifacts caused by AntWiki rows where the `Species` column contains an entire binomial. Only two-token `TaxonName` values are retained:

```text
Genus species
```

Three-token names are treated as trinomials or infraspecific names and omitted. Rows where the genus in `TaxonName` does not match the `Genus` column are treated as probable table errors and omitted with a warning.

### Biogeography-aware TACT completion

If `biogeo = TRUE`, the wrapper assigns species to realms using a country-to-realm table. The table should contain either:

```text
Country    UdvardyRealm
```

or:

```text
Country    Realm
```

For this project, the country table is expected at:

```text
project/tables/country_udvardy_realm.tsv
```

The biogeography-aware workflow uses type-locality country as a proxy for biogeographic affinity. Existing backbone representatives are assigned realms when possible. Missing species are then preferentially assigned to temporary same-realm anchors within their genus. If no same-realm anchor exists, the species is assigned among all available anchors for that genus. This is especially useful for globally distributed genera where missing species should not be grafted indiscriminately across distant biogeographic clades.

### Example: biogeography-aware TACT run

```r
res_tact <- run_tact_grafting(
  backbone_tree = here::here(
    "project", "results", "grafted",
    "final_tree_monophyletic_4183sp.tre"
  ),
  taxonomy = here::here(
    "project", "tables",
    "antwiki-valid-species-8Mar2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "tact",
    "Formicidae_complete_tact_biogeo"
  ),
  taxonomy_format = "antwiki",
  tact_runner = "docker",
  docker_image = "jonchang/tact",
  outgroups = NULL,
  seed = 42,
  nonmono = "split",
  nonmono_allocation = "proportional",
  genus_only = "replace_random_species",
  species_code = "keep",
  drop_code_species_after_tact = TRUE,
  enforce_taxonomy_tip_count = FALSE,
  biogeo = TRUE,
  country_realm_map = here::here(
    "project", "tables",
    "country_udvardy_realm.tsv"
  ),
  biogeo_unknown_realm = "Unknown",
  biogeo_apply_to_all_genera = TRUE,
  write_biogeo_labelled_trees = TRUE,
  keep_temp = TRUE
)
```

### Important arguments

- `backbone_tree`: current tree to complete.
- `taxonomy`: AntWiki or other taxonomy table.
- `out_prefix`: prefix for all output files.
- `tact_runner`: `"docker"`, `"system"`, or `"none"`.
- `seed`: controls wrapper-level stochastic choices, including allocation of missing species to temporary anchors.
- `nonmono`: use `"split"` for normal treatment of non-monophyletic genera.
- `nonmono_allocation`: use `"proportional"`, `"equal"`, or `"random"`.
- `drop_code_species_after_tact`: removes code-like and placeholder scaffold taxa after TACT.
- `enforce_taxonomy_tip_count`: if `FALSE`, writes reports and warns rather than stopping on tree/taxonomy mismatches.
- `biogeo`: enables biogeography-aware allocation.
- `country_realm_map`: country-to-realm table.
- `biogeo_apply_to_all_genera`: if `TRUE`, applies realm-aware temporary anchors to all genera with backbone representatives and missing species.
- `write_biogeo_labelled_trees`: writes inspection trees with realm names appended to tip labels.
- `keep_temp`: keeps the TACT work directory for debugging.

### Main outputs

The main cleaned final tree is:

```text
<out_prefix>_tacted_cleaned.tre
```

Important audit outputs are:

```text
<out_prefix>_validation.tsv
<out_prefix>_taxonomy_mismatch_summary.tsv
<out_prefix>_tree_not_taxonomy.tsv
<out_prefix>_taxonomy_not_tree.tsv
<out_prefix>_dropped_code_species.tsv
<out_prefix>_biogeo_taxonomy_realms.tsv
<out_prefix>_nonmono_temp_name_map.tsv
<out_prefix>_nonmono_genera.tsv
<out_prefix>_skipped_taxa.tsv
<out_prefix>_excluded_taxonomy.tsv
<out_prefix>_excluded_tips.txt
```

If `write_biogeo_labelled_trees = TRUE`, the wrapper also writes:

```text
<out_prefix>_backbone_biogeo_labels.tre
<out_prefix>_tacted_cleaned_biogeo_labels.tre
```

These realm-labelled trees are intended for visual inspection and should not be treated as the primary clean analysis tree.

### Checking the TACT output

```r
ape::Ntip(res_tact$tree)
read.delim(res_tact$paths$validation)
read.delim(res_tact$paths$taxonomy_mismatch_summary)
read.delim(res_tact$paths$biogeo_taxonomy_realms)

grep("TACTTMP|TACTEXCL", res_tact$tree$tip.label, value = TRUE)
```

The last command should return `character(0)` for the cleaned tree.

### Randomized TACT replicates

The wrapper contains stochastic steps, and TACT itself performs stochastic grafting. To generate multiple plausible completed trees, run the TACT step with different seeds and different output prefixes. For stronger wrapper-level randomization, use `nonmono_allocation = "random"`.

```r
seeds <- 1:20

for (seed_i in seeds) {
  run_tact_grafting(
    backbone_tree = here::here(
      "project", "results", "grafted",
      "final_tree_monophyletic_4183sp.tre"
    ),
    taxonomy = here::here(
      "project", "tables",
      "antwiki-valid-species-8Mar2026.txt"
    ),
    out_prefix = here::here(
      "project", "results", "tact",
      sprintf("Formicidae_complete_tact_biogeo_rep%03d", seed_i)
    ),
    taxonomy_format = "antwiki",
    tact_runner = "docker",
    docker_image = "jonchang/tact",
    seed = seed_i,
    nonmono = "split",
    nonmono_allocation = "random",
    genus_only = "replace_random_species",
    species_code = "keep",
    drop_code_species_after_tact = TRUE,
    enforce_taxonomy_tip_count = FALSE,
    biogeo = TRUE,
    country_realm_map = here::here(
      "project", "tables",
      "country_udvardy_realm.tsv"
    ),
    write_biogeo_labelled_trees = FALSE,
    keep_temp = FALSE
  )
}
```

For exploratory visual checks, use `write_biogeo_labelled_trees = TRUE` and `keep_temp = TRUE`. For large replicate batches, set both to `FALSE`.

### TACT troubleshooting

#### Docker cannot run TACT

Confirm Docker works:

```bash
docker run --rm jonchang/tact tact_add_taxa --help
```

#### Taxonomy and tree tip counts do not match

This does not necessarily mean the tree was not created. Inspect:

```text
<out_prefix>_taxonomy_mismatch_summary.tsv
<out_prefix>_tree_not_taxonomy.tsv
<out_prefix>_taxonomy_not_tree.tsv
```

These reports distinguish scaffold labels, omitted table-error rows, and true name mismatches.

#### Many species have `Unknown` realm

This usually means the type-locality string is not an exact match to a row in the country-to-realm map. Add rows for historical or regional strings such as `Tropical Africa`, `South America`, or `East Indies` if those should map to a realm.


## Replicating the ant macrophylogeny workflow from the paper

The full ant-tree workflow used in the paper can be rerun from the files provided in this repository. The goal is to make each major step explicit, reproducible, and updateable when new phylogenies or taxonomy files become available.

### Step 0. Clone repository and load package

```bash
git clone https://github.com/marekborowiec/macro-phylo-maker-project.git
cd macro-phylo-maker-project
```

Then in R:

```r
install.packages(c(
  "devtools",
  "ape",
  "phytools",
  "phangorn",
  "stringr",
  "progress",
  "igraph",
  "MonoPhy",
  "here",
  "readr"
))

devtools::load_all("MacroPhyloMaker")
```

If using the ChronoSTA-enabled final grafting step, also set up Python dependencies. Use an explicit Python executable so that setup, checks, and the final run all use the same environment:

```r
py <- "/home/marek/anaconda3/bin/python3"
setup_chronosta_env(python = py)
check_chronosta_python(py)
```

### Step 1. Graft missing genera onto the genus-level backbone

```r
run_tip_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path     = here::here("project", "tables", "grafted_genera.tsv"),
  out_prefix    = here::here("project", "results", "grafted", "genus"),
  seed_mode     = 42,
  ingroup_anchors = c("Martialis", "Camponotus")
)
```

This creates an updated genus-level backbone by inserting genera missing from the starting chronogram according to the table in `project/tables/grafted_genera.tsv`.

### Step 2. Graft published clade-level phylogenies

```r
run_clade_grafting(
  backbone_path = here::here("project", "backbones", "genus.tre"),
  plan_path = here::here("project", "tables", "clades-to-graft-clean.tsv"),
  authority = here::here(
    "project", "tables",
    "antwiki-valid-species-8Mar2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "backbone_clade_grafted_24Jun2026"
  ),
  seed_mode = 42,
  chronos_select = "auto",
  ultrametric_final = "none",
  plot_pdf = TRUE,
  pdf_auto = TRUE,
  plot_cex = 0.35
)
```

This step reads the clade grafting table, filters donor phylogenies against the authority file, converts phylograms to chronograms when needed, prepares graftable templates, and inserts donor clades into the backbone.

### Step 3. Reconstitute genera removed during clade grafting

Some genus-level terminals may be removed or overwritten when larger clades are grafted wholesale. These can be restored using a second tip-grafting step:

```r
run_tip_grafting(
  backbone_path = here::here(
    "project", "results", "grafted",
    "backbone_clade_grafted_24Jun2026.tre"
  ),
  plan_path = here::here(
    "project", "tables",
    "reconstitute-with-genera-post-grafting.tsv"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-8Jul2026"
  ),
  seed_mode = 42
)
```

This produces the reference tree for the final ChronoSTA-enabled species grafting step.

### Step 4. Add remaining species using ChronoSTA-enabled time-informed grafting

```r
res <- run_chronosta_grafting(
  reference_tree = here::here(
    "project", "results", "grafted",
    "genus-reconstituted-8Jul2026.tre"
  ),
  donor_tree_dir = here::here(
    "project", "chronosta", "source_trees"
  ),
  out_prefix = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled"
  ),
  split_seed = 19982018,
  ref_weight = 5,
  recalibrate = TRUE,
  split_gen = TRUE,
  split_sbt = TRUE,
  prefuse = TRUE,
  monoph_restore = TRUE,
  paraph_exc = c(
    "Monomorium",
    "Camponotus",
    "Rogeria",
    "Syllophopsis",
    "Neoponera"
  ),
  python = py
)
```

This step searches additional donor trees for species absent from the reference tree, calibrates them to the reference timescale where possible, divides them into smaller grafting units, uses ChronoSTA to merge missing species with the relevant reference subtrees, rescales the outputs, and grafts them back into the reference tree.

### Step 5. Complete the tree with biogeography-aware TACT

After ChronoSTA-enabled grafting, use TACT to add remaining species from the AntWiki taxonomy. This step retains valid binomials from `TaxonName`, omits malformed AntWiki rows where the genus fields disagree, optionally assigns species to biogeographic realms using type-locality country, and runs TACT to complete the tree.

```r
res_tact <- run_tact_grafting(
  backbone_tree = here::here(
    "project", "results", "grafted",
    "chronosta_gapfilled_final_tree_monophyletic.nwk"
  ),
  taxonomy = here::here(
    "project", "tables",
    "antwiki-valid-species-8Mar2026.txt"
  ),
  out_prefix = here::here(
    "project", "results", "tact",
    "Formicidae_complete_tact_biogeo"
  ),
  taxonomy_format = "antwiki",
  tact_runner = "docker",
  docker_image = "jonchang/tact",
  outgroups = NULL,
  seed = 42,
  nonmono = "split",
  nonmono_allocation = "proportional",
  genus_only = "replace_random_species",
  species_code = "keep",
  drop_code_species_after_tact = TRUE,
  enforce_taxonomy_tip_count = FALSE,
  biogeo = TRUE,
  country_realm_map = here::here(
    "project", "tables",
    "country_udvardy_realm.tsv"
  ),
  biogeo_unknown_realm = "Unknown",
  biogeo_apply_to_all_genera = TRUE,
  write_biogeo_labelled_trees = TRUE,
  keep_temp = TRUE
)
```

The primary TACT-completed tree is:

```r
res_tact$paths$cleaned_tree
```

The realm-labelled inspection trees are:

```r
res_tact$paths$backbone_biogeo_labels
res_tact$paths$cleaned_tree_biogeo_labels
```

### Step 6. Inspect outputs

The final outputs are written to `project/results/grafted/` and include:

```text
chronosta_gapfilled_final_tree.nwk
chronosta_gapfilled_final_tree_monophyletic.nwk
chronosta_gapfilled_overlap.tsv
chronosta_gapfilled_splits.tsv
chronosta_gapfilled_prefusion_overlap.tsv
chronosta_gapfilled_prefusion_groups.tsv
chronosta_gapfilled_monophyly_removed.tsv
logs/chronosta_grafting_YYYYMMDD_HHMMSS.log
```

The most commonly used final tree is:

```r
res$tree
```

or the corresponding file:

```text
project/results/grafted/chronosta_gapfilled_final_tree_monophyletic.nwk
```

### Re-running the workflow after adding new data

To update the tree when new source phylogenies become available:

1. Add the new source tree to the appropriate folder under `project/published/` or `project/chronosta/source_trees/`.
2. If it should be grafted directly as a clade, add or update its row in:

   ```text
   project/tables/clades-to-graft-clean.tsv
   ```

3. If it should be used only in the ChronoSTA gap-filling step, place it in:

   ```text
   project/chronosta/source_trees/
   ```

4. Re-run the workflow from Step 1 or from the earliest affected step.
5. Compare the new logs, overlap tables, tip counts, and final trees.

For reproducible paper-style runs, do not use `setwd()` or absolute paths. Use `here::here(...)` throughout, and keep seeds fixed for stochastic steps.