#!/usr/bin/env python3
import re
import sys
from pathlib import Path

import pandas as pd

# ---- User-configurable file names ----
MAP_PATH  = Path("pheidole-table.tsv")  # tab-separated text file
TREE_PATH = Path("BEAST_mcc_tree.tre")         # NEWICK tree file
OUT_TREE  = Path("pheidole-processed-chronogram.nwk") # output NEWICK
OUT_CSV   = Path("pheidole-chronogram-mapping-used.csv")    # audit of mapping actually used

# ---- Column name fallbacks (case-sensitive matching here) ----
CODE_COL_CANDIDATES     = ["Extraction Code", "Code", "ID", "Sample", "TaxonCode", "Short", "Label"]
TAXON_COL_CANDIDATES    = ["Taxon", "Name", "FullName"]
SPECIMEN_COL_CANDIDATES = ["Specimen Code", "Specimen", "Voucher", "Accession"]

def pick_column(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Could not find any of the columns {candidates} in: {list(df.columns)}")
    return None

def main():
    # --- Read mapping table (TSV) ---
    df = pd.read_csv(MAP_PATH, sep="\t", dtype=str, keep_default_na=False)
    # Trim whitespace from column names and cells
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()

    code_col     = pick_column(df, CODE_COL_CANDIDATES)
    taxon_col    = pick_column(df, TAXON_COL_CANDIDATES)
    specimen_col = pick_column(df, SPECIMEN_COL_CANDIDATES)

    # Convert Taxon dots to underscores; build final label
    taxon_us = df[taxon_col].str.replace(".", "_", regex=False)
    df["NewLabel"] = taxon_us + "_" + df[specimen_col]

    # Build mapping dict: code -> NewLabel
    code_to_label = dict(zip(df[code_col], df["NewLabel"]))

    # --- Read NEWICK as raw text ---
    tree_txt = TREE_PATH.read_text(encoding="utf-8", errors="replace")

    # --- Replace labels safely ---
    # Match a taxon label only when it is *immediately* followed by ":" (NEWICK branch length)
    # and is preceded by "(" or "," (start of a tip/internal label).
    # We only attempt to replace exact codes present in the table.
    keys_sorted = sorted(code_to_label, key=len, reverse=True)  # longest first avoids partial overlaps
    if not keys_sorted:
        raise RuntimeError("Mapping table produced an empty code-to-label mapping.")

    pattern = r'(?<=[,(])(' + "|".join(map(re.escape, keys_sorted)) + r')(?=:)'

    replaced = {}  # code -> new label (only those found in tree)
    def _repl(m):
        code = m.group(1)
        new = code_to_label.get(code, code)
        if new != code:
            replaced[code] = new
        return new

    new_tree_txt, n_subs = re.subn(pattern, _repl, tree_txt)

    # --- Sanity checks and reporting ---
    # Collect any code-like labels from the original tree (2+ letters + digits) for diagnostics
    code_like = set(re.findall(r'(?<=[,(])([A-Za-z]{2,}\d+)(?=:)', tree_txt))
    unmapped_in_tree = sorted([c for c in code_like if c not in code_to_label])

    # Parentheses balance check (basic)
    if new_tree_txt.count("(") != new_tree_txt.count(")"):
        raise RuntimeError("Parentheses are unbalanced after replacement. Check input formatting.")

    # Write outputs
    OUT_TREE.write_text(new_tree_txt, encoding="utf-8")
    used_df = df[df[code_col].isin(replaced.keys())].copy()
    used_df = used_df[[code_col, taxon_col, specimen_col, "NewLabel"]]
    used_df.rename(columns={code_col: "Code", taxon_col: "Taxon", specimen_col: "Specimen"}, inplace=True)
    used_df.to_csv(OUT_CSV, index=False)

    # Summary
    print("=== Rename summary ===")
    print(f" Mapping table codes:     {len(code_to_label):,}")
    print(f" Codes seen in tree:      {len(code_like):,}")
    print(f" Labels replaced:         {len(replaced):,} (substitutions performed: {n_subs:,})")
    print(f" Unmapped codes in tree:  {len(unmapped_in_tree):,}")
    if unmapped_in_tree:
        print("  e.g.,", ", ".join(unmapped_in_tree[:20]), "...")
    # Show a few examples
    if replaced:
        ex = list(replaced.items())[:10]
        print("\n Examples:")
        for k, v in ex:
            print(f"  {k} -> {v}")

    print(f"\nWrote: {OUT_TREE} and {OUT_CSV}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("\n[ERROR]", e, file=sys.stderr)
        sys.exit(1)