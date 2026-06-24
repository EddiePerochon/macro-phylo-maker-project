import re
import sys
from pathlib import Path

# Usage: python rename_nexus.py input.nex > output.nex
text = Path(sys.argv[1]).read_text(encoding="utf-8")

# 1) Build mapping: short_id -> NewName (Genus_species_Voucher) from namelong
#    Matches lines like: EGPxxxx[& ... namelong="EGPxxxx_VOUCHER_Genus.species_Country" ...]
taxa_pattern = re.compile(
    r'^\s*([A-Za-z0-9._-]+)\s*\[\s*&[^]]*?\bnamelong="([^"]+)"[^]]*]',
    re.MULTILINE
)

rename = {}
for short_id, namelong in taxa_pattern.findall(text):
    # Expected: SAMPLE_VOUCHER_Genus.species_Country
    parts = namelong.split("_")
    if len(parts) >= 3:
        # Safe parse: last part may be Country (ignored), middle can be multiple underscores if present
        sample = parts[0]
        voucher = parts[1]
        genus_species = parts[2]
        # If extra underscores exist after species (e.g. countries with underscores), ignore them
        # Species is like "Strumigenys.jugis" → "Strumigenys_jugis"
        if "." in genus_species:
            genus, species = genus_species.split(".", 1)
            new_label = f"{genus}_{species}_{voucher}"
            rename[short_id] = new_label

# 2) Replace labels globally (taxa block + tree block) respecting token boundaries
def replace_labels(s, mapping):
    # Sort by length descending to avoid prefix-collision (e.g., EGP001 vs EGP0019)
    for old in sorted(mapping.keys(), key=len, reverse=True):
        new = mapping[old]
        # Token boundary: do not replace substrings in larger tokens
        # Allowed token chars in Newick/Nexus labels often include [A-Za-z0-9._-]
        pattern = re.compile(rf'(?<![A-Za-z0-9._-]){re.escape(old)}(?![A-Za-z0-9._-])')
        s = pattern.sub(new, s)
    return s

out = replace_labels(text, rename)
sys.stdout.write(out)