# process-nylanderia.py
import re
import sys

def reformat_label(lbl: str) -> str:
    # 1) SRR reads: SRR2184160_N_MG01 -> Nylanderia_MG01_SRR2184160
    m = re.match(r'^SRR(\d+)_N_([\w]+)$', lbl)
    if m:
        return f"Nylanderia_{m.group(2)}_SRR{m.group(1)}"

    # 2) Ny<number> with locality + species at the end:
    #    Ny233_N_Vanuatu_sp1 -> Nylanderia_sp1_Vanuatu_Ny233
    m = re.match(r'^Ny(\d+)_N_([A-Za-z]+)_(sp[\w]+)$', lbl)
    if m:
        return f"Nylanderia_{m.group(3)}_{m.group(2)}_Ny{m.group(1)}"

    # 3) Ny<number> with "sp..." only (no locality):
    #    Ny220_N_sp1 -> Nylanderia_sp1_Ny220
    m = re.match(r'^Ny(\d+)_N_(sp[\w]+)$', lbl)
    if m:
        return f"Nylanderia_{m.group(2)}_Ny{m.group(1)}"

    # 4) Ny<number> with species or code-like tokens (FG1, CR1, Ecu2, Bra1, Ven3, etc.)
    #    Ny209_N_nuggeti      -> Nylanderia_nuggeti_Ny209
    #    Ny108_N_FG1          -> Nylanderia_FG1_Ny108
    #    Ny129_N_CR1          -> Nylanderia_CR1_Ny129
    m = re.match(r'^Ny(\d+)_N_([A-Za-z][\w]*)$', lbl)
    if m:
        return f"Nylanderia_{m.group(2)}_Ny{m.group(1)}"

    # 5) NylaM with species (including cf_, sp_05, etc.):
    #    NylaM53_flavipes     -> Nylanderia_flavipes_M53
    #    NylaM58_cf_fulva     -> Nylanderia_cf_fulva_M58
    #    NylaM20_sp_05        -> Nylanderia_sp_05_M20
    m = re.match(r'^NylaM(\d+)_([A-Za-z][\w]*)$', lbl)
    if m:
        return f"Nylanderia_{m.group(2)}_M{m.group(1)}"

    # 6) Fallbacks for missing species:
    #    NylaM101             -> Nylanderia_sp_M101
    m = re.match(r'^NylaM(\d+)$', lbl)
    if m:
        return f"Nylanderia_sp_M{m.group(1)}"

    #    Ny188_N              -> Nylanderia_sp_Ny188
    m = re.match(r'^Ny(\d+)_N$', lbl)
    if m:
        return f"Nylanderia_sp_Ny{m.group(1)}"

    # else unchanged
    return lbl

def fix_newick(text: str) -> str:
    # Replace only labels directly followed by :, ) or , (NEWICK-safe)
    return re.sub(r'([^,:()]+)(?=[:),])', lambda m: reformat_label(m.group(1)), text)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process-nylanderia.py input.tre output.tre")
        sys.exit(1)
    with open(sys.argv[1], 'r', encoding='utf-8') as f:
        s = f.read()
    out = fix_newick(s)
    with open(sys.argv[2], 'w', encoding='utf-8') as f:
        f.write(out)