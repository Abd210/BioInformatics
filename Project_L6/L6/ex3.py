import math
import random
import pathlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

FASTA_PATH = r"C:\Users\Abd\Desktop\sdm\BioInformatics\lab5\corona.fasta"
RANDOM_SEED = None     
SAVE_PNG = None        

#we put more enzymes and we choose from them
ENZYMES = {
    "EcoRI":  {"site": "GAATTC",   "cut_index": 1},  # G^AATTC
    "BamHI":  {"site": "GGATCC",   "cut_index": 1},  # G^GATCC
    "HindIII":{"site": "AAGCTT",   "cut_index": 1},  # A^AGCTT
    "NotI":   {"site": "GCGGCCGC", "cut_index": 2},  # GC^GGCCGC
    "XhoI":   {"site": "CTCGAG",   "cut_index": 1},  # C^TCGAG
    "PstI":   {"site": "CTGCAG",   "cut_index": 5},  # CTGCA^G
    "SacI":   {"site": "GAGCTC",   "cut_index": 5},  # GAGCT^C
    "KpnI":   {"site": "GGTACC",   "cut_index": 5},  # GGTAC^C
    "NcoI":   {"site": "CCATGG",   "cut_index": 1},  # C^CATGG
    "NheI":   {"site": "GCTAGC",   "cut_index": 1},  # G^CTAGC
    "SpeI":   {"site": "ACTAGT",   "cut_index": 1},  # A^CTAGT
    "XbaI":   {"site": "TCTAGA",   "cut_index": 1},  # T^CTAGA
}

def read_fasta_one_sequence(path):

    header = "sequence"
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_chunks:
                    break  # first sequence done
                header = line[1:].strip() or "sequence"
            else:
                seq_chunks.append(line.upper())
    seq = "".join(seq_chunks)
    seq = "".join(ch for ch in seq if ch in "ACGT")
    if not seq:
        raise ValueError(f"No A/C/G/T bases found in FASTA: {path}")
    return header, seq

def find_sites(seq, motif):
    sites = []
    m = len(motif)
    i = seq.find(motif, 0)
    while i != -1:
        sites.append(i)
        i = seq.find(motif, i + 1)
    return sites

def enzyme_cut_positions(seq, site, cut_index):

    positions = []
    for start in find_sites(seq, site):
        positions.append(start + cut_index)
    return positions

def pick_five_enzymes(seq, db, n=5):
    cutters = []
    non_cutters = []
    for name, info in db.items():
        cuts = enzyme_cut_positions(seq, info["site"], info["cut_index"])
        (cutters if cuts else non_cutters).append(name)
    random.shuffle(cutters)
    random.shuffle(non_cutters)
    chosen = cutters[:n]
    if len(chosen) < n:
        chosen += non_cutters[: (n - len(chosen))]
    return chosen

def digest_sequence(seq, enzymes, db):
    all_cuts = set()
    per_enzyme = {}
    for name in enzymes:
        site = db[name]["site"]
        cut_index = db[name]["cut_index"]
        positions = enzyme_cut_positions(seq, site, cut_index)
        per_enzyme[name] = sorted(positions)
        all_cuts.update(positions)

    cut_positions = sorted(all_cuts)
    coords = [0] + cut_positions + [len(seq)]
    frags = [coords[i+1] - coords[i] for i in range(len(coords)-1)]
    return sorted(frags, reverse=True), per_enzyme, cut_positions


def migration_positions_bp(lengths):
    logs = [math.log10(L) for L in lengths]
    min_log, max_log = min(logs), max(logs)
    if max_log == min_log:
        return [0.5] * len(lengths)
    top_margin, bottom_margin = 0.1, 0.9
    usable = bottom_margin - top_margin
    return [top_margin + ((max_log - lg) / (max_log - min_log)) * usable for lg in logs]

def plot_single_lane_gel(fragment_lengths, lane_label="Digest", title=None, save_path=None):

    positions = migration_positions_bp(fragment_lengths)
    fig, ax = plt.subplots(figsize=(4.5, 6.5))

    ax.add_patch(Rectangle((0.18, 0.05), 0.64, 0.9, fill=False, linewidth=1.5))

    lane_x_center = 0.5
    lane_width = 0.28
    x_left = lane_x_center - lane_width / 2
    x_right = lane_x_center + lane_width / 2

    ax.add_patch(Rectangle((x_left, 0.08), lane_width, 0.84, fill=False, linewidth=1.0))

    for y, L in sorted(zip(positions, fragment_lengths), key=lambda t: t[0]):
        ax.plot([x_left + 0.01, x_right - 0.01], [y, y], linewidth=4)
        ax.text(x_right + 0.03, y, f"{L} bp", va="center", fontsize=9)

    ax.text((x_left + x_right) / 2, 0.03, lane_label, ha="center", va="center", fontsize=10)
    if title:
        ax.set_title(title, fontsize=11, pad=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)
    else:
        random.seed()

    fasta_path = pathlib.Path(FASTA_PATH)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    header, seq = read_fasta_one_sequence(fasta_path)

    chosen = pick_five_enzymes(seq, ENZYMES, n=5)
    fragments, per_enzyme_sites, all_cut_positions = digest_sequence(seq, chosen, ENZYMES)

    print(f"Sequence: {header} (length {len(seq)} bp)")
    print("Chosen enzymes (5):", ", ".join(chosen))
    for name in chosen:
        site = ENZYMES[name]["site"]
        sites = per_enzyme_sites[name]
        print(f"  - {name} ({site}) -> {len(sites)} site(s): {sites if sites else '—'}")
    print(f"Total unique cut positions: {len(all_cut_positions)} -> fragments: {len(fragments)}")
    print("Fragment sizes (bp, largest first):", fragments)

    plot_single_lane_gel(
        fragment_lengths=fragments,
        lane_label="5-enzyme digest",
        title=f"Digest of {header}",
        save_path=SAVE_PNG
    )
