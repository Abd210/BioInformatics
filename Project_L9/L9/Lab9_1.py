import re
import math
from dataclasses import dataclass
from typing import List, Dict

import matplotlib.pyplot as plt


FASTA_PATH = r"C:\BioInformatics\BioInformatics\Project_L8\L8\corona.fasta" 


@dataclass
class RestrictionEnzyme:
    name: str
    site: str          
    cut_offset: int   

ENZYMES = [
    RestrictionEnzyme("EcoRI",  "GAATTC", 1),  # G|AATTC
    RestrictionEnzyme("BamHI",  "GGATCC", 1),  # G|GATCC
    RestrictionEnzyme("HindIII","AAGCTT", 1),  # A|AGCTT
    RestrictionEnzyme("TaqI",   "TCGA",   1),  # T|CGA
    RestrictionEnzyme("HaeIII", "GGCC",   2),  # GG|CC
]



def read_fasta_sequence(path: str) -> str:
    """Read the first sequence from a FASTA file and return it as a single uppercase string."""
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line)
    sequence = "".join(seq_lines).upper()
    return sequence



def find_cut_positions(seq: str, enzyme: RestrictionEnzyme) -> List[int]:

    pattern = "(?=" + re.escape(enzyme.site) + ")"
    cut_positions = []
    for m in re.finditer(pattern, seq):
        start = m.start()
        cut_pos = start + enzyme.cut_offset
        # Only consider cuts strictly within the sequence
        if 0 < cut_pos < len(seq):
            cut_positions.append(cut_pos)
    return sorted(cut_positions)


def fragments_from_cuts(seq_length: int, cuts: List[int]) -> List[int]:
    if not cuts:
        return [seq_length]

    cuts_sorted = sorted(cuts)
    fragment_lengths = []

    # from start to first cut
    fragment_lengths.append(cuts_sorted[0] - 0)

    # between internal cuts
    for i in range(len(cuts_sorted) - 1):
        fragment_lengths.append(cuts_sorted[i+1] - cuts_sorted[i])

    # last cut to end
    fragment_lengths.append(seq_length - cuts_sorted[-1])

    return fragment_lengths



def simulate_gel(fragments_dict: Dict[str, List[int]], height: int = 200) -> None:
    all_sizes = [size for frags in fragments_dict.values() for size in frags]
    if not all_sizes:
        print("No fragments to display in gel.")
        return

    min_size = max(min(all_sizes), 1)
    max_size = max(all_sizes)

    if max_size == min_size:
        max_size += 1

    log_min = math.log10(min_size)
    log_max = math.log10(max_size)

    lane_names = list(fragments_dict.keys())
    num_lanes = len(lane_names)

 
    gel = [[0.0 for _ in range(num_lanes)] for _ in range(height)]

    for lane_idx, name in enumerate(lane_names):
        sizes = fragments_dict[name]
        for size in sizes:
            size = max(size, 1)
            log_s = math.log10(size)
            norm = (log_s - log_min) / (log_max - log_min)  
            row = int((1.0 - norm) * (height - 1))

            for r in range(max(0, row - 2), min(height, row + 3)):
                gel[r][lane_idx] = 1.0

    plt.figure(figsize=(6, 8))
    plt.imshow(gel, aspect="auto", cmap="gray_r", interpolation="nearest")
    plt.gca().invert_yaxis()  # top = wells

    plt.xticks(range(num_lanes), lane_names, rotation=45, ha="right")
    plt.yticks([])

    plt.xlabel("Enzyme / digest")
    plt.title("Simulated agarose gel")
    plt.tight_layout()
    plt.show()



def main():
    seq = read_fasta_sequence(FASTA_PATH)
    seq_len = len(seq)
    print(f"Sequence length: {seq_len} bp\n")

    all_digest_fragments: Dict[str, List[int]] = {}  
    combined_cuts = set()

    for enzyme in ENZYMES:
        print(f"=== {enzyme.name} ===")
        print(f"Recognition site: {enzyme.site}")
        print(f"Cut offset from 5' end of site: {enzyme.cut_offset}\n")

        cuts = find_cut_positions(seq, enzyme)
        combined_cuts.update(cuts)

        cuts_1based = [c + 1 for c in cuts]

        print(f"Number of cleavages: {len(cuts)}")
        if cuts_1based:
            print(f"Cleavage positions (1-based): {cuts_1based}")
        else:
            print("No cleavage sites found.")

        fragments = fragments_from_cuts(seq_len, cuts)
        all_digest_fragments[enzyme.name] = fragments

        print(f"Fragment lengths (bp): {fragments}\n")

    print("=== Combined digest (all enzymes) ===")
    combined_cuts_sorted = sorted(combined_cuts)
    combined_fragments = fragments_from_cuts(seq_len, combined_cuts_sorted)
    all_digest_fragments["All"] = combined_fragments

    print(f"Number of total cuts: {len(combined_cuts_sorted)}")
    print(f"All cleavage positions (1-based): {[c+1 for c in combined_cuts_sorted]}")
    print(f"Combined fragment lengths (bp): {combined_fragments}\n")

    simulate_gel(all_digest_fragments, height=200)


if __name__ == "__main__":
    main()
