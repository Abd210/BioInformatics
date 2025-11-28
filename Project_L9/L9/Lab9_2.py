'''
Download influenza genomes & use the restriction enzymes from your set to 
simulate the electrophorese=is job for these viral genomes.

Eliminate all commonalities from the 10 electrophoresis gel genomes and make 1 general
elecrophoresis gel that shows only the differences between these genomes
'''

import re
import math
from dataclasses import dataclass
from typing import List, Dict

import matplotlib.pyplot as plt


# ========= CONFIGURATION =========

# Put your 10 influenza genome FASTA files here:
GENOME_FILES = {
    "Flu1":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\Boebeisivirus.fasta",
    "Flu2":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\acute.fasta",
    "Flu3":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\HepatitisB.fasta",
    "Flu4":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\HepatitisC.fasta",
    "Flu5":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\papillomavirus.fasta",
    "Flu6":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\Nile.fasta",
    "Flu7":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\Measles.fasta",
    "Flu8":  r"C:\BioInformatics\BioInformatics\Project_L7\L7\Zika.fasta",
    "Flu9":  r"C:\BioInformatics\BioInformatics\Project_L8\L8\corona.fasta",
    "Flu10": r"C:\BioInformatics\BioInformatics\Project_L7\L7\Enterobacteria.fasta",
}



@dataclass
class RestrictionEnzyme:
    name: str
    site: str          # recognition sequence (5'→3')
    cut_offset: int    # position (in bases) from the start of site where the cut occurs


# Enzymes (cut offsets chosen to match your diagrams: G|AATTC etc.)
ENZYMES = [
    RestrictionEnzyme("EcoRI",  "GAATTC", 1),  # G|AATTC
    RestrictionEnzyme("BamHI",  "GGATCC", 1),  # G|GATCC
    RestrictionEnzyme("HindIII","AAGCTT", 1),  # A|AGCTT
    RestrictionEnzyme("TaqI",   "TCGA",   1),  # T|CGA
    RestrictionEnzyme("HaeIII", "GGCC",   2),  # GG|CC
]


# ========= FASTA PARSING =========

def read_fasta_sequence(path: str) -> str:
    """Read the first sequence from a FASTA file and return it as a single uppercase string."""
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # header line – skip, assume single-sequence FASTA
                continue
            seq_lines.append(line)
    sequence = "".join(seq_lines).upper()
    return sequence


# ========= RESTRICTION DIGEST LOGIC =========

def find_cut_positions(seq: str, enzyme: RestrictionEnzyme) -> List[int]:
    """
    Find all cut positions (0-based indices) for a given enzyme in the sequence.
    Uses lookahead regex so overlapping sites are also found.
    Cut coordinate = start of site + cut_offset.
    """
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
    """Given sequence length and sorted cut positions, return fragment lengths."""
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


def combined_digest_fragments(seq: str, enzymes: List[RestrictionEnzyme]) -> List[int]:
    """
    Digest a sequence with ALL enzymes at once (combined digest).
    Returns the list of fragment sizes.
    """
    seq_len = len(seq)
    all_cuts = set()

    for enzyme in enzymes:
        cuts = find_cut_positions(seq, enzyme)
        all_cuts.update(cuts)

    all_cuts_sorted = sorted(all_cuts)
    fragments = fragments_from_cuts(seq_len, all_cuts_sorted)
    return fragments


# ========= GEL SIMULATION (PLOT) =========

def simulate_gel(fragments_dict: Dict[str, List[int]], height: int = 200, title: str = "Simulated agarose gel") -> None:
    """
    Plot a simple agarose gel using matplotlib.
    Each key in fragments_dict is a lane (name), value is list of fragment sizes.
    Larger fragments stay near the top, smaller go toward the bottom.
    """
    # Flatten all fragment sizes to determine scaling
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

    # Prepare an empty "image" for the gel (height x num_lanes)
    # 0 = background, 1 = band
    gel = [[0.0 for _ in range(num_lanes)] for _ in range(height)]

    # For each lane, map fragment sizes to row indices and mark bands
    for lane_idx, name in enumerate(lane_names):
        sizes = fragments_dict[name]
        for size in sizes:
            size = max(size, 1)
            log_s = math.log10(size)
            norm = (log_s - log_min) / (log_max - log_min)  # 0..1
            # Big fragments (high log_s) stay near top (row 0)
            row = int((1.0 - norm) * (height - 1))

            # Draw a small vertical "thick band" around this row
            for r in range(max(0, row - 2), min(height, row + 3)):
                gel[r][lane_idx] = 1.0

    # Plot the gel
    plt.figure(figsize=(7, 8))
    plt.imshow(gel, aspect="auto", cmap="gray_r", interpolation="nearest")
    plt.gca().invert_yaxis()  # top = wells

    plt.xticks(range(num_lanes), lane_names, rotation=45, ha="right")
    plt.yticks([])

    plt.xlabel("Genome")
    plt.title(title)
    plt.tight_layout()
    plt.show()


# ========= FILTERING COMMON BANDS =========

def bin_size(size: int, bin_width: int) -> int:
    """
    Convert a fragment size (bp) to a bin index, so we can consider
    bands within ~bin_width bp as equivalent.
    """
    return size // bin_width


def compute_difference_fragments(
    fragments_per_genome: Dict[str, List[int]],
    bin_width: int = 20
) -> Dict[str, List[int]]:
    """
    Remove bands that are 'common' to ALL genomes (within a size bin).
    Returns new dict where each genome only has bands that are not in the
    global intersection of bins.

    bin_width: resolution of the gel in bp (e.g. 20 bp).
    """
    # For each genome, compute set of bins it occupies
    genome_bins: Dict[str, List[int]] = {}
    bin_sets: Dict[str, set] = {}

    for name, frags in fragments_per_genome.items():
        bins = [bin_size(s, bin_width) for s in frags]
        genome_bins[name] = bins
        bin_sets[name] = set(bins)

    # Bins that appear in ALL genomes (common bands)
    all_names = list(fragments_per_genome.keys())
    common_bins = set(bin_sets[all_names[0]])
    for name in all_names[1:]:
        common_bins &= bin_sets[name]

    print(f"Common band bins across all genomes (bin_width={bin_width} bp): {sorted(common_bins)}")

    # Now remove fragments that fall into common bins
    diff_fragments: Dict[str, List[int]] = {}
    for name, frags in fragments_per_genome.items():
        bins = genome_bins[name]
        keep_sizes = [
            size for size, b in zip(frags, bins)
            if b not in common_bins
        ]
        diff_fragments[name] = keep_sizes

    return diff_fragments


# ========= MAIN SCRIPT =========

def main():
    # 1. Load all influenza genomes and digest them with the enzyme set
    fragments_per_genome: Dict[str, List[int]] = {}

    for genome_name, fasta_path in GENOME_FILES.items():
        seq = read_fasta_sequence(fasta_path)
        print(f"{genome_name}: sequence length = {len(seq)} bp")

        fragments = combined_digest_fragments(seq, ENZYMES)
        fragments.sort(reverse=True)  # sort big-to-small, optional
        fragments_per_genome[genome_name] = fragments

        print(f"  Number of fragments: {len(fragments)}")
        print(f"  First 10 fragment sizes (bp): {fragments[:10]}\n")

    # 2. Gel with ALL bands for all genomes
    simulate_gel(
        fragments_per_genome,
        height=200,
        title="Influenza genomes - combined digest (all bands)"
    )

    # 3. Remove common bands across all genomes (within tolerance)
    # Adjust bin_width to change how strict "same band" is
    diff_fragments = compute_difference_fragments(
        fragments_per_genome,
        bin_width=20   # ~20 bp resolution
    )

    # 4. Gel that shows ONLY differences between genomes
    simulate_gel(
        diff_fragments,
        height=200,
        title="Influenza genomes - differences only (common bands removed)"
    )


if __name__ == "__main__":
    main()
