'''
 Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. 
 The results of gel electrophoresis indicate the relative sizes of fragments, 
 which is useful for restriction mapping and analyzing PCR fragments.
 	1.	Take an arbitrary DNA sequence from the NCBI, between 1000 and 3000 nucleotides.
 	2.	Take 10 random samples from this sequence, between 100–3000 bases.
 	3.	Store these samples in an array.
 	4.	Simulate the migration of these DNA fragments on the electrophoresis gel, based on their molecular weights – however, their length should be sufficient for this scale (show a visual representation).

NOTE:
 Short DNA fragments meet small friction and travel faster through the electrophoresis gel. Long DNA fragments exhibit a high friction force and travel slowly through the electrophoresis gel.
'''

import math
import random
import pathlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


FASTA_PATH = r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\corona.fasta"
NUM_SAMPLES = 10      
MIN_LEN = 100         
MAX_LEN = 3000        
LANE_LABEL = "Corona DNA"
RANDOM_SEED = None    
SAVE_PNG = None       

#1
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
                    break 
                header = line[1:].strip() or "sequence"
            else:
                seq_chunks.append(line.upper())
    seq = "".join(seq_chunks)
    seq = "".join(ch for ch in seq if ch in "ACGT")
    if not seq:
        raise ValueError("No A/C/G/T bases found in the provided FASTA.")
    return header, seq

#2
def take_random_samples(sequence, num_samples=10, min_len=100, max_len=3000):
    actual_max = min(max_len, len(sequence))
    if actual_max < min_len:
        raise ValueError(
            f"Parent sequence too short for requested minimum length "
            f"({len(sequence)} < {min_len})."
        )
    samples = []
    for _ in range(num_samples):
        L = random.randint(min_len, actual_max)
        start = random.randint(0, len(sequence) - L)
        samples.append(sequence[start:start + L])
    return samples

#3
def migration_positions_bp(lengths):

    logs = [math.log10(L) for L in lengths]
    min_log, max_log = min(logs), max(logs)
    if max_log == min_log:
        return [0.5] * len(lengths) 

    top_margin, bottom_margin = 0.1, 0.9
    usable = bottom_margin - top_margin
    positions = [top_margin + ((max_log - lg) / (max_log - min_log)) * usable for lg in logs]
    return positions

def plot_gel(positions, lengths, lane_label="Lane 1", title=None,
             lane_x_center=0.5, lane_width=0.28, save_path=None):

    fig, ax = plt.subplots(figsize=(4.5, 6.5))

    ax.add_patch(Rectangle((0.18, 0.05), 0.64, 0.9, fill=False, linewidth=1.5))

    x_left = lane_x_center - lane_width / 2
    x_right = lane_x_center + lane_width / 2

    ax.add_patch(Rectangle((x_left, 0.08), lane_width, 0.84, fill=False, linewidth=1.0))

    for y, L in sorted(zip(positions, lengths), key=lambda t: t[0]):  # top -> bottom
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

    header, full_seq = read_fasta_one_sequence(fasta_path)

    samples = take_random_samples(
        full_seq,
        num_samples=NUM_SAMPLES,
        min_len=MIN_LEN,
        max_len=MAX_LEN
    )

    sample_lengths = [len(s) for s in samples]
    band_positions = migration_positions_bp(sample_lengths)

    title = f"Gel Simulation — {header}"
    plot_gel(
        positions=band_positions,
        lengths=sample_lengths,
        lane_label=LANE_LABEL,
        title=title,
        save_path=SAVE_PNG
    )
