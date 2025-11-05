import math
import random
import pathlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# fasta files
FASTA_FILES = [
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\corona.fasta",  
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\acute.fasta",
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\Boebeisivirus.fasta",
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\Enterobacteria.fasta",
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\Project_L6\L6\HepatitisB.fasta",
]

NUM_SAMPLES_PER_LANE = 10  
MIN_LEN = 100
MAX_LEN = 3000
RANDOM_SEED = None   
SAVE_PNG = None    


def read_fasta_one_sequence(path):
    """
    Read the FIRST sequence from a FASTA file.
    Returns (header, sequence) with sequence filtered to A/C/G/T.
    """
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
        raise ValueError(f"No A/C/G/T bases found in FASTA: {path}")
    return header, seq


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


def migration_positions_bp(lengths):
    logs = [math.log10(L) for L in lengths]
    min_log, max_log = min(logs), max(logs)
    if max_log == min_log:
        return [0.5] * len(lengths)  # all same size -> same position
    top_margin, bottom_margin = 0.1, 0.9
    usable = bottom_margin - top_margin
    positions = [top_margin + ((max_log - lg) / (max_log - min_log)) * usable for lg in logs]
    return positions


def plot_multi_gel(positions_list, lengths_list, lane_labels, save_path=None, title="Gel Simulation"):

    n_lanes = len(positions_list)
    assert n_lanes == len(lengths_list) == len(lane_labels)

    fig, ax = plt.subplots(figsize=(1.8 + 1.1 * n_lanes, 6.5))

    gel_left, gel_bottom, gel_width, gel_height = 0.06, 0.05, 0.88, 0.9
    ax.add_patch(Rectangle((gel_left, gel_bottom), gel_width, gel_height, fill=False, linewidth=1.5))

    lane_total_width = 0.7  
    lane_x_start = gel_left + 0.08
    lane_x_end = lane_x_start + lane_total_width
    lane_centers = [
        lane_x_start + (i + 0.5) * (lane_total_width / n_lanes) for i in range(n_lanes)
    ]
    lane_width = (lane_total_width / n_lanes) * 0.6  
    label_y = gel_bottom + 0.02

    for idx, (cent, band_positions, band_lengths) in enumerate(zip(lane_centers, positions_list, lengths_list)):
        x_left = cent - lane_width / 2
        x_right = cent + lane_width / 2

        ax.add_patch(Rectangle((x_left, gel_bottom + 0.03), lane_width, gel_height - 0.06, fill=False, linewidth=1.0))

        for y, L in sorted(zip(band_positions, band_lengths), key=lambda t: t[0]):
            ax.plot([x_left + 0.01, x_right - 0.01], [y, y], linewidth=4)
            ax.text(x_right + 0.015, y, f"{L} bp", va="center", fontsize=9)

        ax.text(cent, label_y, lane_labels[idx], ha="center", va="center", fontsize=10, rotation=0)

    if title:
        ax.set_title(title, fontsize=12, pad=12)

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

    if not FASTA_FILES:
        raise ValueError("Add 1–5 FASTA paths to FASTA_FILES.")

    lane_labels = []
    positions_list = []
    lengths_list = []

    for path in FASTA_FILES:
        p = pathlib.Path(path)
        if not p.exists():
            raise FileNotFoundError(f"FASTA file not found: {p}")

        header, seq = read_fasta_one_sequence(p)
        samples = take_random_samples(
            seq,
            num_samples=NUM_SAMPLES_PER_LANE,
            min_len=MIN_LEN,
            max_len=MAX_LEN
        )

        lens = [len(s) for s in samples]
        pos = migration_positions_bp(lens)

        lane_labels.append(p.stem)  
        positions_list.append(pos)
        lengths_list.append(lens)

    plot_multi_gel(
        positions_list=positions_list,
        lengths_list=lengths_list,
        lane_labels=lane_labels,
        save_path=SAVE_PNG,
        title="Multi-lane Gel Simulation"
    )
