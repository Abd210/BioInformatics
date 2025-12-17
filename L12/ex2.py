# ex2.py

'''
Download 10 influenza genomes, adapt your app from the previous assignment in order to scan each genome for possible
motifs. For each genome, make a chart that shows the signal with most likely locations of real functional motif.
'''

import os
import math
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------
# A) Motif instances (from Ex1 / your screenshot)
# ----------------------------
seqs_all = [
    "GAGGTAAAC",  # row 1
    "TCCGTAAGT",  # row 2
    "CAGGTTGGA",  # row 3
    "ACAGTCAGT",  # row 4
    "TAGGTCATT",  # row 5
    "TAGGTACTG",  # row 6
    "ATGGTAACT",  # row 7
    "CAGGTATAC",  # row 8
    "TGTGTGAGT",  # row 9
    "AAGGTAAGT",  # row 10
]

USE_FIRST_N = 9          # change to 10 if your instructor wants all 10
pseudocount = 1.0        # smoothing; set to 0.0 only if required (risk of log(0))
alphabet = ["A", "C", "G", "T"]
null = {b: 0.25 for b in alphabet}

seqs = seqs_all[:USE_FIRST_N]
if not seqs:
    raise ValueError("No motif sequences provided.")
L = len(seqs[0])
if any(len(s) != L for s in seqs):
    raise ValueError("Motif sequences must be aligned and same length.")
if any(any(ch not in alphabet for ch in s) for s in seqs):
    raise ValueError("Motif sequences must contain only A/C/G/T.")

# ----------------------------
# B) Build log-likelihood matrix
# ----------------------------
counts = {b: [0] * L for b in alphabet}
for s in seqs:
    for j, ch in enumerate(s):
        counts[ch][j] += 1

count_df = pd.DataFrame(counts, index=range(1, L + 1)).T
count_df.columns = list(range(1, L + 1))  # <-- IMPORTANT: use int columns like your colleague

N = len(seqs)
denom = N + pseudocount * len(alphabet)
weight_df = (count_df + pseudocount) / denom

loglike_df = weight_df.copy()
for b in alphabet:
    loglike_df.loc[b] = [math.log(p / null[b]) for p in weight_df.loc[b]]

print("=== Log-likelihood Matrix ===")
print(loglike_df.round(3))
print()

# ----------------------------
# C) FASTA reading
# ----------------------------
def read_fasta(path: str) -> str:
    """Reads the first FASTA record and returns its sequence (uppercase, no spaces)."""
    seq_lines = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()

# ----------------------------
# D) Scan genome -> DataFrame(Position, Score)  (like your colleague)
# ----------------------------
def scan_sequence_df(seq: str, loglike: pd.DataFrame, motif_len: int) -> pd.DataFrame | None:
    if seq is None or len(seq) < motif_len:
        return None

    records = []
    for start in range(0, len(seq) - motif_len + 1):
        w = seq[start:start + motif_len]

        # skip ambiguous windows (N, etc.)
        if any(ch not in alphabet for ch in w):
            continue

        score = 0.0
        for pos, base in enumerate(w, start=1):
            score += float(loglike.loc[base, pos])  # pos is int column

        records.append((start + 1, score))  # 1-based position

    if not records:
        return None

    return pd.DataFrame(records, columns=["Position", "Score"])

# ----------------------------
# E) Plot in colleague style (percentile threshold + red points)
# ----------------------------
def plot_motif_signal(scores_df: pd.DataFrame, genome_name: str, threshold_percentile: int = 95):
    if scores_df is None or len(scores_df) == 0:
        print(f"No data to plot for {genome_name}")
        return

    finite_scores = scores_df["Score"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(finite_scores) == 0:
        print(f"No finite scores for {genome_name}")
        return

    threshold = np.percentile(finite_scores, threshold_percentile)
    high_scoring = scores_df[scores_df["Score"] >= threshold]

    plt.figure(figsize=(14, 6))

    # blue thin line
    plt.plot(
        scores_df["Position"],
        scores_df["Score"],
        linewidth=0.5,
        alpha=0.6,
        color="blue",
        label="Score signal"
    )

    # red points = candidates
    plt.scatter(
        high_scoring["Position"],
        high_scoring["Score"],
        color="red",
        s=20,
        alpha=0.7,
        zorder=5,
        label=f"High scores (≥{threshold_percentile}th percentile)"
    )

    # orange threshold line
    plt.axhline(
        y=threshold,
        color="orange",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label=f"Threshold ({threshold:.2f})"
    )

    plt.xlabel("Genome Position")
    plt.ylabel("Log-likelihood Score")
    plt.title(f"Motif Signal for {genome_name}\n({len(high_scoring)} candidate motif locations)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plot_filename = f"{genome_name.replace(' ', '_')}_motif_signal.png"
    plt.savefig(plot_filename, dpi=150)
    plt.close()

    print(f"Saved plot: {plot_filename}")

    print(f"\n=== Top 10 motif candidates for {genome_name} ===")
    print(scores_df.nlargest(10, "Score").to_string(index=False))
    print()

# ----------------------------
# F) Run on your 10 influenza genomes
# ----------------------------
def natural_key(filename: str):
    import re
    m = re.search(r"\((\d+)\)", filename)
    return int(m.group(1)) if m else filename

# safest for your naming
fasta_files = sorted(glob.glob("sequence (*.fasta"), key=natural_key)

if not fasta_files:
    fasta_files = [f"sequence ({i}).fasta" for i in range(1, 11) if os.path.exists(f"sequence ({i}).fasta")]

if not fasta_files:
    raise FileNotFoundError("No FASTA files found. Put 'sequence (1).fasta' ... 'sequence (10).fasta' in this folder.")

print(f"Motif length: {L}, sequences used: {N}, pseudocount: {pseudocount}")
print(f"Found {len(fasta_files)} FASTA files.\n")

summary_rows = []
for i, path in enumerate(fasta_files, 1):
    print(f"\n{'='*60}")
    print(f"Processing Genome {i}: {path}")
    print(f"{'='*60}")

    seq = read_fasta(path)
    print(f"Genome length: {len(seq)} bp")

    scores_df = scan_sequence_df(seq, loglike_df, L)

    genome_name = f"Influenza_Genome_{i}"
    plot_motif_signal(scores_df, genome_name, threshold_percentile=95)

    # best hit summary
    if scores_df is not None and len(scores_df) > 0:
        best_row = scores_df.loc[scores_df["Score"].idxmax()]
        best_pos = int(best_row["Position"])
        best_score = float(best_row["Score"])
        best_window = seq[best_pos - 1: best_pos - 1 + L]
    else:
        best_pos = None
        best_score = None
        best_window = None

    summary_rows.append({
        "file": path,
        "genome_length": len(seq),
        "best_start_1based": best_pos,
        "best_window": best_window,
        "best_score": best_score,
        "plot": f"{genome_name}_motif_signal.png"
    })

summary_df = pd.DataFrame(summary_rows)
print(f"\n{'='*60}")
print("SUMMARY")
print(f"{'='*60}")
print(summary_df.to_string(index=False))
