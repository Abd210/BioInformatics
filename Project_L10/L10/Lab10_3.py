'''
Download 10 influenza genomes and plot their objective digital stains. Plot the the center of each digital stain on a 2nd plot and 
label the slots with their name/version
'''

import os
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# YOUR GENOME PATHS
# ------------------------------------------------------------
GENOME_FILES = {
    "Flu1":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Boebeisivirus.fasta",
    "Flu2":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\acute.fasta",
    "Flu3":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\HepatitisB.fasta",
    "Flu4":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\HepatitisC.fasta",
    "Flu5":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\papillomavirus.fasta",
    "Flu6":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Nile.fasta",
    "Flu7":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Measles.fasta",
    "Flu8":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Zika.fasta",
    "Flu9":  r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L8\L8\corona.fasta",
    "Flu10": r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L7\L7\Enterobacteria.fasta",
}

WINDOW = 30


# ------------------------------------------------------------
# FASTA READER
# ------------------------------------------------------------
def read_fasta(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return "".join(seq)


# ------------------------------------------------------------
# CG% (PromKappa)
# ------------------------------------------------------------
def compute_cpg_content(seq):
    CG = seq.count("C") + seq.count("G")
    return (CG / len(seq)) * 79.92  # scaled for PromKappa look


# ------------------------------------------------------------
# Index of Coincidence (PromKappa)
# ------------------------------------------------------------
def index_of_coincidence(seq):
    L = len(seq)
    freqs = {n: seq.count(n) for n in "ACGT"}
    numerator = sum(freqs[n] * (freqs[n] - 1) for n in freqs)
    return (numerator / (L * (L - 1))) * 100


# ------------------------------------------------------------
# Compute full pattern for genome
# ------------------------------------------------------------
def compute_pattern(sequence, window):
    CG_values = []
    IC_values = []
    centers = []

    for i in range(len(sequence) - window + 1):
        win = sequence[i:i+window]
        CG_values.append(compute_cpg_content(win))
        IC_values.append(index_of_coincidence(win))
        centers.append(i + window/2)

    return np.array(centers), np.array(CG_values), np.array(IC_values)


# ------------------------------------------------------------
# MAIN EXECUTION
# ------------------------------------------------------------
center_points = []

plt.figure(figsize=(16, 10))
plt.title("Objective Digital Stains for 10 Genomes")
plt.xlabel("Genome Position")
plt.ylabel("ODS Value (CG%, IC%)")

for name, path in GENOME_FILES.items():
    print(f"Loading {name} from {path}")
    seq = read_fasta(path)

    centers, CG_vals, IC_vals = compute_pattern(seq, WINDOW)

    # Plot only CG% (or both if you want)
    plt.plot(centers, CG_vals, label=name, linewidth=2)

    # center of weight for this genome
    Cw = np.sum(centers * CG_vals) / np.sum(CG_vals)
    center_points.append((name, Cw))

plt.legend()
plt.grid(True)
plt.show()


# ------------------------------------------------------------
# PLOT THE CENTERS (2nd Plot)
# ------------------------------------------------------------
plt.figure(figsize=(12, 6))
plt.title("Centers of Weight of the 10 Influenza Genomes")
plt.xlabel("Genome Index (1–10)")
plt.ylabel("Center Weight Position")

for i, (name, cw) in enumerate(center_points, start=1):
    plt.scatter(i, cw, s=100)
    plt.text(i + 0.05, cw, name, fontsize=10)

plt.grid(True)
plt.show()
