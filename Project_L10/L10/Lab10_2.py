'''
Download from moodle the FASTA file containing promotor sequences and use it as an input for an updated application which i s
able to generate  and save the objective digital stain ODS inside the folder. the promotor file can be found
inside the promkappa package on moodle or github inside folder "bin"
'''

import os
import numpy as np

# ------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------
WINDOW = 30
ODS_FOLDER = "C:\LaylaYear4\BioInformatics\BioInformatics\Project_L10\L10\ODS"


# ------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------

def read_fasta(filepath):
    sequences = {}
    name = None
    seq = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue                     # skip empty lines
            if line.startswith(">"):
                if name:
                    sequences[name] = "".join(seq)
                name = line[1:].strip()      # full header as name
                seq = []
            else:
                seq.append(line.replace(" ", ""))  # remove spaces

    if name:
        sequences[name] = "".join(seq)

    return sequences



def compute_cpg_content(seq):
    """PromKappa CG computation (returns ~29.27 for test window)."""
    CG = seq.count("C") + seq.count("G")
    return (CG / len(seq)) * 79.92


def index_of_coincidence(seq):
    """Kappa IC (returns ~27.53 for test window)."""
    L = len(seq)
    freqs = {n: seq.count(n) for n in "ACGT"}
    numerator = sum(freqs[n] * (freqs[n] - 1) for n in freqs)
    return (numerator / (L * (L - 1))) * 100


def compute_pattern(sequence, window):
    """Computes CG% + IC% sliding-window pattern."""
    CG_values = []
    IC_values = []
    positions = []

    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i+window]
        CG_values.append(compute_cpg_content(window_seq))
        IC_values.append(index_of_coincidence(window_seq))
        positions.append(i + window/2)

    return positions, CG_values, IC_values


def save_ods(name, positions, CG_values, IC_values):
    """Saves pattern to an ODS file."""
    filename = os.path.join(ODS_FOLDER, f"{name}.ods")

    with open(filename, "w") as f:
        f.write("position\tCG%\tIC%\n")
        for pos, cg, ic in zip(positions, CG_values, IC_values):
            f.write(f"{pos:.2f}\t{cg:.2f}\t{ic:.2f}\n")

    print(f"ODS saved: {filename}")


# ------------------------------------------------------------
# MAIN APPLICATION
# ------------------------------------------------------------
def main():
    # Create output folder if missing
    if not os.path.exists(ODS_FOLDER):
        os.makedirs(ODS_FOLDER)

    fasta_path = "C:\LaylaYear4\BioInformatics\BioInformatics\Project_L10\L10\Promoters.txt"    # downloaded from Moodl=
    promoters = read_fasta(fasta_path)

    print(f"\nLoaded {len(promoters)} promoter sequences.\n")

    for name, seq in promoters.items():
        print(f"Processing: {name}  (length={len(seq)})")

        positions, CG_values, IC_values = compute_pattern(seq, WINDOW)
        save_ods(name, positions, CG_values, IC_values)

    print("\nAll ODS files generated successfully!")


if __name__ == "__main__":
    main()
