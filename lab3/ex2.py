# Sliding Window Melting Temperature (Tm) Plotter
# Window = 9
# Calculates both:
#   Temp1 = 4(G + C) + 2(A + T)
#   Temp2 = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/length
# and plots both vs sequence position.

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
import math

WINDOW = 9  # sliding window size
NA_CONC = 0.001  # [Na+] mol/L

def read_fasta(path):
    """Read the first record from a FASTA file."""
    header = None
    seq_parts = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is None:
                    header = line[1:].strip()
                else:
                    break
            else:
                seq_parts.append(line)
    if header is None:
        raise ValueError("No FASTA header found.")
    return header, ''.join(seq_parts).upper()

def tm_wallace(seq):
    """Temp1: Wallace rule."""
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    return 4*(g + c) + 2*(a + t)

def tm_salt_adjusted(seq, na=NA_CONC):
    """Temp2: Salt-adjusted formula."""
    length = len(seq)
    gc_percent = (seq.count('G') + seq.count('C')) / length * 100
    return -(81.5 + 16.6*math.log10(na) + 0.41*gc_percent - 600/length)

def sliding_window_tm(seq, k=WINDOW):
    """Compute Temp1 and Temp2 for each sliding window."""
    n = len(seq)
    if n < k:
        raise ValueError(f"Sequence too short for window size ({k}).")
    positions = []
    temp1_vals = []
    temp2_vals = []
    for i in range(n - k + 1):
        win = seq[i:i + k]
        t1 = tm_wallace(win)
        t2 = tm_salt_adjusted(win)
        center = i + k//2 + 1  # 1-based center
        positions.append(center)
        temp1_vals.append(t1)
        temp2_vals.append(t2)
    return positions, temp1_vals, temp2_vals

def plot_tm(positions, temp1, temp2, title):
    plt.figure(figsize=(8,5))
    plt.plot(positions, temp1, label='Temp1 (Wallace)', color='red')
    plt.plot(positions, temp2, label='Temp2 (Salt-adjusted)', color='blue')
    plt.xlabel('Sequence position (nt)')
    plt.ylabel('Melting Temperature (°C)')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

def main():
    root = tk.Tk()
    root.withdraw()
    fasta_path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not fasta_path:
        messagebox.showinfo("Cancelled", "No file selected.")
        return

    try:
        header, seq = read_fasta(fasta_path)
        positions, temp1, temp2 = sliding_window_tm(seq)
        plot_tm(positions, temp1, temp2, f"Melting Temperature Profile — {header}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    main()
