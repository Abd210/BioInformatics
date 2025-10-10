

#!/usr/bin/env python3
"""
FASTA sliding-window nucleotide frequency viewer

- GUI (Tkinter) to select a FASTA file
- Sliding window (default 30; user editable)
- Computes relative frequencies of A, C, G, T within each window
  (denominator counts only A/C/G/T; ambiguous chars like N are ignored)
- Plots four signals (A, C, G, T) using Matplotlib embedded in Tkinter
"""

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import os
import re
from collections import Counter

# Matplotlib in Tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

DNA_ALPHABET = ("A", "C", "G", "T")

def read_fasta(filepath: str) -> str:
    """
    Minimal FASTA reader: concatenates all sequence lines, ignores headers.
    Returns uppercase sequence with whitespace removed.
    """
    seq_parts = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(re.sub(r"\s+", "", line))
    return "".join(seq_parts).upper()

def sliding_window_frequencies(seq: str, k: int, alphabet=DNA_ALPHABET):
    """
    Compute per-window relative frequencies for each symbol in `alphabet`.

    Denominator for each window = number of characters inside the window
    that are in `alphabet` (so ambiguous bases like N are ignored).
    If denominator is 0 for a window, all frequencies are 0 for that window.

    Returns:
        xs: list of window indices (0-based start positions)
        freq_vectors: dict {symbol: list[float]} with one value per window
    """
    n = len(seq)
    if n < k:
        return [], {a: [] for a in alphabet}

    xs = list(range(0, n - k + 1))
    freq_vectors = {a: [] for a in alphabet}

    for start in xs:
        window = seq[start:start + k]
        # count only alphabet symbols
        filtered = [ch for ch in window if ch in alphabet]
        denom = len(filtered)

        if denom == 0:
            for a in alphabet:
                freq_vectors[a].append(0.0)
            continue

        counts = Counter(filtered)
        for a in alphabet:
            freq_vectors[a].append(counts.get(a, 0) / denom)

    return xs, freq_vectors

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-Window Nucleotide Frequencies (FASTA)")
        self.geometry("1000x700")

        self.sequence = ""
        self.filepath = None

        # --- Top controls ---
        top = ttk.Frame(self, padding=10)
        top.pack(side=tk.TOP, fill=tk.X)

        self.file_label_var = tk.StringVar(value="No file selected")
        ttk.Button(top, text="Open FASTA…", command=self.open_fasta).pack(side=tk.LEFT)
        ttk.Label(top, textvariable=self.file_label_var).pack(side=tk.LEFT, padx=8)

        ttk.Label(top, text="Window size:").pack(side=tk.LEFT, padx=(20, 5))
        self.window_var = tk.StringVar(value="30")
        self.window_entry = ttk.Entry(top, textvariable=self.window_var, width=6)
        self.window_entry.pack(side=tk.LEFT)

        ttk.Button(top, text="Analyze & Plot", command=self.analyze_and_plot).pack(side=tk.LEFT, padx=12)

        # --- Figure area ---
        self.fig = Figure(figsize=(10, 5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Relative nucleotide frequencies per sliding window")
        self.ax.set_xlabel("Window start index (0-based)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # --- Status bar ---
        self.status_var = tk.StringVar(value="Open a FASTA file to begin.")
        status = ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor="w", padding=5)
        status.pack(side=tk.BOTTOM, fill=tk.X)

    def open_fasta(self):
        path = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA files", "*.fa *.fasta *.fna *.ffn *.faa *.frn"), ("All files", "*.*")]
        )
        if not path:
            return
        try:
            seq = read_fasta(path)
            if not seq:
                raise ValueError("No sequence content found in the file.")
            self.sequence = seq
            self.filepath = path
            self.file_label_var.set(os.path.basename(path))
            self.status_var.set(f"Loaded sequence length: {len(seq)} bases.")
        except Exception as e:
            messagebox.showerror("Error reading FASTA", str(e))
            self.status_var.set("Failed to load FASTA.")

    def analyze_and_plot(self):
        if not self.sequence:
            messagebox.showwarning("No sequence", "Please open a FASTA file first.")
            return

        # Validate window size
        try:
            k = int(self.window_var.get())
            if k <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid window size", "Window size must be a positive integer.")
            return

        if len(self.sequence) < k:
            messagebox.showwarning(
                "Sequence too short",
                f"Sequence length ({len(self.sequence)}) is shorter than window size ({k})."
            )
            return

        # Compute frequency vectors first (as requested), then plot
        xs, freq_vectors = sliding_window_frequencies(self.sequence, k, DNA_ALPHABET)

        # --- Plot ---
        self.ax.clear()
        self.ax.set_title("Relative nucleotide frequencies per sliding window")
        self.ax.set_xlabel("Window start index (0-based)")
        self.ax.set_ylabel("Relative frequency")
        self.ax.set_ylim(0, 1)

        # Plot in fixed A, C, G, T order to ensure four signals
        for base in DNA_ALPHABET:
            self.ax.plot(xs, freq_vectors[base], label=base)

        self.ax.legend(title="Nucleotide")
        self.ax.grid(True, linestyle="--", alpha=0.4)

        self.canvas.draw_idle()

        # Status update
        self.status_var.set(
            f"Plotted {len(xs)} windows (size {k}). "
            f"Averages — A:{avg(freq_vectors['A']):.3f} C:{avg(freq_vectors['C']):.3f} "
            f"G:{avg(freq_vectors['G']):.3f} T:{avg(freq_vectors['T']):.3f}"
        )

def avg(vals):
    return sum(vals) / len(vals) if vals else 0.0

if __name__ == "__main__":
    # Nice platform-consistent theming
    try:
        import tkinter.ttk as ttk  # already imported; ensures availability
        # Use "clam" or "vista"/"xpnative" if available
        root = App()
        root.mainloop()
    except Exception as e:
        # Fallback plain run if theming errors occur
        app = App()
        app.mainloop()

