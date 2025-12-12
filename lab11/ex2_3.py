import os
import itertools
import threading
from bisect import bisect_left
from array import array

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

DNA_ALLOWED = set("ACGTN")
NEG_INF = -1_000_000_000


# ------------------------- FASTA IO -------------------------

def read_fasta_first_record(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(line)
    seq = "".join(seq_parts).upper()
    seq = "".join(c for c in seq if c in DNA_ALLOWED)
    return seq


# ------------------------- Needleman–Wunsch (small regions) -------------------------

def nw_small(s1: str, s2: str, match=1, mismatch=-1, gap=-1):
    n, m = len(s1), len(s2)
    trace = [bytearray(m + 1) for _ in range(n + 1)]
    prev = array('i', [0] * (m + 1))
    curr = array('i', [0] * (m + 1))

    trace[0][0] = 0
    for j in range(1, m + 1):
        prev[j] = prev[j - 1] + gap
        trace[0][j] = 2

    for i in range(1, n + 1):
        curr[0] = prev[0] + gap
        trace[i][0] = 1
        c1 = s1[i - 1]
        for j in range(1, m + 1):
            c2 = s2[j - 1]
            sub = match if (c1 == c2 and c1 != 'N') else mismatch

            diag = prev[j - 1] + sub
            up = prev[j] + gap
            left = curr[j - 1] + gap

            best = diag
            d = 0
            if up > best:
                best = up
                d = 1
            if left > best:
                best = left
                d = 2

            curr[j] = best
            trace[i][j] = d

        prev, curr = curr, prev

    final_score = prev[m]

    i, j = n, m
    a1, a2 = [], []
    matches = 0
    while i > 0 or j > 0:
        if i == 0:
            d = 2
        elif j == 0:
            d = 1
        else:
            d = trace[i][j]

        if d == 0:
            c1 = s1[i - 1]
            c2 = s2[j - 1]
            a1.append(c1)
            a2.append(c2)
            if c1 == c2 and c1 != 'N':
                matches += 1
            i -= 1
            j -= 1
        elif d == 1:
            a1.append(s1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(s2[j - 1])
            j -= 1

    a1.reverse()
    a2.reverse()
    return "".join(a1), "".join(a2), final_score, matches


# ------------------------- SAFE banded fallback -------------------------

def banded_nw(s1: str, s2: str, match=1, mismatch=-1, gap=-1, band=800):
    n, m = len(s1), len(s2)
    band = max(band, abs(n - m) + 300)

    j_start = [0] * (n + 1)
    j_end = [0] * (n + 1)
    trace_rows = [None] * (n + 1)

    for i in range(n + 1):
        js = max(0, i - band)
        je = min(m, i + band)
        if js > je:
            raise ValueError("Band too small. Increase band.")
        j_start[i] = js
        j_end[i] = je
        trace_rows[i] = bytearray(je - js + 1)

    prev_js = j_start[0]
    prev_je = j_end[0]
    prev = array('i', [NEG_INF] * (prev_je - prev_js + 1))
    for j in range(prev_js, prev_je + 1):
        prev[j - prev_js] = gap * j
        trace_rows[0][j - prev_js] = 2

    for i in range(1, n + 1):
        curr_js = j_start[i]
        curr_je = j_end[i]
        curr = array('i', [NEG_INF] * (curr_je - curr_js + 1))
        c1 = s1[i - 1]

        for j in range(curr_js, curr_je + 1):
            idx = j - curr_js
            if j == 0:
                curr[idx] = gap * i
                trace_rows[i][idx] = 1
                continue

            diag = NEG_INF
            pj = j - 1
            if prev_js <= pj <= prev_je:
                c2 = s2[j - 1]
                sub = match if (c1 == c2 and c1 != 'N') else mismatch
                diag = prev[pj - prev_js] + sub

            up = NEG_INF
            if prev_js <= j <= prev_je:
                up = prev[j - prev_js] + gap

            left = NEG_INF
            cj = j - 1
            if curr_js <= cj <= curr_je:
                left = curr[cj - curr_js] + gap

            best = diag
            d = 0
            if up > best:
                best = up
                d = 1
            if left > best:
                best = left
                d = 2

            curr[idx] = best
            trace_rows[i][idx] = d

        prev, prev_js, prev_je = curr, curr_js, curr_je

    final_score = prev[m - prev_js]

    i, j = n, m
    a1, a2 = [], []
    matches = 0
    while i > 0 or j > 0:
        js = j_start[i]
        d = trace_rows[i][j - js]
        if i == 0:
            d = 2
        elif j == 0:
            d = 1

        if d == 0:
            c1 = s1[i - 1]
            c2 = s2[j - 1]
            a1.append(c1)
            a2.append(c2)
            if c1 == c2 and c1 != 'N':
                matches += 1
            i -= 1
            j -= 1
        elif d == 1:
            a1.append(s1[i - 1])
            a2.append('-')
            i -= 1
        else:
            a1.append('-')
            a2.append(s2[j - 1])
            j -= 1

    a1.reverse()
    a2.reverse()
    return "".join(a1), "".join(a2), final_score, matches


# ------------------------- Anchors -------------------------

def build_kmer_index(seq: str, k: int, occ_limit: int):
    index = {}
    counts = {}
    L = len(seq)
    for pos in range(0, L - k + 1):
        kmer = seq[pos:pos + k]
        c = counts.get(kmer, 0) + 1
        counts[kmer] = c
        if c > occ_limit:
            continue
        index.setdefault(kmer, []).append(pos)

    for kmer, c in list(counts.items()):
        if c > occ_limit and kmer in index:
            del index[kmer]
    return index


def lis_chain(matches):
    if not matches:
        return []
    matches.sort(key=lambda x: (x[0], -x[1]))
    tails = []
    tails_idx = []
    prev = [-1] * len(matches)

    for idx, (_, j) in enumerate(matches):
        pos = bisect_left(tails, j)
        if pos == len(tails):
            tails.append(j)
            tails_idx.append(idx)
        else:
            tails[pos] = j
            tails_idx[pos] = idx
        if pos > 0:
            prev[idx] = tails_idx[pos - 1]

    chain = []
    k = tails_idx[-1]
    while k != -1:
        chain.append(matches[k])
        k = prev[k]
    chain.reverse()
    return chain


def extend_exact_run(s1: str, s2: str, i: int, j: int, k: int):
    n, m = len(s1), len(s2)
    i0, j0 = i, j
    i1, j1 = i + k, j + k

    while i0 > 0 and j0 > 0 and s1[i0 - 1] == s2[j0 - 1]:
        i0 -= 1
        j0 -= 1
    while i1 < n and j1 < m and s1[i1] == s2[j1]:
        i1 += 1
        j1 += 1
    return i0, j0, (i1 - i0)


def find_anchors(s1: str, s2: str, k=13, min_run=30, occ_limit=40, max_anchors=120):
    if len(s1) < k or len(s2) < k:
        return []
    if len(s2) <= len(s1):
        idx_seq, scan_seq, swapped = s2, s1, False
    else:
        idx_seq, scan_seq, swapped = s1, s2, True

    idx = build_kmer_index(idx_seq, k, occ_limit=occ_limit)

    hits = []
    for i in range(0, len(scan_seq) - k + 1):
        kmer = scan_seq[i:i + k]
        pos_list = idx.get(kmer)
        if not pos_list:
            continue
        for j in pos_list:
            hits.append((j, i) if swapped else (i, j))

    chain = lis_chain(hits)
    if not chain:
        return []

    anchors = []
    last_i_end = -1
    last_j_end = -1
    for (i, j) in chain:
        i0, j0, run = extend_exact_run(s1, s2, i, j, k)
        if run < min_run:
            continue
        if i0 <= last_i_end or j0 <= last_j_end:
            continue
        anchors.append((i0, j0, run))
        last_i_end = i0 + run
        last_j_end = j0 + run
        if len(anchors) >= max_anchors:
            break
    return anchors


# ------------------------- Recursive step-by-step alignment -------------------------

def align_region(s1: str, s2: str,
                 match: int, mismatch: int, gap: int,
                 k: int, min_anchor: int, occ_limit: int,
                 max_cells: int, max_len: int,
                 band: int,
                 depth: int = 0, max_depth: int = 2):
    n, m = len(s1), len(s2)

    if n == 0:
        return "-" * m, s2, gap * m, 0
    if m == 0:
        return s1, "-" * n, gap * n, 0

    if n * m <= max_cells and max(n, m) <= max_len:
        return nw_small(s1, s2, match=match, mismatch=mismatch, gap=gap)

    if depth < max_depth:
        anchors = find_anchors(s1, s2, k=k, min_run=min_anchor, occ_limit=occ_limit)
        if anchors:
            out1, out2 = [], []
            total_score, total_matches = 0, 0
            p1 = 0
            p2 = 0

            for (i0, j0, run) in anchors:
                a1, a2, sc, mt = align_region(
                    s1[p1:i0], s2[p2:j0],
                    match, mismatch, gap,
                    k, min_anchor, occ_limit,
                    max_cells, max_len,
                    band,
                    depth=depth + 1, max_depth=max_depth
                )
                out1.append(a1); out2.append(a2)
                total_score += sc; total_matches += mt

                block1 = s1[i0:i0 + run]
                block2 = s2[j0:j0 + run]
                out1.append(block1); out2.append(block2)
                total_score += run * match
                total_matches += sum(1 for c in block1 if c != 'N')

                p1 = i0 + run
                p2 = j0 + run

            a1, a2, sc, mt = align_region(
                s1[p1:], s2[p2:],
                match, mismatch, gap,
                k, min_anchor, occ_limit,
                max_cells, max_len,
                band,
                depth=depth + 1, max_depth=max_depth
            )
            out1.append(a1); out2.append(a2)
            total_score += sc; total_matches += mt

            return "".join(out1), "".join(out2), total_score, total_matches

    return banded_nw(s1, s2, match=match, mismatch=mismatch, gap=gap, band=band)


def align_two_genomes(s1: str, s2: str,
                      match=1, mismatch=-1, gap=-1,
                      k=13, min_anchor=30, occ_limit=40,
                      max_cells=2_000_000, max_len=2500,
                      band=800):
    if s1 == s2:
        matches = sum(1 for c in s1 if c != 'N')
        aln_len = len(s1)
        score = aln_len * match
        identity = (matches / aln_len) if aln_len else 0.0
        return s1, s2, score, matches, aln_len, identity

    a1, a2, score, matches = align_region(
        s1, s2,
        match, mismatch, gap,
        k, min_anchor, occ_limit,
        max_cells, max_len,
        band
    )
    aln_len = len(a1)
    identity = (matches / aln_len) if aln_len else 0.0
    return a1, a2, score, matches, aln_len, identity


# ------------------------- Bars + Console -------------------------

def make_midline(a1: str, a2: str) -> str:
    out = []
    for x, y in zip(a1, a2):
        if x == '-' or y == '-':
            out.append(' ')
        elif x == y and x != 'N':
            out.append('|')
        else:
            out.append('.')
    return "".join(out)


def print_alignment_console(nameA: str, nameB: str, a1: str, a2: str,
                            score: int, matches: int, identity: float,
                            width: int = 100, max_blocks: int = 10):
    mid = make_midline(a1, a2)
    print("\n" + "=" * 110)
    print(f"ALIGNMENT: {nameA}  VS  {nameB}")
    print(f"Score={score} | Matches={matches} | Length={len(a1)} | Similarity={identity*100:.2f}%")
    print("-" * 110)

    blocks = 0
    for i in range(0, len(a1), width):
        print(a1[i:i + width])
        print(mid[i:i + width])
        print(a2[i:i + width])
        print()
        blocks += 1
        if max_blocks > 0 and blocks >= max_blocks and (i + width) < len(a1):
            print("... (truncated console output; full alignment is visible in the app) ...\n")
            break


def format_alignment_blocks(a1: str, a2: str, width: int = 120) -> str:
    mid = make_midline(a1, a2)
    out = []
    for i in range(0, len(a1), width):
        out.append(a1[i:i + width])
        out.append(mid[i:i + width])
        out.append(a2[i:i + width])
        out.append("")
    return "\n".join(out)


def extract_match_segments(a1: str, a2: str, min_run: int = 50):
    segs = []
    posA = 0
    posB = 0
    run_startA = None
    run_startB = None
    run_len = 0

    for x, y in zip(a1, a2):
        advA = (x != '-')
        advB = (y != '-')
        is_match = (x != '-' and y != '-' and x == y and x != 'N')

        if is_match:
            if run_len == 0:
                run_startA = posA
                run_startB = posB
            run_len += 1
        else:
            if run_len >= min_run:
                segs.append((run_startA, run_startB, run_len))
            run_len = 0
            run_startA = None
            run_startB = None

        if advA:
            posA += 1
        if advB:
            posB += 1

    if run_len >= min_run:
        segs.append((run_startA, run_startB, run_len))

    return segs


# ------------------------- GUI -------------------------

class AlignmentApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("HIV Pairwise Alignment Application (NW step-by-step)")
        self.geometry("1400x780")

        self.genomes = []
        self.pairs = []
        self.results = {}
        self.row_to_pair = {}

        self.sim_matrix = None
        self._worker = None

        top = ttk.Frame(self, padding=8)
        top.pack(fill="x")

        self.folder_var = tk.StringVar(value="")
        ttk.Label(top, text="Folder:").pack(side="left")
        ttk.Entry(top, textvariable=self.folder_var, width=70).pack(side="left", padx=6)
        ttk.Button(top, text="Browse", command=self.browse_folder).pack(side="left")

        self.k_var = tk.IntVar(value=13)
        self.min_anchor_var = tk.IntVar(value=30)
        self.band_var = tk.IntVar(value=800)
        self.min_run_var = tk.IntVar(value=50)
        self.console_width_var = tk.IntVar(value=100)
        self.console_blocks_var = tk.IntVar(value=10)

        params = ttk.Frame(top)
        params.pack(side="left", padx=14)

        ttk.Label(params, text="k").grid(row=0, column=0, padx=4)
        ttk.Spinbox(params, from_=9, to=25, width=4, textvariable=self.k_var).grid(row=0, column=1)

        ttk.Label(params, text="minAnchor").grid(row=0, column=2, padx=4)
        ttk.Spinbox(params, from_=10, to=200, width=6, textvariable=self.min_anchor_var).grid(row=0, column=3)

        ttk.Label(params, text="band").grid(row=0, column=4, padx=4)
        ttk.Spinbox(params, from_=200, to=8000, width=7, textvariable=self.band_var).grid(row=0, column=5)

        ttk.Label(params, text="bars>=").grid(row=0, column=6, padx=4)
        ttk.Spinbox(params, from_=10, to=500, width=6, textvariable=self.min_run_var).grid(row=0, column=7)

        ttk.Label(params, text="consoleWidth").grid(row=1, column=0, padx=4, pady=(6, 0))
        ttk.Spinbox(params, from_=60, to=200, width=6, textvariable=self.console_width_var).grid(row=1, column=1, pady=(6, 0))

        ttk.Label(params, text="consoleBlocks").grid(row=1, column=2, padx=4, pady=(6, 0))
        ttk.Spinbox(params, from_=0, to=9999, width=7, textvariable=self.console_blocks_var).grid(row=1, column=3, pady=(6, 0))
        ttk.Label(params, text="(0=print ALL)").grid(row=1, column=4, columnspan=2, sticky="w", pady=(6, 0))

        ttk.Button(top, text="Run All Pairwise (45)", command=self.run_all_pairs).pack(side="right", padx=6)

        prog = ttk.Frame(self, padding=(8, 0, 8, 8))
        prog.pack(fill="x")
        self.status_var = tk.StringVar(value="Select folder with 10 HIV FASTA files.")
        ttk.Label(prog, textvariable=self.status_var).pack(side="left")
        self.pbar = ttk.Progressbar(prog, length=420, mode="determinate")
        self.pbar.pack(side="right")

        main = ttk.Frame(self, padding=8)
        main.pack(fill="both", expand=True)

        # LEFT: pair list table
        left = ttk.Frame(main)
        left.pack(side="left", fill="y")

        ttk.Label(left, text="Similarity Table (click row to view)").pack(anchor="w")

        cols = ("pair", "similarity", "score", "matches", "alnlen")
        self.tree = ttk.Treeview(left, columns=cols, show="headings", height=28)
        for c, t, w in [
            ("pair", "Pair", 360),
            ("similarity", "Similarity %", 90),
            ("score", "NW Score", 80),
            ("matches", "Matches", 80),
            ("alnlen", "Align Len", 80),
        ]:
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor=("w" if c == "pair" else "center"))

        ysb = ttk.Scrollbar(left, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=ysb.set)
        self.tree.pack(side="left", fill="y")
        ysb.pack(side="right", fill="y")
        self.tree.bind("<<TreeviewSelect>>", self.on_table_select)

        # RIGHT: bars + blocks + alignment
        right = ttk.Frame(main)
        right.pack(side="left", fill="both", expand=True, padx=(12, 0))

        bars_frame = ttk.LabelFrame(right, text="Bar Visualization (similar regions)")
        bars_frame.pack(fill="x")
        self.canvas = tk.Canvas(bars_frame, height=230, bg="white")
        self.canvas.pack(fill="x")

        blocks_frame = ttk.LabelFrame(right, text="Match Blocks")
        blocks_frame.pack(fill="x", pady=(10, 0))

        self.blocks_list = tk.Listbox(blocks_frame, height=7)
        self.blocks_list.pack(fill="x", expand=False, padx=6, pady=6)
        self.blocks_list.bind("<<ListboxSelect>>", self.on_block_select)

        self.block_info = tk.StringVar(value="")
        ttk.Label(blocks_frame, textvariable=self.block_info).pack(anchor="w", padx=8, pady=(0, 8))

        text_frame = ttk.LabelFrame(right, text="Alignment (full shown here; console prints blocks)")
        text_frame.pack(fill="both", expand=True, pady=(10, 0))

        self.text = tk.Text(text_frame, wrap="none")
        self.text.pack(fill="both", expand=True)

        xscroll = ttk.Scrollbar(text_frame, orient="horizontal", command=self.text.xview)
        yscroll = ttk.Scrollbar(text_frame, orient="vertical", command=self.text.yview)
        self.text.configure(xscrollcommand=xscroll.set, yscrollcommand=yscroll.set)
        xscroll.pack(side="bottom", fill="x")
        yscroll.pack(side="right", fill="y")

        # Button to show matrix after finishing (also auto pops up)
        ttk.Button(left, text="Show Similarity Matrix (10×10)", command=self.show_matrix_window).pack(
            fill="x", pady=(8, 0)
        )

    def browse_folder(self):
        folder = filedialog.askdirectory(title="Select folder with HIV FASTA files")
        if not folder:
            return
        self.folder_var.set(folder)
        self.load_genomes(folder)

    def load_genomes(self, folder: str):
        files = sorted([f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))])
        genomes = []
        for fn in files:
            seq = read_fasta_first_record(os.path.join(folder, fn))
            if seq:
                genomes.append((fn, seq))

        if len(genomes) < 2:
            messagebox.showerror("Error", "Need at least 2 valid FASTA files.")
            return

        self.genomes = genomes[:10]
        n = len(self.genomes)
        self.sim_matrix = [[None for _ in range(n)] for _ in range(n)]
        for i in range(n):
            self.sim_matrix[i][i] = 100.0

        self.results.clear()
        self.pairs.clear()
        self.row_to_pair.clear()

        for item in self.tree.get_children():
            self.tree.delete(item)

        row_id = 0
        for i, j in itertools.combinations(range(n), 2):
            self.pairs.append((i, j))
            nameA = self.genomes[i][0]
            nameB = self.genomes[j][0]
            iid = str(row_id)
            self.row_to_pair[iid] = (i, j)
            self.tree.insert("", "end", iid=iid, values=(f"{nameA} VS {nameB}", "-", "-", "-", "-"))
            row_id += 1

        self.pbar["value"] = 0
        self.pbar["maximum"] = len(self.pairs)
        self.status_var.set(f"Loaded {n} genomes. Ready to run 45 alignments.")

    def run_all_pairs(self):
        if not self.genomes:
            messagebox.showwarning("No genomes", "Load folder first.")
            return
        if self._worker and self._worker.is_alive():
            messagebox.showinfo("Running", "Already running.")
            return

        self.results.clear()
        self.pbar["value"] = 0
        self.status_var.set("Running alignments...")

        params = dict(
            k=int(self.k_var.get()),
            min_anchor=int(self.min_anchor_var.get()),
            band=int(self.band_var.get()),
            min_run=int(self.min_run_var.get()),
            console_width=int(self.console_width_var.get()),
            console_blocks=int(self.console_blocks_var.get()),
        )

        self._worker = threading.Thread(target=self._run_worker, args=(params,), daemon=True)
        self._worker.start()

    def _run_worker(self, params):
        total = len(self.pairs)
        for done, (i, j) in enumerate(self.pairs, start=1):
            nameA, seqA = self.genomes[i]
            nameB, seqB = self.genomes[j]

            alA, alB, score, matches, aln_len, ident = align_two_genomes(
                seqA, seqB,
                match=1, mismatch=-1, gap=-1,
                k=params["k"],
                min_anchor=params["min_anchor"],
                occ_limit=40,
                max_cells=2_000_000,
                max_len=2500,
                band=params["band"],
            )

            print_alignment_console(
                nameA, nameB, alA, alB, score, matches, ident,
                width=params["console_width"],
                max_blocks=params["console_blocks"],
            )

            segs = extract_match_segments(alA, alB, min_run=params["min_run"])
            segs_with_seq = []
            for (sA, sB, L) in segs:
                if 0 <= sA < len(seqA) and sA + L <= len(seqA):
                    segs_with_seq.append((sA, sB, L, seqA[sA:sA + L]))

            self.results[(i, j)] = {
                "nameA": nameA, "nameB": nameB,
                "lenA": len(seqA), "lenB": len(seqB),
                "alA": alA, "alB": alB,
                "score": score, "matches": matches,
                "aln_len": aln_len, "ident": ident,
                "segments": segs_with_seq,
                "min_run": params["min_run"],
            }

            sim_percent = ident * 100.0
            self.sim_matrix[i][j] = sim_percent
            self.sim_matrix[j][i] = sim_percent

            self.after(0, self._update_progress, done, total, i, j, score, matches, aln_len, ident)

        self.after(0, self._done_running)

    def _update_progress(self, done, total, i, j, score, matches, aln_len, ident):
        self.pbar["value"] = done
        nameA = self.genomes[i][0]
        nameB = self.genomes[j][0]
        self.status_var.set(f"Done {done}/{total}: {nameA} vs {nameB}")

        iid = str(done - 1)
        self.tree.item(iid, values=(
            f"{nameA} VS {nameB}",
            f"{ident*100:.2f}",
            str(score),
            str(matches),
            str(aln_len),
        ))

        self.tree.selection_set(iid)
        self.tree.see(iid)
        self.show_pair(i, j)

    def _done_running(self):
        self.status_var.set("Finished all 45 alignments.")
        # Auto pop-up the matrix table when finished
        self.show_matrix_window()
        messagebox.showinfo("Done", "Finished all pairwise alignments.\nMatrix table opened.")

    def on_table_select(self, event):
        sel = self.tree.selection()
        if not sel:
            return
        iid = sel[0]
        pair = self.row_to_pair.get(iid)
        if not pair:
            return
        self.show_pair(pair[0], pair[1])

    def show_pair(self, i: int, j: int):
        res = self.results.get((i, j))
        self.blocks_list.delete(0, tk.END)
        self.block_info.set("")
        self.canvas.delete("all")
        self.text.delete("1.0", tk.END)

        if not res:
            self.text.insert(tk.END, "Not computed yet.")
            return

        header = (
            f"{res['nameA']} VS {res['nameB']}\n"
            f"NW Score={res['score']} | Matches={res['matches']} | AlignmentLen={res['aln_len']} | "
            f"Similarity={res['ident']*100:.2f}%\n\n"
        )
        self.text.insert(tk.END, header)
        self.text.insert(tk.END, format_alignment_blocks(res["alA"], res["alB"], width=120))

        self.draw_bars(res)

        for k, (sA, sB, L, _dna) in enumerate(res["segments"], start=1):
            self.blocks_list.insert(tk.END, f"Block {k}: A[{sA}..{sA+L-1}]  B[{sB}..{sB+L-1}]  len={L}")

        self.block_info.set("Click a block above to see DNA preview." if res["segments"]
                            else "No match blocks with current bars>=. Try lower bars>=.")

    def on_block_select(self, event):
        if not self.blocks_list.curselection():
            return
        # find selected row
        sel = self.tree.selection()
        if not sel:
            return
        pair = self.row_to_pair.get(sel[0])
        if not pair:
            return
        res = self.results.get(pair)
        if not res:
            return

        bidx = int(self.blocks_list.curselection()[0])
        sA, sB, L, dna = res["segments"][bidx]
        preview = dna if L <= 160 else (dna[:160] + "...")
        self.block_info.set(
            f"Selected Block {bidx+1}: A[{sA}..{sA+L-1}]  B[{sB}..{sB+L-1}]  len={L} | DNA: {preview}"
        )

    def draw_bars(self, res):
        self.canvas.delete("all")
        self.update_idletasks()

        W = max(950, self.canvas.winfo_width())
        H = 230
        margin = 70
        track_w = W - 2 * margin
        track_h = 14
        yA = 80
        yB = 160

        lenA = res["lenA"]
        lenB = res["lenB"]
        max_len = max(lenA, lenB)

        def sx(pos):
            return margin + (pos / max_len) * track_w

        def sw(length):
            return max(1.0, (length / max_len) * track_w)

        title = (f"{res['nameA']} vs {res['nameB']} | Similarity={res['ident']*100:.2f}% | "
                 f"bars>= {res['min_run']} (count={len(res['segments'])})")
        self.canvas.create_text(10, 10, anchor="nw", text=title, font=("Segoe UI", 10, "bold"))

        self.canvas.create_text(margin, yA - 22, anchor="w", text=f"{res['nameA']} (len={lenA})")
        self.canvas.create_rectangle(margin, yA, margin + track_w, yA + track_h, outline="black", fill="#f0f0f0")

        self.canvas.create_text(margin, yB - 22, anchor="w", text=f"{res['nameB']} (len={lenB})")
        self.canvas.create_rectangle(margin, yB, margin + track_w, yB + track_h, outline="black", fill="#f0f0f0")

        for (sA, sB, L, _dna) in res["segments"]:
            xA = sx(sA); w = sw(L)
            xB = sx(sB)
            self.canvas.create_rectangle(xA, yA, xA + w, yA + track_h, outline="", fill="#2E86C1")
            self.canvas.create_rectangle(xB, yB, xB + w, yB + track_h, outline="", fill="#28B463")

    # --------- MATRIX WINDOW AFTER FINISH ---------

    def show_matrix_window(self):
        if not self.genomes or self.sim_matrix is None:
            messagebox.showwarning("Matrix", "No matrix yet. Run alignments first.")
            return

        win = tk.Toplevel(self)
        win.title("Similarity Matrix (10×10) - % Identity")
        win.geometry("900x500")

        names = [fn for (fn, _s) in self.genomes]
        n = len(names)

        table = ttk.Treeview(win, show="headings")
        cols = ["Genome"] + [str(i + 1) for i in range(n)]
        table["columns"] = cols

        table.heading("Genome", text="Genome")
        table.column("Genome", width=240, anchor="w")

        for j in range(n):
            table.heading(str(j + 1), text=str(j + 1))
            table.column(str(j + 1), width=60, anchor="center")

        for i in range(n):
            row = [f"{i+1}: {names[i]}"]
            for j in range(n):
                v = self.sim_matrix[i][j]
                row.append("-" if v is None else f"{v:.2f}")
            table.insert("", "end", values=row)

        ysb = ttk.Scrollbar(win, orient="vertical", command=table.yview)
        xsb = ttk.Scrollbar(win, orient="horizontal", command=table.xview)
        table.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)

        table.pack(side="left", fill="both", expand=True)
        ysb.pack(side="right", fill="y")
        xsb.pack(side="bottom", fill="x")

        # optional heatmap button if matplotlib installed
        def open_heatmap():
            try:
                import matplotlib.pyplot as plt
            except Exception:
                messagebox.showerror("Heatmap", "matplotlib not installed. Install it or skip heatmap.")
                return

            import numpy as np
            data = np.zeros((n, n), dtype=float)
            for i in range(n):
                for j in range(n):
                    v = self.sim_matrix[i][j]
                    data[i, j] = 0.0 if v is None else v

            plt.figure(figsize=(8, 6))
            im = plt.imshow(data, aspect='auto')
            plt.title("Pairwise Genome Similarity Heatmap")
            plt.xticks(range(n), [f"s{j+1}" for j in range(n)], rotation=45, ha="right")
            plt.yticks(range(n), [f"s{i+1}" for i in range(n)])
            for i in range(n):
                for j in range(n):
                    plt.text(j, i, f"{data[i,j]:.2f}", ha="center", va="center", fontsize=7)
            plt.colorbar(im, label="Similarity (%)")
            plt.tight_layout()
            plt.show()

        ttk.Button(win, text="Open Heatmap (optional)", command=open_heatmap).pack(pady=8)


def main():
    app = AlignmentApp()
    app.mainloop()


if __name__ == "__main__":
    main()
