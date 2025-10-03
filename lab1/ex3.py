#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import time
import threading
import gzip
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

class Canceled(Exception):
    pass

# ---------- STREAMED FASTA ANALYZER WITH PROGRESS CALLBACK ----------
def analyze_fasta(path, ignore_case=True, progress_cb=None, cancel_event=None):
    """
    Stream through a FASTA file (supports .gz) and compute:
    - first header,
    - number of sequences,
    - total symbols,
    - alphabet (first-seen order, no sets),
    - counts and relative frequencies per symbol.

    progress_cb(bytes_read, bytes_total, records, symbols) is called periodically.
    cancel_event (threading.Event) aborts cleanly if set.
    """
    opener = gzip.open if path.lower().endswith(".gz") else open
    is_gz = path.lower().endswith(".gz")
    bytes_total = os.path.getsize(path)

    counts = {}   # insertion-ordered dict (preserves first-seen order)
    total_symbols = 0
    first_header = None
    n_records = 0

    def add_char(ch):
        nonlocal total_symbols
        if ignore_case:
            ch = ch.upper()
        if ch not in counts:
            counts[ch] = 0
        counts[ch] += 1
        total_symbols += 1

    def bytes_read_now(fh):
        # Try to get an accurate byte position for progress.
        try:
            if is_gz:
                # Text wrapper -> BufferedReader -> GzipFile -> underlying fileobj
                return fh.buffer.raw.fileobj.tell()
            else:
                # For plain text, the buffered binary stream position is in bytes.
                return fh.buffer.tell()
        except Exception:
            return None

    update_interval = 0.2  # seconds
    next_update = time.time() + update_interval

    with opener(path, mode="rt", encoding="utf-8", errors="ignore", newline="") as fh:
        in_seq = False
        for line in fh:
            if cancel_event and cancel_event.is_set():
                raise Canceled("Canceled by user")

            if line.startswith(">"):
                n_records += 1
                if first_header is None:
                    first_header = line.strip()[1:]
                in_seq = True
            elif in_seq:
                for ch in line.strip():
                    if ch.isalpha() or ch == "*":
                        add_char(ch)

            # periodic progress update
            if progress_cb and time.time() >= next_update:
                br = bytes_read_now(fh)
                if br is None:
                    br = 0
                progress_cb(br, bytes_total, n_records, total_symbols)
                next_update += update_interval

    # final progress update (100%)
    if progress_cb:
        br = bytes_total
        progress_cb(br, bytes_total, n_records, total_symbols)

    alphabet = list(counts.keys())
    freqs = {sym: counts[sym] / total_symbols for sym in counts} if total_symbols > 0 else {}

    return {
        "file": os.path.basename(path),
        "size_bytes": bytes_total,
        "first_header": first_header or "(none found)",
        "records": n_records,
        "length": total_symbols,
        "alphabet": alphabet,
        "counts": counts,
        "freqs": freqs,
        "progress_basis": "compressed bytes" if is_gz else "file bytes",
    }

# ---------- TKINTER GUI ----------
class FastaGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Analyzer (detailed progress)")
        self.geometry("820x620")
        self.minsize(760, 560)

        self.app_dir = os.path.dirname(os.path.abspath(__file__))  # lab1 folder
        self.cancel_event = None
        self.start_time = None

        # Style: make the progress bar THICK and visible
        style = ttk.Style(self)
        style.theme_use(style.theme_use())
        style.configure(
            "Big.Horizontal.TProgressbar",
            thickness=24,
        )

        # ---- Top controls
        top = ttk.Frame(self, padding=12)
        top.pack(fill="x")

        self.ignore_case = tk.BooleanVar(value=True)

        self.pick_btn = ttk.Button(top, text="Choose FASTA…", command=self.choose_file)
        self.pick_btn.pack(side="left")

        ttk.Checkbutton(top, text="Ignore case (a/A same)",
                        variable=self.ignore_case).pack(side="left", padx=(12, 0))

        self.path_label = ttk.Label(top, text="No file chosen")
        self.path_label.pack(side="left", padx=12)

        # ---- Progress area (big + detailed)
        prog = ttk.Frame(self, padding=(12, 0, 12, 0))
        prog.pack(fill="x")

        self.prog_title = ttk.Label(prog, text="Progress", font=("Segoe UI", 12, "bold"))
        self.prog_title.pack(anchor="w", pady=(4, 2))

        self.progress = ttk.Progressbar(
            prog, mode="determinate", style="Big.Horizontal.TProgressbar", length=780, maximum=100
        )
        self.progress.pack(fill="x", pady=(0, 6))

        self.prog_detail = ttk.Label(
            prog,
            text="Waiting…",
            font=("Consolas", 10),
            justify="left"
        )
        self.prog_detail.pack(anchor="w", pady=(0, 8))

        btns = ttk.Frame(prog)
        btns.pack(fill="x", pady=(0, 10))
        self.abort_btn = ttk.Button(btns, text="Abort", command=self.abort, state="disabled")
        self.abort_btn.pack(side="right")

        # ---- Output
        out_frame = ttk.Frame(self, padding=12)
        out_frame.pack(fill="both", expand=True)

        self.output = tk.Text(out_frame, wrap="word", height=20)
        self.output.pack(fill="both", expand=True)
        self.output.configure(state="disabled")

        # Slight DPI bump (if supported)
        try:
            self.tk.call("tk", "scaling", 1.1)
        except tk.TclError:
            pass

    # ---- File picking
    def choose_file(self):
        path = filedialog.askopenfilename(
            title="Select a FASTA file",
            initialdir=self.app_dir,  # starts in your lab1 folder
            filetypes=[
                ("FASTA files", "*.fa *.fna *.fasta *.faa *.fas *.fa.gz *.fna.gz *.fasta.gz"),
                ("All files", "*.*"),
            ],
        )
        if not path:
            return

        self.path_label.config(text=path)
        self.prepare_progress("Opening file…")
        self.pick_btn.config(state="disabled")
        self.abort_btn.config(state="normal")
        self.cancel_event = threading.Event()
        self.start_time = time.time()

        t = threading.Thread(
            target=self._analyze_worker, args=(path, self.ignore_case.get()), daemon=True
        )
        t.start()

    def prepare_progress(self, msg):
        self.progress.config(mode="determinate", maximum=100, value=0)
        self.prog_detail.config(text=msg)

    def abort(self):
        if self.cancel_event:
            self.cancel_event.set()
            self.abort_btn.config(state="disabled")
            self.prog_detail.config(text="Canceling…")

    # ---- Background worker
    def _analyze_worker(self, path, ignore_case):
        def ui_progress(bytes_read, bytes_total, records, symbols):
            # compute stats
            now = time.time()
            elapsed = max(1e-6, now - (self.start_time or now))
            mb_read = bytes_read / (1024 * 1024)
            mb_total = bytes_total / (1024 * 1024) if bytes_total else 0.0
            speed = mb_read / elapsed
            pct = 0.0 if not bytes_total else min(100.0, (bytes_read / bytes_total) * 100.0)

            detail = (
                f"{pct:6.2f}%  |  "
                f"{mb_read:,.2f} MB / {mb_total:,.2f} MB  |  "
                f"{speed:,.2f} MB/s  |  "
                f"records: {records:,}  |  symbols: {symbols:,}"
            )
            # send to main thread
            self.after(0, lambda: self.update_progress(pct, detail))

        try:
            results = analyze_fasta(
                path,
                ignore_case=ignore_case,
                progress_cb=ui_progress,
                cancel_event=self.cancel_event,
            )
            self.after(0, lambda: self.show_results(results))
        except Canceled:
            self.after(0, self.show_canceled)
        except Exception as e:
            self.after(0, lambda: self.show_error(e))

    # ---- UI updates (main thread)
    def update_progress(self, pct, detail_text):
        # Make the bar super visible and update text
        self.progress["value"] = pct
        self.prog_detail.config(text=detail_text)

    def show_canceled(self):
        self.pick_btn.config(state="normal")
        self.abort_btn.config(state="disabled")
        self.progress["value"] = 0
        self.prog_detail.config(text="Canceled.")
        # keep whatever was in the output area

    def show_results(self, res):
        self.pick_btn.config(state="normal")
        self.abort_btn.config(state="disabled")
        self.progress["value"] = 100

        size_str = f"{res['size_bytes'] / (1024*1024):.2f} MB" if res.get("size_bytes") else "n/a"
        alpha_str = " ".join(res["alphabet"]) if res["alphabet"] else "(empty)"

        lines = []
        lines.append(f"File: {res['file']}  |  Size: {size_str}  |  Progress basis: {res['progress_basis']}")
        lines.append(f"First header: {res['first_header']}")
        lines.append(f"Sequences (records): {res['records']}")
        lines.append(f"Total symbols counted: {res['length']}")
        lines.append(f"Alphabet (first-seen order): {alpha_str}")
        lines.append("")
        lines.append("Symbol\tCount\tRelative frequency")
        lines.append("------------------------------------")
        for sym, cnt in res["counts"].items():
            freq = res["freqs"].get(sym, 0.0)
            lines.append(f"{sym}\t{cnt}\t{freq:.6f}")

        self.output_set("\n".join(lines))

    def show_error(self, err):
        self.pick_btn.config(state="normal")
        self.abort_btn.config(state="disabled")
        messagebox.showerror("Error", f"Failed to read FASTA:\n\n{err}")

    def output_set(self, text):
        self.output.configure(state="normal")
        self.output.delete("1.0", "end")
        self.output.insert("end", text)
        self.output.configure(state="disabled")

if __name__ == "__main__":
    app = FastaGUI()
    app.mainloop()
