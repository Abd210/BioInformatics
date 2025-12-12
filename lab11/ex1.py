

import tkinter as tk
from tkinter import ttk, messagebox

DEFAULT_S1 = "ACCGTGAAGCCAATAC"
DEFAULT_S2 = "AGCGTGCAGCCAATAC"


def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=0):
    n, m = len(s1), len(s2)

    score = [[0]*(m+1) for _ in range(n+1)]
    trace = [[None]*(m+1) for _ in range(n+1)]  # 'D', 'U', 'L'

    trace[0][0] = 'E'
    for i in range(1, n+1):
        score[i][0] = score[i-1][0] + gap
        trace[i][0] = 'U'
    for j in range(1, m+1):
        score[0][j] = score[0][j-1] + gap
        trace[0][j] = 'L'

    def sub(a, b):
        return match if a == b else mismatch

    # Fill
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = score[i-1][j-1] + sub(s1[i-1], s2[j-1])
            up   = score[i-1][j]   + gap
            left = score[i][j-1]   + gap

            best = diag
            d = 'D'
            # tie-break: D > U > L (like many demos)
            if up > best:
                best = up
                d = 'U'
            if left > best:
                best = left
                d = 'L'

            score[i][j] = best
            trace[i][j] = d

    # Traceback
    i, j = n, m
    a1, a2, mid = [], [], []
    matches = 0
    path = [(i, j)]
    while not (i == 0 and j == 0):
        d = trace[i][j]
        if d == 'D':
            c1, c2 = s1[i-1], s2[j-1]
            a1.append(c1); a2.append(c2)
            if c1 == c2:
                mid.append('|'); matches += 1
            else:
                mid.append('.')
            i -= 1; j -= 1
        elif d == 'U':
            c1 = s1[i-1]
            a1.append(c1); a2.append('-'); mid.append(' ')
            i -= 1
        elif d == 'L':
            c2 = s2[j-1]
            a1.append('-'); a2.append(c2); mid.append(' ')
            j -= 1
        else:
            raise RuntimeError("Bad traceback direction.")
        path.append((i, j))

    a1 = ''.join(reversed(a1))
    a2 = ''.join(reversed(a2))
    mid = ''.join(reversed(mid))
    path.reverse()

    aln_len = len(a1)
    identity = (matches / aln_len) if aln_len else 0.0
    final_score = score[n][m]

    return score, trace, a1, mid, a2, final_score, matches, aln_len, identity, path


# --------------------- Drawing Helpers ---------------------

def clamp(x, lo, hi):
    return lo if x < lo else hi if x > hi else x

def lerp(a, b, t):
    return int(a + (b - a) * t)

def rgb_hex(r, g, b):
    return f"#{r:02x}{g:02x}{b:02x}"

def heat_color(val, vmin, vmax):
    """
    Simple gradient: dark blue -> purple -> red
    """
    if vmax == vmin:
        t = 0.5
    else:
        t = (val - vmin) / (vmax - vmin)
    t = clamp(t, 0.0, 1.0)

    # two-stage gradient for nicer look
    # 0..0.5: dark blue -> purple
    # 0.5..1: purple -> red
    if t <= 0.5:
        tt = t / 0.5
        r = lerp(10, 120, tt)
        g = lerp(10, 0, tt)
        b = lerp(60, 160, tt)
    else:
        tt = (t - 0.5) / 0.5
        r = lerp(120, 230, tt)
        g = lerp(0, 30, tt)
        b = lerp(160, 40, tt)

    return rgb_hex(r, g, b)

def draw_heatmap(canvas, matrix, show_grid=True):
    canvas.delete("all")
    rows = len(matrix)
    cols = len(matrix[0]) if rows else 0
    if rows == 0 or cols == 0:
        return

    w = int(canvas.winfo_width())
    h = int(canvas.winfo_height())
    if w <= 1 or h <= 1:
        w, h = 360, 280

    cell_w = w / cols
    cell_h = h / rows

    vmin = min(min(r) for r in matrix)
    vmax = max(max(r) for r in matrix)

    for i in range(rows):
        for j in range(cols):
            x0 = j * cell_w
            y0 = i * cell_h
            x1 = x0 + cell_w
            y1 = y0 + cell_h
            fill = heat_color(matrix[i][j], vmin, vmax)
            if show_grid:
                canvas.create_rectangle(x0, y0, x1, y1, fill=fill, outline="#202020", width=1)
            else:
                canvas.create_rectangle(x0, y0, x1, y1, fill=fill, outline=fill)

def draw_path_grid(canvas, path_coords, rows, cols, show_grid=True):
    """
    rows = n+1, cols = m+1 for DP grid
    """
    canvas.delete("all")
    w = int(canvas.winfo_width())
    h = int(canvas.winfo_height())
    if w <= 1 or h <= 1:
        w, h = 360, 280

    cell_w = w / cols
    cell_h = h / rows

    path_set = set(path_coords)

    bg = "#fff7b3"   # pale yellow
    pathc = "#cc2b2b"  # red blocks

    for i in range(rows):
        for j in range(cols):
            x0 = j * cell_w
            y0 = i * cell_h
            x1 = x0 + cell_w
            y1 = y0 + cell_h
            fill = pathc if (i, j) in path_set else bg
            if show_grid:
                canvas.create_rectangle(x0, y0, x1, y1, fill=fill, outline="#202020", width=1)
            else:
                canvas.create_rectangle(x0, y0, x1, y1, fill=fill, outline=fill)


# --------------------- GUI App ---------------------

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Needleman–Wunsch Alignment (Native Python)")
        self.geometry("1100x650")
        self.minsize(980, 600)

        self._build_ui()

        # run initial alignment
        self.after(50, self.align)

    def _build_ui(self):
        outer = ttk.Frame(self, padding=8)
        outer.pack(fill="both", expand=True)

        outer.columnconfigure(0, weight=0)
        outer.columnconfigure(1, weight=1)
        outer.rowconfigure(0, weight=1)

        # Left panel
        left = ttk.Frame(outer)
        left.grid(row=0, column=0, sticky="ns", padx=(0, 10))
        left.columnconfigure(0, weight=1)

        seq_frame = ttk.LabelFrame(left, text="Sequences", padding=8)
        seq_frame.grid(row=0, column=0, sticky="ew", pady=(0, 10))

        ttk.Label(seq_frame, text="Sq 1 =").grid(row=0, column=0, sticky="w")
        self.s1_var = tk.StringVar(value=DEFAULT_S1)
        self.s1_entry = ttk.Entry(seq_frame, textvariable=self.s1_var, width=28)
        self.s1_entry.grid(row=0, column=1, sticky="ew", padx=(6, 0))

        ttk.Label(seq_frame, text="Sq 2 =").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.s2_var = tk.StringVar(value=DEFAULT_S2)
        self.s2_entry = ttk.Entry(seq_frame, textvariable=self.s2_var, width=28)
        self.s2_entry.grid(row=1, column=1, sticky="ew", padx=(6, 0), pady=(6, 0))

        seq_frame.columnconfigure(1, weight=1)

        param_frame = ttk.LabelFrame(left, text="Parameters", padding=8)
        param_frame.grid(row=1, column=0, sticky="ew", pady=(0, 10))

        ttk.Label(param_frame, text="Gap =").grid(row=0, column=0, sticky="w")
        ttk.Label(param_frame, text="Match =").grid(row=1, column=0, sticky="w", pady=(6,0))
        ttk.Label(param_frame, text="Mismatch =").grid(row=2, column=0, sticky="w", pady=(6,0))

        self.gap_var = tk.IntVar(value=0)
        self.match_var = tk.IntVar(value=1)
        self.mismatch_var = tk.IntVar(value=-1)

        self.gap_spin = ttk.Spinbox(param_frame, from_=-10, to=10, textvariable=self.gap_var, width=6)
        self.match_spin = ttk.Spinbox(param_frame, from_=-10, to=10, textvariable=self.match_var, width=6)
        self.mismatch_spin = ttk.Spinbox(param_frame, from_=-10, to=10, textvariable=self.mismatch_var, width=6)

        self.gap_spin.grid(row=0, column=1, sticky="w", padx=(6, 0))
        self.match_spin.grid(row=1, column=1, sticky="w", padx=(6, 0), pady=(6,0))
        self.mismatch_spin.grid(row=2, column=1, sticky="w", padx=(6, 0), pady=(6,0))

        opt_frame = ttk.LabelFrame(left, text="Options", padding=8)
        opt_frame.grid(row=2, column=0, sticky="ew", pady=(0, 10))

        self.grid_var = tk.BooleanVar(value=True)
        self.auto_align_var = tk.BooleanVar(value=False)

        ttk.Checkbutton(opt_frame, text="Plot grid", variable=self.grid_var, command=self.redraw).grid(row=0, column=0, sticky="w")
        ttk.Checkbutton(opt_frame, text="Auto align on typing", variable=self.auto_align_var).grid(row=1, column=0, sticky="w", pady=(6,0))

        self.align_btn = ttk.Button(left, text="Align", command=self.align)
        self.align_btn.grid(row=3, column=0, sticky="ew", pady=(0, 10))

        presets = ttk.LabelFrame(left, text="Presets", padding=8)
        presets.grid(row=4, column=0, sticky="ew")
        presets.columnconfigure(0, weight=1)
        presets.columnconfigure(1, weight=1)

        # 4 simple preset scoring sets (keep it clean)
        ttk.Button(presets, text="Setting 1\n(0,1,-1)", command=lambda: self.set_params(0, 1, -1)).grid(row=0, column=0, sticky="ew", padx=(0,6), pady=(0,6))
        ttk.Button(presets, text="Setting 2\n(-1,1,-1)", command=lambda: self.set_params(-1, 1, -1)).grid(row=0, column=1, sticky="ew", pady=(0,6))
        ttk.Button(presets, text="Setting 3\n(-2,2,-1)", command=lambda: self.set_params(-1, 2, -2)).grid(row=1, column=0, sticky="ew", padx=(0,6))
        ttk.Button(presets, text="Setting 4\n(-2,1,-1)", command=lambda: self.set_params(-1, 1, -2)).grid(row=1, column=1, sticky="ew")

        # Right panel (plots + output)
        right = ttk.Frame(outer)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)
        right.columnconfigure(1, weight=1)

        heat_frame = ttk.LabelFrame(right, text="Graphic representation of the alignment matrix", padding=6)
        heat_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 8), pady=(0, 8))
        heat_frame.rowconfigure(0, weight=1)
        heat_frame.columnconfigure(0, weight=1)

        path_frame = ttk.LabelFrame(right, text="Traceback path (cells on chosen path)", padding=6)
        path_frame.grid(row=0, column=1, sticky="nsew", pady=(0, 8))
        path_frame.rowconfigure(0, weight=1)
        path_frame.columnconfigure(0, weight=1)

        self.heat_canvas = tk.Canvas(heat_frame, bg="black", highlightthickness=0)
        self.heat_canvas.grid(row=0, column=0, sticky="nsew")

        self.path_canvas = tk.Canvas(path_frame, bg="white", highlightthickness=0)
        self.path_canvas.grid(row=0, column=0, sticky="nsew")

        out_frame = ttk.LabelFrame(right, text="Show Alignment", padding=6)
        out_frame.grid(row=1, column=0, columnspan=2, sticky="nsew")
        out_frame.rowconfigure(0, weight=1)
        out_frame.columnconfigure(0, weight=1)

        self.out_text = tk.Text(out_frame, wrap="none", height=12)
        self.out_text.grid(row=0, column=0, sticky="nsew")

        yscroll = ttk.Scrollbar(out_frame, orient="vertical", command=self.out_text.yview)
        yscroll.grid(row=0, column=1, sticky="ns")
        self.out_text.configure(yscrollcommand=yscroll.set)

        xscroll = ttk.Scrollbar(out_frame, orient="horizontal", command=self.out_text.xview)
        xscroll.grid(row=1, column=0, sticky="ew")
        self.out_text.configure(xscrollcommand=xscroll.set)

        # Re-render on resize
        self.heat_canvas.bind("<Configure>", lambda e: self.redraw())
        self.path_canvas.bind("<Configure>", lambda e: self.redraw())

        # optional auto-align
        for ent in (self.s1_entry, self.s2_entry):
            ent.bind("<KeyRelease>", self._maybe_auto_align)

        for spin in (self.gap_spin, self.match_spin, self.mismatch_spin):
            spin.bind("<KeyRelease>", self._maybe_auto_align)
            spin.bind("<<Increment>>", self._maybe_auto_align)
            spin.bind("<<Decrement>>", self._maybe_auto_align)

        self._last_result = None

    def _maybe_auto_align(self, _event=None):
        if self.auto_align_var.get():
            self.align()

    def set_params(self, gap, match, mismatch):
        self.gap_var.set(gap)
        self.match_var.set(match)
        self.mismatch_var.set(mismatch)
        self.align()

    def align(self):
        s1 = self.s1_var.get().strip().upper()
        s2 = self.s2_var.get().strip().upper()

        if not s1 or not s2:
            messagebox.showwarning("Input error", "Please enter both sequences.")
            return

        try:
            gap = int(self.gap_var.get())
            match = int(self.match_var.get())
            mismatch = int(self.mismatch_var.get())
        except Exception:
            messagebox.showwarning("Input error", "Parameters must be integers.")
            return

        try:
            result = needleman_wunsch(s1, s2, match=match, mismatch=mismatch, gap=gap)
        except Exception as e:
            messagebox.showerror("Alignment error", str(e))
            return

        self._last_result = (s1, s2, gap, match, mismatch, result)
        self.redraw()
        self._render_output()

    def redraw(self):
        if not self._last_result:
            return
        s1, s2, gap, match, mismatch, result = self._last_result
        score, trace, a1, mid, a2, final_score, matches, aln_len, identity, path = result

        show_grid = self.grid_var.get()
        draw_heatmap(self.heat_canvas, score, show_grid=show_grid)
        draw_path_grid(self.path_canvas, path, rows=len(score), cols=len(score[0]), show_grid=show_grid)

    def _render_output(self):
        s1, s2, gap, match, mismatch, result = self._last_result
        score, trace, a1, mid, a2, final_score, matches, aln_len, identity, path = result

        self.out_text.delete("1.0", "end")

        self.out_text.insert("end", f"{a1}\n")
        self.out_text.insert("end", f"{mid}\n")
        self.out_text.insert("end", f"{a2}\n\n")

        self.out_text.insert("end", f"Matches = {matches}\n")
        self.out_text.insert("end", f"Length  = {aln_len}\n")
        self.out_text.insert("end", f"Similarity = {identity*100:.0f} %\n")
        self.out_text.insert("end", f"Final score = {final_score}\n")
        self.out_text.insert("end", f"Tracing back: M[{len(s1)},{len(s2)}]\n")

        # Optional: show a few path coords
        if len(path) <= 80:
            self.out_text.insert("end", f"Path cells: {path}\n")


if __name__ == "__main__":
    app = App()
    app.mainloop()
