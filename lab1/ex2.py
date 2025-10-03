def relative_frequencies(seq: str, ignore_case: bool = True, strip_ws: bool = True):

    if strip_ws:
        seq = "".join(seq.split())
    if not seq:
        return {}
    if ignore_case:
        seq = seq.upper()

    counts = {}  
    for ch in seq:
        if ch not in counts:
            counts[ch] = 0
        counts[ch] += 1

    n = len(seq)
    return {ch: counts[ch] / n for ch in counts}



S = "ATTTCGCCGATA"
freqs = relative_frequencies(S)
for ch, f in freqs.items():
    print(f"{ch}: {f:.3f} ({f*100:.1f}%)")
