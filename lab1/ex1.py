def alphabet(seq: str):
    seq = seq.strip()
    result = []
    for ch in seq:
        if ch not in result:
            result.append(ch)
    return result

# demo
s = "ATTTCGCCGATA"
print(alphabet(s))  # -> ['A', 'T', 'C', 'G']
