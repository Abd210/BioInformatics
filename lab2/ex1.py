S = "ATTGTCCCAATCTGTTG"
bases = ["A", "C", "G", "T"]

#2 letters
dinucl = []
for nr1 in bases:
    for nr2 in bases:
        dinucl.append(nr1 + nr2)

#3 letters
trinucl = []
for nr1 in bases:
    for nr2 in bases:
        for nr3 in bases:
            trinucl.append(nr1 + nr2 + nr3)

# k = 2  -> total windows = len(S)-2+1
# k = 3  -> total windows = len(S)-3+1
# General: total windows = len(S) - k + 1

def brute_force_engine(S, leng, combos):
    total = len(S) - leng + 1
    for c in combos:
        count = S.count(c)
        percent = (count / total) * 100
        print(f"{c}\t{count}\t{percent:.2f}%")

print("=== Dinucleotides (k=2) ===")
brute_force_engine(S, 2, dinucl)

print("\n=== Trinucleotides (k=3) ===")
brute_force_engine(S, 3, trinucl)