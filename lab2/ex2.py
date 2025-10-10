
S = "ATTGTCCCAATCTGTTG"

# Dinucleotides
print("Dinucleotides:")
dinucl_counts = {}
for i in range(len(S)-1):  
    dinucl = S[i:i+2]
    if dinucl in dinucl_counts:
        dinucl_counts[dinucl] += 1
    else:
        dinucl_counts[dinucl] = 1

total_dinucl = sum(dinucl_counts.values())
for dinucl, count in dinucl_counts.items():
    percent = (count / total_dinucl) * 100
    print(dinucl, count, f"{percent:.2f}%")

# Trinucleot
print("\nTrinucleotides:")
trinucl_counts = {}
for i in range(len(S)-2):  
    trinucl = S[i:i+3]
    if trinucl in trinucl_counts:
        trinucl_counts[trinucl] += 1
    else:
        trinucl_counts[trinucl] = 1

total_trinuc = sum(trinucl_counts.values())
for trinuc, count in trinucl_counts.items():
    percent = (count / total_trinuc) * 100
    print(trinuc, count, f"{percent:.2f}%")



