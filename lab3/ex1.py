# DNA melting temperature calculator
# Formula 1: Tm = 4(G + C) + 2(A + T)
# Formula 2: Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/length

import math

# Input DNA sequence
dna = input("Enter DNA sequence: ").upper()

# Count bases
A = dna.count('A')
T = dna.count('T')
G = dna.count('G')
C = dna.count('C')

Tm1 = 4*(G + C) + 2*(A + T)

Na = 0.001
length = len(dna)
GC_percent = ((G + C) / length) * 100
Tm2 = -(81.5 + 16.6*math.log10(Na) + 0.41*GC_percent - 600/length)

# Output
print(f"Simple formula Tm = {Tm1:.2f} °C")
print(f"Salt-adjusted formula Tm = {Tm2:.2f} °C")
