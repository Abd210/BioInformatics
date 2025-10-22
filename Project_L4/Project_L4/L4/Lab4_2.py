import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
from Bio.Seq import Seq
import numpy as np  # For side-by-side bar positions

# Function to read FASTA file and extract sequence
def read_fasta(file_path):
    with open(file_path, "r") as file:
        seq_record = SeqIO.read(file, "fasta")
    return str(seq_record.seq)

# Convert DNA (T) to RNA (U)
def dna_to_rna(dna_sequence):
    return dna_sequence.replace("T", "U")

# Calculate codon frequencies
def calculate_codon_frequencies(sequence):
    codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
    return Counter(codons)

# Plot Top N Codons
def plot_top_codons(codon_counts, title, top_n=10):
    codons, counts = zip(*codon_counts.most_common(top_n))
    plt.figure(figsize=(10, 6))
    plt.bar(codons, counts, color='maroon')
    plt.title(title)
    plt.xlabel('Codon')
    plt.ylabel('Frequency')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

# Read genome sequences
covid_sequence = read_fasta("Project_L4/L4/sequenceCOVID.fasta")
flu_sequence = read_fasta("Project_L4/L4/sequenceINFLUENZA.fasta")

# Convert DNA to RNA if needed
covid_sequence_rna = dna_to_rna(covid_sequence)
flu_sequence_rna = dna_to_rna(flu_sequence)

# Calculate codon frequencies for RNA sequences
covid_codon_counts = calculate_codon_frequencies(covid_sequence_rna)
flu_codon_counts = calculate_codon_frequencies(flu_sequence_rna)

# Plot top 10 codons for each genome
plot_top_codons(covid_codon_counts, "Top 10 Codons in SARS-CoV-2 Genome")
plot_top_codons(flu_codon_counts, "Top 10 Codons in Influenza A Genome")

# Compare common codons between both genomes
common_codons = set(covid_codon_counts.keys()) & set(flu_codon_counts.keys())
print("Common Codons:", common_codons)

# Translate to amino acids and display top 3 for each genome
def translate_to_amino_acids(codon_counts):
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        if codon != 'Stop':
            aa = str(Seq(codon).translate())
            aa_counts[aa] += count
    return aa_counts

covid_aa_counts = translate_to_amino_acids(covid_codon_counts)
flu_aa_counts = translate_to_amino_acids(flu_codon_counts)

print("Top 3 Amino Acids in SARS-CoV-2:", covid_aa_counts.most_common(3))
print("Top 3 Amino Acids in Influenza A:", flu_aa_counts.most_common(3))

# 🔹 NEW: Compare and plot most frequent shared codons side by side
def compare_common_codons_side_by_side(covid_counts, flu_counts, top_n=10):
    # Find common codons
    common = set(covid_counts.keys()) & set(flu_counts.keys())

    # Get top N common codons by combined frequency
    combined_counts = Counter()
    for codon in common:
        combined_counts[codon] = covid_counts[codon] + flu_counts[codon]
    top_codons = [codon for codon, _ in combined_counts.most_common(top_n)]

    # Get counts for COVID and Influenza
    covid_values = [covid_counts[codon] for codon in top_codons]
    flu_values = [flu_counts[codon] for codon in top_codons]

    # Bar positions
    x = np.arange(len(top_codons))
    width = 0.35

    # Plot side-by-side bars
    plt.figure(figsize=(12, 6))
    plt.bar(x - width/2, covid_values, width, label='COVID', color='darkmagenta')
    plt.bar(x + width/2, flu_values, width, label='Influenza', color='green')
    plt.title(f"Top {top_n} Shared Codons: COVID vs Influenza")
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.xticks(x, top_codons, rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.show()

# Plot the third comparison chart
compare_common_codons_side_by_side(covid_codon_counts, flu_codon_counts, top_n=10)
