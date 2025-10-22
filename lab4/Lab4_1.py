
S="UGCAUGUCGCCGAUAUAAUAUCAU"


# Genetic Code Dictionary
genetic_code = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

def dna_to_rna(dna_sequence):
    return dna_sequence.replace("T", "U")


dna_to_rna(S)

def translate_to_protein(sequence):
    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

    protein = []
    started = False
    for codon in codons:
        if codon == "AUG":  # Start translation at 'AUG' (Methionine)
            started = True
            protein.append(genetic_code.get(codon, "Invalid Codon"))
            continue
        if started:
            if codon in ["UAA", "UAG", "UGA"]:
                protein.append(genetic_code.get(codon, "Invalid Codon"))
                 # Stop translation at a stop codon
                break
            protein.append(genetic_code.get(codon, "Invalid Codon"))
    
    return protein

dna_to_rna(S)

S = "UGC"+"AUG"+"UCG"+"CCG"+"AUA"+"UAA"+"UAU"+"CAU"

protein_sequence = translate_to_protein(S)

print("Amino Acid Sequence: ", protein_sequence)





'''
Implement an app that converts the coding region of a gene into an aminoacid sequence. 
use the genetic code tabe from moodle

'''