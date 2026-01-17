'''
1. Use a random DNA sequence of about 50 letters. Use this sequence to compute the transiotion probabilities between letters.
Your output should be  the transition matrix stored as a JSON file

2. Use a random english text of about 300 letters, that implies spaces, punctuation and compute the transition probabilities between words,
Store the transition matrix as a Json file. for ease of implementation, you can represent each new word by using a symbol of your choice(ascii)

'''

import random
import json

def generate_dna_sequence(length=50, seed=None):
    if seed is not None:
        random.seed(seed)
    letters = ["A", "C", "G", "T"]
    return "".join(random.choice(letters) for _ in range(length))

def transition_matrix_from_sequence(seq, letters=("A", "C", "G", "T")):
    # Count transitions: counts[from][to]
    counts = {a: {b: 0 for b in letters} for a in letters}

    for i in range(len(seq) - 1):
        a, b = seq[i], seq[i + 1]
        if a in counts and b in counts[a]:
            counts[a][b] += 1

    # Convert counts to probabilities row-wise
    probs = {}
    for a in letters:
        row_total = sum(counts[a].values())
        if row_total == 0:
            # If this letter never appeared as a "from" state, set all probs to 0
            probs[a] = {b: 0.0 for b in letters}
        else:
            probs[a] = {b: counts[a][b] / row_total for b in letters}

    return probs

def save_matrix_to_json(matrix, filename="transition_matrix.json"):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(matrix, f, indent=2)

if __name__ == "__main__":
    dna = generate_dna_sequence(length=50, seed=42)  
    matrix = transition_matrix_from_sequence(dna)

    print("DNA sequence:", dna)
    print("Transition probabilities:")
    for from_letter, row in matrix.items():
        print(from_letter, row)

    save_matrix_to_json(matrix, "transition_matrix.json")
    print("\nSaved to: transition_matrix.json")
