import json, random, textwrap
from pathlib import Path

DNA_JSON = Path(r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L13\L13\transition_matrix.json")          # <-- change if your DNA file name differs
WORD_JSON = Path("C:\LaylaYear4\BioInformatics\BioInformatics\Project_L13\L13\word_transition_matrix.json")    # from the earlier word task

def wrap(s, width=90):
    return "\n".join(textwrap.wrap(s, width=width))

def sample_from_probs(prob_dict):
    """Sample a key from {key: prob}. Works even if probs aren't perfectly normalized."""
    items = list(prob_dict.items())
    if not items:
        return None
    keys = [k for k, _ in items]
    probs = [float(p) for _, p in items]
    total = sum(probs)
    if total <= 0:
        return random.choice(keys)

    r = random.random() * total
    c = 0.0
    for k, p in zip(keys, probs):
        c += p
        if r <= c:
            return k
    return keys[-1]

def synthesize_markov(matrix, states, length, start_state=None):
    if not states:
        raise ValueError("No states available to synthesize from.")
    if start_state is None:
        start_state = random.choice(states)

    seq = [start_state]
    for _ in range(length - 1):
        cur = seq[-1]
        row = matrix.get(cur, {})
        nxt = sample_from_probs(row) if row else None
        if nxt is None:
            nxt = random.choice(states)
        seq.append(nxt)
    return seq

def top_transitions(matrix, state, k=4):
    row = matrix.get(state, {})
    items = sorted(row.items(), key=lambda x: float(x[1]), reverse=True)
    return [(n, float(p)) for n, p in items[:k] if float(p) > 0]

# -------------------- DNA: BEFORE & AFTER --------------------
# Supports BOTH formats:
#  (1) plain dict matrix: {"A":{"C":0.2,...}, ...}
#  (2) wrapped dict: {"transition_matrix": {...}, "sequence": "..."} etc.
dna_obj = json.loads(DNA_JSON.read_text(encoding="utf-8"))

if "transition_matrix" in dna_obj:
    dna_matrix = dna_obj["transition_matrix"]
    dna_seq_before = dna_obj.get("sequence") or dna_obj.get("dna") or ""
else:
    dna_matrix = dna_obj
    dna_seq_before = ""

dna_states = list(dna_matrix.keys())

print("\n" + "="*90)
print("DNA MARKOV MODEL (letters)")
print("="*90)

if dna_seq_before:
    print("\nBEFORE (original DNA used to build matrix)")
    print(f"Length: {len(dna_seq_before)}")
    print(wrap(dna_seq_before, 90))
else:
    print("\nBEFORE")
    print("(No DNA sequence stored in this JSON — only the transition matrix.)")

print("\nWhat the model learned (top transitions per base)")
for s in dna_states:
    tops = top_transitions(dna_matrix, s, k=4)
    if not tops:
        print(f"  {s} -> (no outgoing transitions observed)")
    else:
        pretty = ", ".join([f"{n}:{p:.2f}" for n, p in tops])
        print(f"  {s} -> {pretty}")

print("\nAFTER (synthesized DNA samples)")
for i in range(1, 4):
    new_dna = "".join(synthesize_markov(dna_matrix, dna_states, length=80))
    print(f"\nSample #{i} (80 letters):")
    print(wrap(new_dna, 90))

# -------------------- TEXT: BEFORE & AFTER --------------------
words_obj = json.loads(WORD_JSON.read_text(encoding="utf-8"))

# From OUR earlier JSON:
text_before = words_obj.get("random_text", "")
symbol_to_word = words_obj["symbol_to_word"]
symbols = list(symbol_to_word.keys())
word_matrix = words_obj["transition_matrix"]

print("\n" + "="*90)
print("TEXT MARKOV MODEL (words)")
print("="*90)

print("\nBEFORE (original text used to build matrix)")
if text_before:
    print(f"Chars: {len(text_before)}")
    print(wrap(text_before, 90))
else:
    print("(No original text stored in this JSON.)")

# Choose a start symbol (try starting from a word that was capitalized in the original text)
# Our stored mapping is lowercased words, so this is just a simple random start:
start_sym = random.choice(symbols) if symbols else None

print("\nWhat the model learned (a few example transitions)")
# Pick a few random symbols to preview transitions:
for sym in random.sample(symbols, k=min(7, len(symbols))):
    tops = top_transitions(word_matrix, sym, k=6)
    if not tops:
        continue
    pretty = ", ".join([f"{symbol_to_word[n]}:{p:.2f}" for n, p in tops])
    print(f"  '{symbol_to_word[sym]}' -> {pretty}")

print("\nAFTER (synthesized text samples)")
for i in range(1, 4):
    sym_seq = synthesize_markov(word_matrix, symbols, length=60, start_state=start_sym)
    tok_seq = [symbol_to_word[s] for s in sym_seq]
    new_text = " ".join(tok_seq)
    print(f"\nSample #{i} (60 words):")
    print(wrap(new_text, 90))

print("\n" + "="*90)
print("Done.")
print("="*90)
