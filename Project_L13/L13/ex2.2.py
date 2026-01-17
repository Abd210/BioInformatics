import random
import re
import json
import string
from collections import defaultdict

def generate_random_english_text(target_chars=300, seed=42):
    """
    Generates ~target_chars of random English-ish text that includes spaces + punctuation.
    """
    random.seed(seed)

    vocab = [
        "the","quick","brown","fox","jumps","over","lazy","dog","today","maybe","because","however",
        "science","music","coffee","dream","cloud","river","city","night","light","small","bright",
        "people","think","feel","walk","run","talk","laugh","write","read","learn","build","create",
        "simple","tiny","funny","strange","wonder","future","memory","pattern","system","signal"
    ]
    punctuation = [".", ",", "!", "?", ";", ":"]

    parts = []
    while len(" ".join(parts)) < target_chars:
        w = random.choice(vocab)

        # Randomly attach punctuation to the word sometimes
        if random.random() < 0.25:
            w += random.choice(punctuation)

        # Randomly capitalize sometimes
        if random.random() < 0.10:
            w = w.capitalize()

        parts.append(w)

    text = " ".join(parts)

    # Trim nicely to <= target_chars (cut at last space so we don't slice a word)
    if len(text) > target_chars:
        cut = text.rfind(" ", 0, target_chars)
        text = text[:cut] if cut != -1 else text[:target_chars]

    # Ensure it ends with punctuation for realism
    if text and text[-1] not in ".!?":
        text += random.choice([".", "!", "?"])

    return text

def tokenize_words(text):
    """
    Extract words (ignoring punctuation for word transitions).
    Keeps apostrophes inside words if any (e.g., don't -> don't).
    """
    return [w.lower() for w in re.findall(r"[A-Za-z']+", text)]

def make_word_symbols(words):
    """
    Map each unique word to an ASCII symbol (string).
    Uses A-Z, a-z, 0-9 first, then falls back to W<number> (still ASCII).
    """
    unique = list(dict.fromkeys(words))  # preserve order of first appearance
    symbol_pool = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)  # 62 symbols

    word_to_sym = {}
    sym_to_word = {}

    for i, w in enumerate(unique):
        sym = symbol_pool[i] if i < len(symbol_pool) else f"W{i}"  # still ASCII text
        word_to_sym[w] = sym
        sym_to_word[sym] = w

    return word_to_sym, sym_to_word

def transition_probabilities_between_words(words, word_to_sym):
    """
    Build transition probabilities between word symbols.
    Returns dict: from_sym -> {to_sym: prob}
    """
    counts = defaultdict(lambda: defaultdict(int))

    for i in range(len(words) - 1):
        a, b = words[i], words[i + 1]
        sa, sb = word_to_sym[a], word_to_sym[b]
        counts[sa][sb] += 1

    probs = {}
    for sa, row in counts.items():
        total = sum(row.values())
        probs[sa] = {sb: c / total for sb, c in row.items()} if total else {}

    return probs

def save_to_json(payload, filename="word_transition_matrix.json"):
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=True)

if __name__ == "__main__":
    text = generate_random_english_text(target_chars=300, seed=7)
    words = tokenize_words(text)

    word_to_sym, sym_to_word = make_word_symbols(words)
    transition_matrix = transition_probabilities_between_words(words, word_to_sym)

    output = {
        "random_text": text,
        "tokenized_words": words,
        "word_to_symbol": word_to_sym,
        "symbol_to_word": sym_to_word,
        "transition_matrix": transition_matrix
    }

    print("Random text:\n", text)
    print("\n#words:", len(words), "| unique:", len(word_to_sym))
    print("\nSaved JSON -> word_transition_matrix.json")

    save_to_json(output, "word_transition_matrix.json")
