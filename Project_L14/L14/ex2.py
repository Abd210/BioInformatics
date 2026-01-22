import re
import math
import random
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt


#Ai simulated text
EMINESCU_TEXT = """
luna cade peste codru si peste lac lin tacut; dorul trece prin frunze,
susura stele rare si vremea curge domol in noapte. iubirea e drum lung,
si visul se intoarce mereu la acelasi gand, ca un ecou pe apa.
"""

STANESCU_TEXT = """
cuvintele se indoaie ca lumina, iar sensul fuge in spatele unei idei;
un mar musca din timp, o metafora se joaca in aerul gol. eu spun: clipa,
si clipa raspunde cu o forma noua. realitatea e o intrebare in miscare.
"""


MIHAI_TEXT = None


WORD_RE = re.compile(r"[a-zA-Zăâîșşțţ]+", re.UNICODE)

def tokenize(text: str):
    return WORD_RE.findall(text.lower())


def bigram_counts(tokens):
    counts = defaultdict(Counter)  # counts[w1][w2]
    totals = Counter()             # totals[w1]
    for w1, w2 in zip(tokens, tokens[1:]):
        counts[w1][w2] += 1
        totals[w1] += 1
    return counts, totals




def p_next(counts, totals, vocab_size, w1, w2, alpha=0.5):
    denom = totals.get(w1, 0) + alpha * vocab_size
    if denom <= 0:
        return 1.0 / vocab_size
    num = counts.get(w1, {}).get(w2, 0) + alpha
    return num / denom




def llr_bigram(w1, w2, counts_plus, totals_plus, counts_minus, totals_minus, vocab_size, alpha=0.5):
    p_plus = p_next(counts_plus, totals_plus, vocab_size, w1, w2, alpha=alpha)
    p_minus = p_next(counts_minus, totals_minus, vocab_size, w1, w2, alpha=alpha)
    return math.log(p_plus / p_minus, 2)




def transition_scores(tokens, counts_plus, totals_plus, counts_minus, totals_minus, vocab_size, alpha=0.5):
    scores = []
    labels = []
    for w1, w2 in zip(tokens, tokens[1:]):
        s = llr_bigram(w1, w2, counts_plus, totals_plus, counts_minus, totals_minus, vocab_size, alpha=alpha)
        scores.append(s)
        labels.append(f"{w1} {w2}")  # you can change to f"{w1}→{w2}" if you prefer
    return scores, labels




def label_scores(scores, eps=0.25):
    labs = []
    for s in scores:
        if s > eps:
            labs.append("+")
        elif s < -eps:
            labs.append("-")
        else:
            labs.append("0")
    return labs

def contiguous_segments(score_labels):
    if not score_labels:
        return []
    segs = []
    start = 0
    cur = score_labels[0]
    for i in range(1, len(score_labels)):
        if score_labels[i] != cur:
            segs.append((start, i, cur))
            start = i
            cur = score_labels[i]
    segs.append((start, len(score_labels), cur))
    return segs



def sample_next_word(counts, w1):
    if w1 not in counts or not counts[w1]:
        return None
    items = list(counts[w1].items())
    words = [w for w, c in items]
    weights = [c for w, c in items]
    return random.choices(words, weights=weights, k=1)[0]

def generate_mihai_blend(counts_plus, counts_minus, vocab, n_words=90, seed=11):
    random.seed(seed)
    vocab = list(vocab)
    w = random.choice(vocab)
    out = [w]

    for i in range(n_words - 1):
        # create "zones":
        # first 30% -> mostly plus
        # middle 30% -> mostly minus
        # last 40% -> mostly plus again
        frac = i / max(1, n_words - 1)
        if frac < 0.30:
            p_plus = 0.75
        elif frac < 0.60:
            p_plus = 0.25
        else:
            p_plus = 0.75

        use_plus = (random.random() < p_plus)
        nxt = sample_next_word(counts_plus if use_plus else counts_minus, w)
        if nxt is None:
            w = random.choice(vocab)
        else:
            w = nxt
        out.append(w)

    return " ".join(out) + "."



def plot_like_screenshot(scores, labels, eps=0.25, max_xticks=45):
    n = len(scores)
    x = list(range(n))

    step = 1
    if n > max_xticks:
        step = math.ceil(n / max_xticks)

    y = scores

    plt.figure(figsize=(14, 6))
    ax = plt.gca()

    ax.plot(x, y, marker="o", linewidth=2, color="#7b2cbf")

    ax.axhline(0.0, linestyle="--", linewidth=1, color="black", alpha=0.6)

    ax.axhline(+eps, linestyle=":", linewidth=1, color="gray", alpha=0.6)
    ax.axhline(-eps, linestyle=":", linewidth=1, color="gray", alpha=0.6)

    y_pos = [v if v > eps else 0 for v in y]
    y_neg = [v if v < -eps else 0 for v in y]

    pos_color = "#4c78a8"  # blue
    neg_color = "#e45756"  # red

    ax.fill_between(x, 0, y_pos, where=[v > eps for v in y], alpha=0.35, label="Eminescu Zone (+)", color=pos_color)
    ax.fill_between(x, 0, y_neg, where=[v < -eps for v in y], alpha=0.35, label="Stănescu Zone (-)", color=neg_color)

    ax.set_title("Plagiarism Detection: Eminescu (+) vs. Stănescu (-)", fontsize=14)
    ax.set_ylabel("Log-Likelihood Score (LLR)")
    ax.set_xlabel("Text Transitions (Scanning the Accused Text)")

    ax.set_xticks(x[::step])
    ax.set_xticklabels([labels[i] for i in x[::step]], rotation=45, ha="right")

    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.show()



tok_plus = tokenize(EMINESCU_TEXT)
tok_minus = tokenize(STANESCU_TEXT)

counts_plus, totals_plus = bigram_counts(tok_plus)
counts_minus, totals_minus = bigram_counts(tok_minus)

vocab = set(tok_plus) | set(tok_minus)
V = len(vocab)

if MIHAI_TEXT is None:
    MIHAI_TEXT = generate_mihai_blend(counts_plus, counts_minus, vocab, n_words=85, seed=11)

tok_m = tokenize(MIHAI_TEXT)

ALPHA = 0.5
EPS = 0.25  # "neither" threshold (tune this)

scores, xlabels = transition_scores(tok_m, counts_plus, totals_plus, counts_minus, totals_minus, V, alpha=ALPHA)

total_score = sum(scores)
print("Mihai text preview:\n", MIHAI_TEXT[:300], "...\n")
print("TOTAL SCORE:", total_score)
print("GLOBAL:", "EMINESCU-LIKE (+)" if total_score > 0 else "STANESCU-LIKE (-)")

plot_like_screenshot(scores, xlabels, eps=EPS, max_xticks=40)

labs = label_scores(scores, eps=EPS)
segs = contiguous_segments(labs)

rows = []
for a, b, lab in segs:
    if lab == "+":
        name = "EMINESCU (+)"
    elif lab == "-":
        name = "STANESCU (-)"
    else:
        name = "NEITHER (0)"
    seg_text = " ".join(tok_m[a:b+1])
    rows.append((a, b, name, seg_text[:160] + ("..." if len(seg_text) > 160 else "")))

df_segments = pd.DataFrame(rows, columns=["start_transition", "end_transition", "label", "preview"])
print("\nSEGMENTS (contiguous zones):")
print(df_segments.to_string(index=False))
