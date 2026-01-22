import math
import pandas as pd
import matplotlib.pyplot as plt



S1 = "ATCGATTCGATATCATACACGTAT"          
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"         
S  = "CAGGTTGGAAACGTAA"                  

ALPH = "ACGT"
IDX = {ch: i for i, ch in enumerate(ALPH)}



def transition_counts(seq: str):
    counts = [[0]*4 for _ in range(4)]
    for a, b in zip(seq, seq[1:]):
        counts[IDX[a]][IDX[b]] += 1
    return counts



def transition_probs(counts):
    probs = [[0.0]*4 for _ in range(4)]
    for i in range(4):
        row_sum = sum(counts[i])
        for j in range(4):
            probs[i][j] = counts[i][j] / row_sum if row_sum else 0.0
    return probs



def llr_matrix(p_plus, p_minus, base=2.0):
    ln_base = math.log(base)
    beta = [[0.0]*4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            pp, pm = p_plus[i][j], p_minus[i][j]
            if pp > 0.0 and pm > 0.0:
                beta[i][j] = math.log(pp / pm) / ln_base
            else:
                beta[i][j] = 0.0
    return beta


def score_sequence(seq, beta):
    total = 0.0
    rows = []
    for t, (a, b) in enumerate(zip(seq, seq[1:]), start=1):
        v = beta[IDX[a]][IDX[b]]
        total += v
        rows.append((t, a, b, f"{a}→{b}", v, total))
    df = pd.DataFrame(rows, columns=["step", "from", "to", "pair", "beta", "cumulative"])
    return total, df


def plot_matrix(mat_df: pd.DataFrame, title: str):
    plt.figure(figsize=(5, 4))
    ax = plt.gca()
    im = ax.imshow(mat_df.values)  # default colormap
    ax.set_title(title)
    ax.set_xticks(range(4)); ax.set_yticks(range(4))
    ax.set_xticklabels(list(ALPH)); ax.set_yticklabels(list(ALPH))
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for i in range(4):
        for j in range(4):
            val = mat_df.values[i][j]
            ax.text(j, i, f"{val:.3f}", ha="center", va="center")
    plt.tight_layout()
    plt.show()

def plot_cumulative(df_steps: pd.DataFrame, title="Cumulative log-likelihood score along S"):
    plt.figure(figsize=(8, 4))
    ax = plt.gca()
    ax.plot(df_steps["step"], df_steps["cumulative"], marker="o")
    ax.set_title(title)
    ax.set_xlabel("Transition step")
    ax.set_ylabel("Cumulative score")
    ax.set_xticks(df_steps["step"])
    ax.set_xticklabels(df_steps["pair"], rotation=45, ha="right")
    plt.tight_layout()
    plt.show()

def plot_contributions_bar(df_steps: pd.DataFrame, title="Per-transition β contributions (which transitions matter most)"):
    # Sort by absolute contribution so big positives/negatives stand out
    df_sorted = df_steps.copy()
    df_sorted["abs_beta"] = df_sorted["beta"].abs()
    df_sorted = df_sorted.sort_values("abs_beta", ascending=False)

    plt.figure(figsize=(8, 4))
    ax = plt.gca()
    ax.bar(range(len(df_sorted)), df_sorted["beta"].values)
    ax.set_title(title)
    ax.set_xlabel("Transitions (sorted by |β|)")
    ax.set_ylabel("β value")
    ax.set_xticks(range(len(df_sorted)))
    ax.set_xticklabels(df_sorted["pair"].tolist(), rotation=45, ha="right")
    plt.tight_layout()
    plt.show()


c_plus  = transition_counts(S1)
c_minus = transition_counts(S2)

p_plus  = transition_probs(c_plus)
p_minus = transition_probs(c_minus)

beta = llr_matrix(p_plus, p_minus, base=2.0)

# Make DataFrames
df_counts_plus  = pd.DataFrame(c_plus,  index=list(ALPH), columns=list(ALPH))
df_counts_minus = pd.DataFrame(c_minus, index=list(ALPH), columns=list(ALPH))
df_probs_plus   = pd.DataFrame(p_plus,  index=list(ALPH), columns=list(ALPH))
df_probs_minus  = pd.DataFrame(p_minus, index=list(ALPH), columns=list(ALPH))
df_beta         = pd.DataFrame(beta,    index=list(ALPH), columns=list(ALPH))

# Score new sequence
score, df_steps = score_sequence(S, beta)


print("S1 (CpG+):", S1)
print("S2 (CpG-):", S2)
print("S  (new): ", S)
print("\nFinal log-likelihood score for S =", score)
print("Decision:", "CpG ISLAND (+)" if score > 0 else "NON-ISLAND (-)")

print("\nCounts Tr+ (from S1):\n", df_counts_plus)
print("\nCounts Tr- (from S2):\n", df_counts_minus)
print("\nProbabilities Tr+:\n", df_probs_plus.round(3))
print("\nProbabilities Tr-:\n", df_probs_minus.round(3))
print("\nβ matrix = log2(Tr+/Tr-):\n", df_beta.round(3))

print("\nPer-transition contributions:\n", df_steps[["step", "pair", "beta", "cumulative"]].round(6))


plot_matrix(df_probs_plus,  "Tr+ transition probabilities (from S1)")
plot_matrix(df_probs_minus, "Tr- transition probabilities (from S2)")
plot_matrix(df_beta,        "β log-likelihood matrix: log2(Tr+/Tr-)")

plot_cumulative(df_steps)
plot_contributions_bar(df_steps)
