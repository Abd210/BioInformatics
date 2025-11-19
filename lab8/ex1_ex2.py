import random
import textwrap


NUCLEOTIDES = "ACGT"

IR1 = "ATGCA"
IR2 = "TACGT"

MIN_INTERNAL = 8
MAX_INTERNAL = 20


def random_dna(length: int) -> str:
    return "".join(random.choice(NUCLEOTIDES) for _ in range(length))


def build_sequence():

    while True:
        num_tes = random.randint(3, 4) 
        seq_parts = []
        tes_truth = []
        current_len = 0

        for _ in range(num_tes):
            bg_len = random.randint(10, 40)
            bg = random_dna(bg_len)
            seq_parts.append(bg)
            current_len += bg_len
            internal_len = random.randint(MIN_INTERNAL, MAX_INTERNAL)
            internal = random_dna(internal_len)
            te_type = random.choice(["IR1-IR2", "IR2-IR1"])

            if te_type == "IR1-IR2":
                te_seq = IR1 + internal + IR2
            else: 
                te_seq = IR2 + internal + IR1

            start = current_len
            end = current_len + len(te_seq) - 1  

            seq_parts.append(te_seq)
            current_len = end + 1

            tes_truth.append({
                "type": te_type,
                "start0": start,
                "end0": end
            })

        bg_len = random.randint(10, 40)
        bg = random_dna(bg_len)
        seq_parts.append(bg)
        current_len += bg_len

        if 200 <= current_len <= 400:
            seq = "".join(seq_parts)
            return seq, tes_truth



def detect_tes(seq: str):
    """
    Detect positions of transposable elements based on IR1/IR2 patterns.
    A TE is:
      IR1-IR2: IR1 ... IR2  with internal length in [MIN_INTERNAL, MAX_INTERNAL]
      IR2-IR1: IR2 ... IR1  with internal length in [MIN_INTERNAL, MAX_INTERNAL]
    Returns:
      list of dicts like build_sequence() ground truth.
    """
    results = []
    n = len(seq)
    len1 = len(IR1)
    len2 = len(IR2)

    ir1_positions = [i for i in range(n - len1 + 1) if seq[i:i + len1] == IR1]
    ir2_positions = [i for i in range(n - len2 + 1) if seq[i:i + len2] == IR2]

    for i in ir1_positions:
        start_internal = i + len1 + MIN_INTERNAL
        end_internal = min(i + len1 + MAX_INTERNAL, n - len2)
        for j in range(start_internal, end_internal + 1):
            if seq[j:j + len2] == IR2:
                results.append({
                    "type": "IR1-IR2",
                    "start0": i,
                    "end0": j + len2 - 1
                })
                break  

    for i in ir2_positions:
        start_internal = i + len2 + MIN_INTERNAL
        end_internal = min(i + len2 + MAX_INTERNAL, n - len1)
        for j in range(start_internal, end_internal + 1):
            if seq[j:j + len1] == IR1:
                results.append({
                    "type": "IR2-IR1",
                    "start0": i,
                    "end0": j + len1 - 1
                })
                break

    return results


if __name__ == "__main__":
    random.seed() 

    sequence, ground_truth_tes = build_sequence()

    detected_tes = detect_tes(sequence)

    print("DNA length:", len(sequence))
    print("\nDNA sequence (wrapped at 60 bp):")
    print("\n".join(textwrap.wrap(sequence, 60)))

    print("\nGround truth transposable elements (1-based positions):")
    for idx, te in enumerate(ground_truth_tes, start=1):
        print(
            f"  TE{idx} {te['type']}: "
            f"{te['start0'] + 1} – {te['end0'] + 1}"
        )

    print("\nDetected transposable elements (1-based positions):")
    for idx, te in enumerate(detected_tes, start=1):
        print(
            f"  TE{idx} {te['type']}: "
            f"{te['start0'] + 1} – {te['end0'] + 1}"
        )
