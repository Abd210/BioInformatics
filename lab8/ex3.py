import os
import bisect
import matplotlib.pyplot as plt

GENOME_FILES = [
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\lab8\acute.fasta",
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\lab8\Boebeisivirus.fasta",
    r"C:\Users\Abd\Desktop\sdm\BioInformatics\lab8\corona.fasta",
]

MIN_IR_LEN = 4
MAX_IR_LEN = 6
MAX_PRINT_PER_GENOME = 50

BASE_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def reverse_complement(seq: str) -> str:
    return "".join(BASE_COMP.get(b, "N") for b in reversed(seq))


def read_fasta(path: str):
    header = None
    seq_chunks = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                continue
            seq_chunks.append(line)
    seq = "".join(seq_chunks).upper()
    if header is None:
        header = os.path.basename(path)
    return header, seq


def find_inverted_repeats(
    seq: str,
    min_len: int = MIN_IR_LEN,
    max_len: int = MAX_IR_LEN,
    max_distance: int = None,
):
    n = len(seq)
    if max_distance is None:
        max_distance = n  # allow spacer to go all the way to the end

    results = []
    for k in range(min_len, max_len + 1):
        positions_by_kmer = {}
        for i in range(0, n - k + 1):
            kmer = seq[i : i + k]
            if any(b not in BASE_COMP for b in kmer):
                continue
            positions_by_kmer.setdefault(kmer, []).append(i)

        seen_pairs = set()

        for motif, pos_list in positions_by_kmer.items():
            rc = reverse_complement(motif)
            if rc not in positions_by_kmer:
                continue
            rc_positions = sorted(positions_by_kmer[rc])

            for p in pos_list:
                idx = bisect.bisect_right(rc_positions, p)
                while idx < len(rc_positions):
                    q = rc_positions[idx]
                    spacer_len = q - (p + k)
                    if spacer_len < 0:
                        idx += 1
                        continue
                    if spacer_len > max_distance:
                        break
                    pair_key = (k, p, q)
                    if pair_key not in seen_pairs:
                        seen_pairs.add(pair_key)
                        results.append(
                            {
                                "k": k,
                                "left_start": p,
                                "left_end": p + k - 1,
                                "right_start": q,
                                "right_end": q + k - 1,
                                "spacer_len": spacer_len,
                                "motif": motif,
                                "motif_rc": rc,
                            }
                        )
                    idx += 1
    results.sort(key=lambda r: (r["left_start"], r["right_start"], -r["k"]))
    return results


def plot_ir_pairs(header: str, genome_len: int, ir_pairs):
    if not ir_pairs:
        return
    x = [ir["left_start"] for ir in ir_pairs]
    y = [ir["spacer_len"] for ir in ir_pairs]
    plt.figure(figsize=(10, 4))
    plt.scatter(x, y, s=5)
    plt.xlabel("Left IR start position (bp)")
    plt.ylabel("Spacer length between IRs (bp)")
    plt.title(f"Inverted repeat pairs in {header} (len={genome_len} bp)")
    plt.tight_layout()
    plt.show()


def main():
    for path in GENOME_FILES:
        if not os.path.isfile(path):
            print(f"\n[!] File not found: {path}")
            continue

        print("\n" + "=" * 70)
        print(f"Processing genome file: {path}")
        header, seq = read_fasta(path)
        print(f"Genome header: {header}")
        print(f"Genome length: {len(seq)} bp")

        # max_distance is now the full genome length
        ir_pairs = find_inverted_repeats(seq, max_distance=len(seq))
        print(
            f"Total inverted repeat pairs found (length {MIN_IR_LEN}-{MAX_IR_LEN}): "
            f"{len(ir_pairs)}"
        )

        to_show = ir_pairs[:MAX_PRINT_PER_GENOME]
        if not to_show:
            print("No inverted repeats found with current settings.")
            continue

        print(f"\nShowing first {len(to_show)} possible transposons:")
        print("(Positions are 1-based in this output.)")
        print(
            "index\tIR_len\tleft_start\tleft_end\t"
            "right_start\tright_end\tspacer_len\tmotif\tmotif_rc"
        )

        for idx, ir in enumerate(to_show, start=1):
            print(
                f"{idx}\t"
                f"{ir['k']}\t"
                f"{ir['left_start'] + 1}\t"
                f"{ir['left_end'] + 1}\t"
                f"{ir['right_start'] + 1}\t"
                f"{ir['right_end'] + 1}\t"
                f"{ir['spacer_len']}\t"
                f"{ir['motif']}\t"
                f"{ir['motif_rc']}"
            )

        if len(ir_pairs) > MAX_PRINT_PER_GENOME:
            print(
                f"\n... {len(ir_pairs) - MAX_PRINT_PER_GENOME} more pairs not printed; "
                f"increase MAX_PRINT_PER_GENOME if needed."
            )

        plot_ir_pairs(header, len(seq), ir_pairs)


if __name__ == "__main__":
    main()
