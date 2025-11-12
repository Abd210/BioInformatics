import re, argparse
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

def load_multifasta(path):
    out=[]
    for rec in SeqIO.parse(path,"fasta"):
        out.append((rec.id, re.sub(r"[^ACGTN]","",str(rec.seq).upper())))
    if not out: raise RuntimeError("No sequences in FASTA.")
    return out

def find_tandem_repeat_blocks(seq,k_min=2,k_max=6,min_copies=2):
    valid=set("ACGT"); s=seq.upper(); n=len(s)
    blocks=Counter(); i=0
    while i<n:
        advanced=False
        for k in range(k_min,k_max+1):
            if i+k*min_copies>n: continue
            m=s[i:i+k]
            if set(m)-valid: continue
            copies=1; j=i+k
            while j+k<=n and s[j:j+k]==m: copies+=1; j+=k
            if copies>=min_copies:
                blocks[m]+=1; i=j; advanced=True; break
        if not advanced: i+=1
    return blocks

def freqs(blocks, L, min_copies):
    return {m: c/max(1, L-len(m)*min_copies+1) for m,c in blocks.items()}

def main():
    p=argparse.ArgumentParser(description="One-row plots of tandem-repeat frequencies (motif 2..6 bp).")
    p.add_argument("--fasta", default=r"C:\Users\Abd\Desktop\sdm\BioInformatics\sequence.fasta")
    p.add_argument("--kmin", type=int, default=2)         
    p.add_argument("--kmax", type=int, default=6)
    p.add_argument("--min-copies", type=int, default=2)
    p.add_argument("--influenza-n", type=int, default=10)
    p.add_argument("--topn", type=int, default=20)
    a=p.parse_args()

    recs=load_multifasta(a.fasta)
    idx=next((i for i,(_,s) in enumerate(recs) if 100<=len(s)<=3000), None)
    if idx is None: raise RuntimeError("No 100–3000 nt sequence found.")
    seq_list=[recs[idx]]+ [r for i,r in enumerate(recs) if i!=idx][:a.influenza_n]

    panels=[]
    for hdr,seq in seq_list:
        b=find_tandem_repeat_blocks(seq,a.kmin,a.kmax,a.min_copies)
        f=freqs(b,len(seq),a.min_copies)
        items=sorted(f.items(), key=lambda kv:(-kv[1],kv[0]))[:a.topn]
        labels=[m for m,_ in items]; values=[v for _,v in items]
        panels.append((hdr,len(seq),labels,values))

    n=len(panels)
    fig,axes=plt.subplots(1,n, figsize=(max(8,4*n),5), sharey=True)
    if n==1: axes=[axes]
    for ax,(hdr,L,labels,vals) in zip(axes,panels):
        ax.bar(range(len(vals)), vals)
        ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_title(f"{hdr} (len={L})", fontsize=10)
        ax.set_xlabel("Motif 2–6 bp", fontsize=9)
    axes[0].set_ylabel("Frequency (blocks / (len - k*min + 1))", fontsize=9)
    fig.suptitle("Tandem-repeat frequencies (motif 2–6 bp, one row)", fontsize=14)
    fig.tight_layout(rect=[0,0,1,0.93]); plt.show()

if __name__=="__main__": main()
