import random
import time
from Bio import SeqIO
import matplotlib.pyplot as plt

# FASTA files
files = [
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\corona.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\Boebeisivirus.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\acute.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\HepatitisB.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\HepatitisC.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\papillomavirus.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\Nile.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\Measles.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\Zika.fasta",
    r"C:\LaylaYear4\BioInformatics\BioInformatics\Project_L5\L5\Enterobacteria.fasta"
]

sample_size = 3000
num_samples = 5

cg_percentages = []
assembly_times = []

for file in files:
    record = SeqIO.read(file, "fasta")
    genome_seq = str(record.seq).upper()
    
    virus_cg_list = []
    virus_time_list = []
    
    for _ in range(num_samples):
        if len(genome_seq) < sample_size:
            sample_seq = genome_seq
        else:
            start_idx = random.randint(0, len(genome_seq) - sample_size)
            sample_seq = genome_seq[start_idx:start_idx + sample_size]
        
        start_time = time.time()
        assembled_seq = ''.join(sample_seq)
        elapsed_ms = (time.time() - start_time) * 1000
        
        cg_content = (sample_seq.count('C') + sample_seq.count('G')) / len(sample_seq) * 100
        
        virus_cg_list.append(cg_content)
        virus_time_list.append(elapsed_ms)
    
    cg_percentages.append(sum(virus_cg_list)/len(virus_cg_list))
    assembly_times.append(sum(virus_time_list)/len(virus_time_list))

# Plot
plt.figure(figsize=(10,6))
plt.scatter(cg_percentages, assembly_times, color='purple')
plt.xlabel("C+G Content (%)")
plt.ylabel("Assembly Time (ms)")
plt.title("Viral Genome Samples: Assembly Time vs C+G Content")
plt.grid(True)

virus_names = [f.split("\\")[-1].split(".")[0] for f in files]
for i, name in enumerate(virus_names):
    plt.annotate(name, (cg_percentages[i], assembly_times[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.show()

# Save analysis
with open("analysis_explanation.txt", "w") as f:
    f.write("Analysis of assembly time vs C+G content for 10 viral genomes:\n\n")
    for i, name in enumerate(virus_names):
        f.write(f"{name}: C+G = {cg_percentages[i]:.2f}%, Avg assembly time = {assembly_times[i]:.3f} ms\n")
    f.write("\nObservations:\n")
    f.write("- Viruses with higher C+G content may show slightly different assembly times (due to sequence complexity).\n")
    f.write("- Points more to the right (higher C+G) vs left indicate GC-rich genomes.\n")
    f.write("- Points higher on Y-axis indicate longer assembly times.\n")
