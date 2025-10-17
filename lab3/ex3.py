# establish a treshold that acts as a cutoff value for each signal. plot a chart that shows the cutoff values as 
# horizontal bars
# either consider a common treshold or one for each signal

import matplotlib.pyplot as plt
import math

sliding_window_size = 9
step_size = 1   

def read_fasta(path):
    with open(path, "r") as f:
        lines = f.readlines()
    seq = ""
    for line in lines:
        if line.startswith(">"):
            continue
        seq += line.strip()
    return seq

def calculate_melting_temps(S):
    P1 = []
    P2 = []

    n = len(S)
    for i in range(0, n - sliding_window_size + 1, step_size):
        window = S[i:i+ sliding_window_size]
        P1.append(calculate_melting_temp_1(window))
        P2.append(calculate_melting_temp_2(window))
    
    return P1, P2

def calculate_melting_temp_1(S):
    a_number, c_number, g_number, t_number = count_molecules(S)
    total_bases = len(S)
    melting_temp = -(81.5 + 16.6 * math.log10(0.001) + 0.41 * (g_number+c_number )/total_bases*100 - (600/total_bases))
    return melting_temp

def calculate_melting_temp_2(S):
    a_number, c_number, g_number, t_number = count_molecules(S)
    melting_temp = 4 * (g_number + c_number) + 2 * (a_number + t_number)
    return melting_temp

def count_molecules(S):
    a_number = S.count("A")
    c_number = S.count("C")
    g_number = S.count("G")
    t_number = S.count("T")
    return a_number, c_number, g_number, t_number

def plot_threshold_bars(P1, P2, threshold_P1, threshold_P2, mode="above"):
    """
    Plot horizontal bars showing regions where signals are above/below thresholds
    mode: "above" or "below" - show bars where condition is met
    """
    # Create binary indicators (1 where condition is met, 0 where not)
    P1_bars = []
    P2_bars = []
    
    for i in range(len(P1)):
        if mode == "above":
            P1_bars.append(1 if P1[i] >= threshold_P1 else 0)
            P2_bars.append(1 if P2[i] >= threshold_P2 else 0)
        elif mode == "below":
            P1_bars.append(1 if P1[i] < threshold_P1 else 0)
            P2_bars.append(1 if P2[i] < threshold_P2 else 0)
    
    # Plot bars
    plt.figure(figsize=(12, 4))
    
    # Plot P1 bars on top row
    plt.subplot(2, 1, 1)
    for i in range(len(P1_bars)):
        if P1_bars[i] == 1:
            plt.barh(0, 1, left=i, height=0.03, color='blue', alpha=0.8)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0, len(P1))
    plt.title(f'P1 - Regions {mode.capitalize()} Threshold ({threshold_P1:.1f}°C)')
    plt.ylabel('P1')
    plt.yticks([])  # Remove y-axis ticks
    
    # Plot P2 bars on bottom row  
    plt.subplot(2, 1, 2)
    for i in range(len(P2_bars)):
        if P2_bars[i] == 1:
            plt.barh(0, 1, left=i, height=0.03, color='orange', alpha=0.8)
    plt.ylim(-0.2, 0.2)
    plt.xlim(0, len(P2))
    plt.title(f'P2 - Regions {mode.capitalize()} Threshold ({threshold_P2:.1f}°C)')
    plt.ylabel('P2')
    plt.yticks([])  # Remove y-axis ticks
    plt.xlabel('Window Number')
    
    plt.tight_layout()
    plt.show()
    
    # Print statistics
    P1_count = sum(P1_bars)
    P2_count = sum(P2_bars)
    print(f"Mode: {mode}")
    print(f"P1 windows {mode} threshold: {P1_count} out of {len(P1)} ({P1_count/len(P1)*100:.1f}%)")
    print(f"P2 windows {mode} threshold: {P2_count} out of {len(P2)} ({P2_count/len(P2)*100:.1f}%)")

# Calculate P1 and P2 vectors
dna_sequence = read_fasta("dna.fasta")
P1, P2 = calculate_melting_temps(dna_sequence)

# Calculate thresholds
threshold_P1 = sum(P1) / len(P1)  # Average of P1
threshold_P2 = sum(P2) / len(P2)  # Average of P2

print(f"P1 Threshold (average): {threshold_P1:.2f}°C")
print(f"P2 Threshold (average): {threshold_P2:.2f}°C")

# Example usage - you can change "above" to "below"
plot_threshold_bars(P1, P2, threshold_P1, threshold_P2, mode="above")

# Uncomment this line to see regions below thresholds instead:
# plot_threshold_bars(P1, P2, threshold_P1, threshold_P2, mode="below")