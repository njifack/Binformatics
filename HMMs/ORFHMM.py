import numpy as np
import matplotlib.pyplot as plt

# Define parameters
states = ['S', 'C', 'T']
x = 50  # Total emissions

# Start probabilities
k_0 = np.array([1, 0, 0])

# Transition matrix
A = np.array([
    [0.05, 0.94, 0.01],
    [0.0, 0.95, 0.05],
    [0.0, 0.0, 1.0]
])

# Simulate sequences
def generate_sequence():
    sequence = []
    current_state = np.random.choice(states, p=k_0)
    for _ in range(x):
        if current_state == 'S':
            sequence.append('ATG')
            current_state = np.random.choice(states, p=A[0])
        elif current_state == 'C':
            # Assuming dummy codons for simplicity
            sequence.append('XXX')
            current_state = np.random.choice(states, p=A[1])
        elif current_state == 'T':
            sequence.append('TAA')
            current_state = np.random.choice(states, p=A[2])
    return sequence

# Generate sequences and measure ORF lengths
orf_lengths = []
for _ in range(250):
    seq = generate_sequence()
    orf_length = 0
    for codon in seq:
        if codon == 'ATG':
            orf_length = 1
        elif codon == 'TAA':
            if orf_length > 0:
                orf_lengths.append(orf_length)
                break
            orf_length = 0
        elif orf_length > 0:
            orf_length += 1


# Write ORF lengths to a file
with open("Evanessa_ORFoutput.txt", "w") as file:
    for length in orf_lengths:
        file.write(f"{length}\n")


# Plot distribution of ORF lengths
plt.hist(orf_lengths, edgecolor='black', bins=20)
plt.title('Distribution of ORF Lengths')
plt.xlabel('ORF Length')
plt.ylabel('Frequency')
plt.show()
