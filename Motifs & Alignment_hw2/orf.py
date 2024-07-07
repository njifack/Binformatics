def scanSeq(sequence):
    start_positions = []
    orf_lengths = []
    orf_sequences = []

    min_orf_length = 60
    codon_length = 3

    for i in range(len(sequence) - min_orf_length + 1):
        if sequence[i:i+codon_length] in ["ATG", "GTG"]:
            orf_start = i
            orf_length = codon_length

            for j in range(i+codon_length, len(sequence), codon_length):
                codon = sequence[j:j+codon_length]
                if codon in ["TAA", "TAG", "TGA"]:
                    orf_end = j+codon_length
                    orf_sequence = sequence[orf_start:orf_end]

                    if len(orf_sequence) >= min_orf_length:
                        start_positions.append(orf_start)
                        orf_lengths.append(len(orf_sequence))
                        orf_sequences.append(orf_sequence)
                    break
                else:
                    orf_length += codon_length

    return orf_sequences, start_positions, orf_lengths

def scoreMotif(sequence, motif):
    motif_score = []
    window_size = 13

    for i in range(len(sequence) - window_size + 1):
        window_sequence = sequence[i:i + window_size]
        window_score = 0

        for j in range(window_size):
            base = window_sequence[j]
            if base in base_idx:
                window_score += motif[base_idx[base]][j]

        motif_score.append(window_score)

    return motif_score


def identifyORFs(input_file, output_file):


    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        contig_id = None
        contig_sequence = ""

        for line in infile:
            if line.startswith('>'):
                if contig_id is not None:
                    orf_sequences, start_positions, orf_lengths = scanSeq(contig_sequence)
                    for i in range(len(orf_sequences)):
                        motif_scores = scoreMotif(orf_sequences[i], motif)
                        if max(motif_scores) > 7.25:
                            output_line = f"> {contig_id}_ORF{i+1}|Length {orf_lengths[i]}|at position {start_positions[i]}\n{orf_sequences[i]}\n"
                            outfile.write(output_line)

                contig_id = line.strip().split()[0][1:]
                contig_sequence = ""
            else:
                contig_sequence += line.strip()

        orf_sequences, start_positions, orf_lengths = scanSeq(contig_sequence)
        for i in range(len(orf_sequences)):
            motif_scores = scoreMotif(orf_sequences[i], motif)
            if max(motif_scores) > 7.25:
                output_line = f"> {contig_id}_ORF{i+1}|Length {orf_lengths[i]}|at position {start_positions[i]}\n{orf_sequences[i]}\n"
                outfile.write(output_line)

# motif score
base_idx = { 'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3 }
motif = [[.5,.5,.5,.5,0,0,0,0,0,2,-99,-99,.5],
         [0,0,0,0,0,0,0,0,0,-99,2,-99,0],
         [0,0,0,0,0,0,0,0,0,-99,-99,-99,0],
         [.5,.5,.5,.5,0,0,0,0,0,.5,-99,2,0]]

input_file = "spaceSeq.fa"
output_file = "EvanessaNkenfack.fa"
identifyORFs(input_file, output_file)
