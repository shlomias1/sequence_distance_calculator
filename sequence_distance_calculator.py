from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
import math

def fetch_sequence(accession):
    """Fetch sequence from NCBI using Biopython."""
    Entrez.email = "shlomiasi1@gmail.com"
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)

def align_sequences(seq1, seq2):
    """Align two sequences using PairwiseAligner with basic settings."""
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # Use local alignment
    aligner.open_gap_score = -10  # Penalty for opening a gap
    aligner.extend_gap_score = -0.5  # Penalty for extending a gap
    alignments = list(aligner.align(seq1, seq2))
    best_alignment = alignments[0]  # Choose the best alignment
    aligned_seq1 = best_alignment.target  # Aligned sequence from target
    aligned_seq2 = best_alignment.query  # Aligned sequence from query
    return aligned_seq1, aligned_seq2

def preprocess_sequences(seq1, seq2):
    """Remove gaps ('-') from aligned sequences."""
    seq1_cleaned = "".join(base1 for base1, base2 in zip(seq1, seq2) if base1 != '-' and base2 != '-')
    seq2_cleaned = "".join(base2 for base1, base2 in zip(seq1, seq2) if base1 != '-' and base2 != '-')
    return seq1_cleaned, seq2_cleaned

def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length.")
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

def jukes_cantor_distance(seq1, seq2):
    """Calculate Jukes-Cantor distance based on Hamming distance."""
    seq_length = len(seq1)
    dH = hamming_distance(seq1, seq2) / seq_length  # Proportion of differences
    if dH >= 0.75:
        raise ValueError("Jukes-Cantor model is not valid when dH >= 0.75.")
    dJC = -3/4 * math.log(1 - (4/3) * dH)
    return dJC

def calculate_transitions_and_transversions(seq1, seq2):
    """Calculate the number of transitions (P) and transversions (Q)."""
    transitions = 0
    transversions = 0

    purines = {'A', 'G'}  # A and G are purines
    pyrimidines = {'C', 'T'}  # C and T are pyrimidines

    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:  # Mutation occurred
            if {base1, base2} <= purines or {base1, base2} <= pyrimidines:
                transitions += 1
            else:
                transversions += 1

    return transitions, transversions

def kimura_distance(seq1, seq2):
    """Calculate Kimura distance between two sequences."""
    seq_length = len(seq1)
    transitions, transversions = calculate_transitions_and_transversions(seq1, seq2)

    P = transitions / seq_length  # Proportion of transitions
    Q = transversions / seq_length  # Proportion of transversions

    # Check for invalid values
    if 1 - 2 * P - Q <= 0 or 1 - 2 * Q <= 0:
        raise ValueError("Kimura model is not valid for these sequences.")

    dK = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
    return dK

def main():
    """Main function to fetch, align, preprocess, and calculate distances."""
    # Accessions for human and orangutan sequences
    human_accession = "D38112"
    orangutan_accession = "D38115"

    print("Fetching sequences...")
    human_seq = fetch_sequence(human_accession)
    orangutan_seq = fetch_sequence(orangutan_accession)

    # Trim sequences to first 500 bases for alignment
    human_seq = human_seq[:500]
    orangutan_seq = orangutan_seq[:500]

    print(f"Trimmed human_seq: {human_seq[:50]}... (length: {len(human_seq)})")
    print(f"Trimmed orangutan_seq: {orangutan_seq[:50]}... (length: {len(orangutan_seq)})")

    print("Aligning sequences...")
    aligned_human, aligned_orangutan = align_sequences(human_seq, orangutan_seq)

    print("Preprocessing sequences...")
    clean_human, clean_orangutan = preprocess_sequences(aligned_human, aligned_orangutan)

    print("Calculating distances...")
    try:
        dH = hamming_distance(clean_human, clean_orangutan)
        dJC = jukes_cantor_distance(clean_human, clean_orangutan)
        dK = kimura_distance(clean_human, clean_orangutan)
        print(f"Hamming Distance (dH): {dH}")
        print(f"Jukes-Cantor Distance (dJC): {dJC}")
        print(f"Kimura Distance (dK): {dK}")
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()