# **DNA Sequence Distance Calculator**

This Python script calculates evolutionary distances between two DNA sequences using multiple metrics, including Hamming distance, Jukes-Cantor distance, and Kimura distance. The script also aligns sequences, removes gaps, and fetches DNA data directly from the NCBI database using Biopython.

---

## **Features**
- **Sequence Retrieval**: Fetches DNA sequences by accession numbers from the NCBI database.
- **Sequence Alignment**: Performs local alignment using Biopython's `PairwiseAligner`.
- **Distance Metrics**:
  - **Hamming Distance**: Counts the number of mismatched positions between two aligned sequences.
  - **Jukes-Cantor Distance**: Adjusts the Hamming distance to account for multiple mutations at the same position.
  - **Kimura Distance**: Differentiates between transitions (purine-to-purine or pyrimidine-to-pyrimidine changes) and transversions (purine-to-pyrimidine or vice versa).
- **Preprocessing**: Cleans aligned sequences by removing gaps for accurate distance calculations.

---

## **Installation**
1. Install Python (version 3.8 or higher recommended).
2. Install Biopython:
   ```bash
   pip install biopython
   ```

---

## **Usage**
1. Clone or download this repository.
2. Run the script:
   ```bash
   python sequence_distance_calculator.py
   ```

### **Inputs**
- The script fetches DNA sequences using hardcoded accession numbers:
  - **Human**: `D38112`
  - **Orangutan**: `D38115`

### **Outputs**
The script prints:
- Trimmed sequences (first 500 bases).
- Aligned sequences (post-alignment).
- Calculated distances:
  - Hamming Distance (dH).
  - Jukes-Cantor Distance (dJC).
  - Kimura Distance (dK).

---

## **Sample Output**
When running the script, you can expect output similar to the following:
```
Fetching sequences...
Trimmed human_seq: GTTTATGTAGCTTACCTCCTCAAAGCAATACACTGAAAATGTTTAGACGG... (length: 500)
Trimmed orangutan_seq: GTTTATGTAGCTTATTCCATCCAAAGCAATACACTGAAAATGTCTCGATG... (length: 500)
Aligning sequences...
Preprocessing sequences...
Calculating distances...
Hamming Distance (dH): 287
Jukes-Cantor Distance (dJC): 1.0871919086191528
Kimura Distance (dK): 1.08728804256795
```

---

## **Explanation of Distance Metrics**
### **1. Hamming Distance**
The Hamming distance counts the number of mismatches between two sequences of equal length:
\[
d_H = \sum_{i=1}^{n} I(s1_i \neq s2_i)
\]
Where \( I \) is an indicator function, and \( n \) is the sequence length. This metric does not account for the possibility of multiple mutations at a single position.

---

### **2. Jukes-Cantor Distance**
The Jukes-Cantor model extends the Hamming distance by modeling the probability of multiple substitutions at the same position. The formula is:
\[
d_{JC} = -\frac{3}{4} \ln(1 - \frac{4}{3} d_H)
\]
Where \( d_H \) is the proportion of mismatched positions. This model assumes:
- Equal substitution rates between all nucleotides.
- Independent substitutions at each position.

**Limitation**: The Jukes-Cantor model becomes invalid if \( d_H \geq 0.75 \).

---

### **3. Kimura Distance**
The Kimura model improves on the Jukes-Cantor model by distinguishing between:
- **Transitions**: Substitutions between purines (A ↔ G) or pyrimidines (C ↔ T).
- **Transversions**: Substitutions between purines and pyrimidines.

The formula is:
\[
d_K = -\frac{1}{2} \ln(1 - 2P - Q) - \frac{1}{4} \ln(1 - 2Q)
\]
Where:
- \( P \): Proportion of transitions.
- \( Q \): Proportion of transversions.

**Limitation**: The model is invalid if:
\[
1 - 2P - Q \leq 0 \quad \text{or} \quad 1 - 2Q \leq 0
\]

---

## **How It Works**
1. **Fetch Sequences**:
   - Retrieves DNA sequences using NCBI accession numbers.
   - Trims sequences to the first 500 bases for performance.

2. **Align Sequences**:
   - Aligns the two sequences locally using `PairwiseAligner`.

3. **Preprocess Sequences**:
   - Removes gaps (`-`) from the aligned sequences.

4. **Calculate Distances**:
   - **Hamming Distance**: Counts mismatches.
   - **Jukes-Cantor Distance**: Models evolutionary distance accounting for multiple mutations.
   - **Kimura Distance**: Further refines distances by distinguishing between transitions and transversions.

---

## **Dependencies**
- **Python**: Version 3.8 or higher.
- **Biopython**: For sequence retrieval and alignment.

---

## **Extending the Script**
- **Custom Accessions**: Modify `human_accession` and `orangutan_accession` in the `main()` function to fetch different sequences.
- **Adjust Alignment Settings**: Tweak gap penalties (`open_gap_score` and `extend_gap_score`) in `align_sequences()` for specific use cases.
- **Longer Sequences**: Remove or increase the trimming limit (currently set to 500 bases).

---

## **Known Limitations**
1. The script trims sequences to 500 bases for computational efficiency. For longer sequences, adjust the trimming or optimize alignment parameters.
2. The Kimura model and Jukes-Cantor model have constraints. If sequences have high divergence, the models may become invalid, raising errors.

---

## **Contributing**
Feel free to submit issues or pull requests to enhance functionality, fix bugs, or add new features.

---

## **License**
This project is licensed under the MIT License. See the LICENSE file for details.
