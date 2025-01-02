# **Overview**
This project involves creating a Python program to analyze DNA sequences in FASTA format and identify Open Reading Frames (ORFs). An ORF is a sequence of DNA that starts with a start codon and ends with a stop codon, potentially representing a gene.

The program supports the following functionalities:
- Identification of ORFs across all six reading frames.
- Reporting of ORFs meeting customizable criteria (e.g., minimum length, start/stop codons).
- Writing formatted output to a specified file.

---

# **Deliverables**
## **Python Scripts**
- `findORFs.py`: Main script to find ORFs.
- `sequenceAnalysis.py`: Contains supporting modules like `OrfFinder` and `FastAreader`.

## **Output Files**
- Example: `tass2ORFdata-ATG-100.txt`.

---

# **Features**
### **Input**
- Reads DNA sequences in FASTA format via STDIN.
- Example: `python findORFs.py < input.fa`.

### **Output**
- Writes formatted data to STDOUT.
- Example: `tass2ORFdata-ATG-100.txt`.

### **ORF Criteria**
- Minimum gene size: 100 nucleotides (start and stop codons included).
- Start codons: `ATG`.
- Stop codons: `TAA`, `TGA`, `TAG`.
- Reports only the largest ORF within each region.

### **Reading Frames**
- Analyzes all six reading frames (three forward, three reverse).
