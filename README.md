# DNA-Sequence-Analyzer-with-Protein-Translation

This project is a bioinformatics toolkit built in Python that reads DNA sequences from a FASTA file and performs several analyses to generate a comprehensive summary of each sequence. The output is written to a CSV file that includes:
	•	Nucleotide Counts (A, T, C, G)
	•	GC Content
	•	Reverse Complement
	•	Start Codon Positions
	•	Restriction Enzyme Cut Sites (e.g., EcoRI, BamHI, HindIII, NotI, SpeI)
	•	Protein Translation (using the standard genetic code; stops at the first stop codon)

Features
	•	Easy Input: Reads multiple sequences from a FASTA file.
	•	Comprehensive Analysis: Outputs a wide array of bioinformatics metrics.
	•	Extensible: Designed to be simple enough for beginners while providing a solid base for advanced analysis.
	•	Open Source: Everyone is welcome to view, modify, and build on this project.

Requirements
	•	Python 3.x
	•	Biopython
	•	Standard libraries: csv

1.	Prepare Your Data:
Place your FASTA file (e.g., coronavirus.fasta) in the project directory or update the file path in the script.

2. Run the script in your IDE

3.	Review the Output:
The script will generate a CSV file named dna_analysis_results_with_proteins.csv containing all your analyzed data.

How It Works

The script performs the following steps:
	•	Sequence Parsing: Uses Biopython’s SeqIO to read sequences from the FASTA file.
	•	Sequence Analysis: For each sequence, it calculates nucleotide counts, GC content, reverse complement, identifies start codon positions, finds restriction enzyme cut sites, and translates the sequence into a protein.
	•	Result Compilation: All analysis results are compiled into a dictionary for each sequence.
	•	Output: The compiled results are then written to a CSV file.

