import csv
from Bio import SeqIO
from Bio.Seq import Seq

# Define restriction enzymes and their recognition sites
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "SpeI": "ACTAGT"
}

def analyze_sequence(dna_seq, sequence_id):
    
    # Ensure uppercase for consistency
    dna_seq = dna_seq.upper()

    # Validate sequence (only A, T, C, G allowed)
    if not all(nuc in "ATCG" for nuc in dna_seq):
        print(f"Skipping {sequence_id}: Invalid characters found in sequence.")
        return None

    # Nucleotide frequencies
    nucleotide_counts = {nuc: dna_seq.count(nuc) for nuc in "ATCG"}

    # GC Content Calculation
    gc_content = round(((nucleotide_counts["G"] + nucleotide_counts["C"]) / len(dna_seq)) * 100, 2)

    # Reverse Complement Calculation
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_complement = "".join(complement[nuc] for nuc in reversed(dna_seq))

    # Find Start Codon Positions ("ATG")
    start_codon_positions = [i for i in range(len(dna_seq) - 2) if dna_seq[i:i+3] == "ATG"]

    # Find Restriction Enzyme Cut Sites
    cut_sites = {}
    for enzyme, site in RESTRICTION_ENZYMES.items():
        positions = [i for i in range(len(dna_seq) - len(site) + 1) if dna_seq[i:i+len(site)] == site]
        if positions:
            cut_sites[enzyme] = positions

    # Translate the DNA sequence into a protein sequence (stops at the first stop codon)
    protein_translation = str(Seq(dna_seq).translate(to_stop=True))

    # Compile results into a dictionary
    result = {
        "Sequence_ID": sequence_id,
        "Nucleotide_Counts": str(nucleotide_counts),
        "GC_Content": gc_content,
        "Reverse_Complement": reverse_complement,
        "Start_Codon_Positions": str(start_codon_positions),
        "Cut_Sites": str(cut_sites),
        "Protein_Translation": protein_translation
    }
    return result

def read_and_analyze_fasta(file_path):
    
    results = []
    for record in SeqIO.parse(file_path, "fasta"):
        analysis = analyze_sequence(str(record.seq), record.id)
        if analysis is not None:
            results.append(analysis)
    return results

def write_to_csv(results, output_csv):
   
    columns = ["Sequence_ID", "Nucleotide_Counts", "GC_Content", "Reverse_Complement",
               "Start_Codon_Positions", "Cut_Sites", "Protein_Translation"]
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"Results written to {output_csv}")

if __name__ == "__main__":
    fasta_file = "/Users/huntereppley/Documents/Bioinformatics/coronavirus.fasta"  # Input FASTA file
    output_csv = "dna_analysis_results_with_proteins.csv"  # Output CSV file name

    # Analyze the sequences from the FASTA file
    results = read_and_analyze_fasta(fasta_file)
    
    # Write the combined results to a CSV file
    write_to_csv(results, output_csv)
