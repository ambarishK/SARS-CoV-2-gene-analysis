UNFILTERED_GENOMES="gisaid_hcov-19_2020_04_29_22.fasta"
FILTERED_GENOMES="Cleaned_up_genes.fasta"
REFERENCE_GENES="reference/genes.txt"
REFERENCE_GENOME="reference/REFERENCE_GENOME.fasta"
TREE_UNFILTERED="gisaid_china.MCC.trees"
TREE_FILTERED="gisaid_china_mod.MCC.trees"
CALCULATED_GENOMES_DATA="distances3.txt"
EXTRA_COMPARISONS_REQUESTS="extra_comparisons.txt"
EXTRA_COMPARISONS_RESULTS="extra_comparisons_results.txt"
DISTANCES_CSV="distances.csv"
COMPARISONS_CSV="comparisons.csv"
MERGED_CSV="genome_neigh.csv"

import os

def data_path(file):
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", file))