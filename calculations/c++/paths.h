#pragma once
#include <string>

namespace filenames {
constexpr const char* UNFILTERED_GENOMES = "gisaid_hcov-19_2020_04_29_22.fasta";
constexpr const char* FILTERED_GENOMES = "Cleaned_up_genes.fasta";
constexpr const char* REFERENCE_GENES = "reference/genes.txt";
constexpr const char* REFERENCE_GENOME = "reference/REFERENCE_GENOME.fasta";
constexpr const char* TREE_UNFILTERED = "gisaid_china.MCC.trees";
constexpr const char* TREE_FILTERED = "gisaid_china_mod.MCC.trees";
constexpr const char* CALCULATED_GENOMES_DATA = "distances3.txt";
constexpr const char* EXTRA_COMPARISONS_REQUESTS = "extra_comparisons.txt";
constexpr const char* EXTRA_COMPARISONS_RESULTS = "extra_comparisons_results.txt";
constexpr const char* DISTANCES_CSV = "distances.csv";
constexpr const char* COMPARISONS_CSV = "comparisons.csv";
constexpr const char* MERGED_CSV = "genome_neigh.csv";
constexpr const char* PEPTIDES_INPUT = "netmhcpan_peptides.txt";

std::string data_path(const char* file, int sub_level = 0);
}
