#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <string_view>
#include <exception>
#include <unordered_map>
#include <fstream>
#include <cctype>
#include <execution>
#include <algorithm>
#include <iterator>

#include "levenshtein.h"
#include "python_printer.h"
#include "string_tools.h"
#include "paths.h"

using namespace std::string_literals;

constexpr const size_t PRE_SUF_LENGTH = 64;
constexpr const size_t MIN_GENOME_SIZE = 29000;
static_assert(PRE_SUF_LENGTH > 3);

std::string gene_to_protein(std::string_view gene) {
	static const std::unordered_map<std::string, std::string> codons {
		{"AAA", "K"}, {"AAC", "N"}, {"AAG", "K"}, {"AAU", "N"},
		{"ACA", "T"}, {"ACC", "T"}, {"ACG", "T"}, {"ACU", "T"},
		{"AGA", "R"}, {"AGC", "S"}, {"AGG", "R"}, {"AGU", "S"},
		{"AUA", "I"}, {"AUC", "I"}, {"AUG", "M"}, {"AUU", "I"},
		{"CAA", "Q"}, {"CAC", "H"}, {"CAG", "Q"}, {"CAU", "H"},
		{"CCA", "P"}, {"CCC", "P"}, {"CCG", "P"}, {"CCU", "P"},
		{"CGA", "R"}, {"CGC", "R"}, {"CGG", "R"}, {"CGU", "R"},
		{"CUA", "L"}, {"CUC", "L"}, {"CUG", "L"}, {"CUU", "L"},
		{"GAA", "E"}, {"GAC", "D"}, {"GAG", "E"}, {"GAU", "D"},
		{"GCA", "A"}, {"GCC", "A"}, {"GCG", "A"}, {"GCU", "A"},
		{"GGA", "G"}, {"GGC", "G"}, {"GGG", "G"}, {"GGU", "G"},
		{"GUA", "V"}, {"GUC", "V"}, {"GUG", "V"}, {"GUU", "V"},
		{"UAA", "_"}, {"UAC", "Y"}, {"UAG", "_"}, {"UAU", "T"},
		{"UCA", "S"}, {"UCC", "S"}, {"UCG", "S"}, {"UCU", "S"},
		{"UGA", "_"}, {"UGC", "C"}, {"UGG", "W"}, {"UGU", "C"},
		{"UUA", "L"}, {"UUC", "F"}, {"UUG", "L"}, {"UUU", "F"}
	};
	std::string protein;
	protein.reserve(gene.size() / 3 + 1);
	size_t i = 0;
	for(; i + 2 < gene.size(); i += 3) {
		std::string codon(gene.data() + i, 3);
		for(size_t j = 0; j < 3; ++j) {
			if(codon[j] == 'T') {
				codon[j] = 'U';
			}
		}
		auto protein_part = codons.find(codon);
		if(protein_part == codons.end()) {
			protein += 'X';
		} else {
			protein += protein_part->second;
		}
	}
	if(i != gene.size()) {
		protein += 'X';
	}
	return protein;
}

struct Genome;

struct Gene {
	using python_print_type = decltype(python_loader_struct<std::string, size_t, size_t, size_t>::with_names_values("name"_sl, "begin"_sl, "end"_sl,"invalid_nucleotides"_sl));
	std::string name;
	const Genome& genome;
	size_t begin;
	size_t end;
	std::string protein;
	size_t invalid_nucleotides() const;

	struct ComparisonResult {
		using python_print_type = decltype(python_loader_struct<std::string, int, int, std::vector<std::tuple<size_t, char, char>>, std::vector<std::tuple<size_t, char, char>>>::with_names_values("name"_sl, "distance"_sl, "protein_distance"_sl, "mutations"_sl, "protein_mutations"_sl));
		const Gene* compared; //may be null
		const Gene* reference; //may be null
		int distance;
		int protein_distance;
		std::vector<std::tuple<size_t, char, char>> mutations;
		std::vector<std::tuple<size_t, char, char>> protein_mutations;
		ComparisonResult() = delete;
		ComparisonResult(const Gene* compared, const Gene* reference, int distance, int protein_distance, std::vector<std::tuple<size_t, char, char>> mutations, std::vector<std::tuple<size_t, char, char>> protein_mutations) : compared(compared), reference(reference), distance(distance), protein_distance(protein_distance), mutations(std::move(mutations)), protein_mutations(std::move(protein_mutations)) {}
		ComparisonResult(const Gene* compared, const Gene* reference) : compared(compared), reference(reference), distance(0), protein_distance(0) {
			if((bool) compared == (bool) reference){
				throw std::logic_error("Compared and reference same!");
			}
		}

		static ComparisonResult from_parsed_python(const python_print_type& p, const Genome& r_compared, const Genome& r_reference);
	};

	ComparisonResult compare_to(const Gene& reference) const;

	static std::unique_ptr<Gene> from_parsed_python(const python_print_type& p, const Genome& genome);

private:
	friend struct Genome;
	Gene() = delete;
	Gene(const Gene&) = delete;
	Gene(Gene&&) = delete;
	Gene(std::string name, const Genome& genome, size_t begin, size_t end, std::string protein) : name(std::move(name)), genome(genome), begin(begin), end(end), protein(std::move(protein)) {}
};

void print_as_python(std::ostream& o, const Gene& gene) {
	print_as_python_dict(o, "name"s, gene.name, "begin"s, gene.begin, "end"s, gene.end, "invalid_nucleotides"s, gene.invalid_nucleotides());
}

void print_as_python(std::ostream& o, const Gene::ComparisonResult& cr) {
	if(cr.compared && cr.reference) {
		print_as_python_dict(o, "name"s, cr.compared->name, "distance"s, cr.distance, "protein_distance"s, cr.protein_distance, "mutations"s, cr.mutations, "protein_mutations"s, cr.protein_mutations);
	} else {
		print_as_python_dict(o, "name"s, cr.compared ? cr.compared->name : cr.reference->name, "distance"s, -1, "protein_distance"s, cr.protein_distance, "mutations"s, cr.mutations, "protein_mutations"s, cr.protein_mutations);
	}
}

template<typename It1, typename It2, typename F>
void transform_with_progress(It1 input_begin, It1 input_end, It2 output_begin, F&& function, std::ostream& output) {
	std::atomic<size_t> done = 0;
	size_t total = std::distance(input_begin, input_end);
	std::transform(std::execution::par_unseq, input_begin, input_end, output_begin, [function(std::move(function)), &done, &output, total](auto arg){
		auto result = function(arg);
		int d = done.fetch_add(1);
		if(d * 100 / total != (d + 1) * 100 / total) {
			output << std::to_string(((d + 1) * 100 / total)) + "%\n";
			output.flush();
		}
		return result;
	});
}

struct Genome {
	using python_print_type = decltype(python_loader_struct<std::string, std::vector<Gene::python_print_type>>::with_names_values("header"_sl, "genes"_sl));
	std::string header;
	std::string data;
	std::vector<std::unique_ptr<Gene>> genes;

	struct ComparisonResult {
		using python_print_type = decltype(python_loader_struct<std::string, std::string, std::vector<Gene::ComparisonResult::python_print_type>, int>::with_names_values("compared"_sl, "reference"_sl, "genes"_sl, "distance"_sl));
		const Genome& compared;
		const Genome& reference;
		std::vector<Gene::ComparisonResult> genes;
		int distance;
		ComparisonResult() = delete;
		ComparisonResult(const ComparisonResult&) = delete;
		ComparisonResult(ComparisonResult&&) = default;
		ComparisonResult(const Genome& compared, const Genome& reference, std::vector<Gene::ComparisonResult> genes, int distance) : compared(compared), reference(reference), genes(std::move(genes)), distance(distance) {}
		static ComparisonResult from_parsed_python(const python_print_type& p, const Genome& compared, const Genome& reference) {
			if(compared.header != std::get<0>(p.as_tuple()) || reference.header != std::get<1>(p.as_tuple())) {
				throw std::runtime_error("Invalid parsed comparison.");
			}
			std::vector<Gene::ComparisonResult> gcr;
			for(const Gene::ComparisonResult::python_print_type& gene : std::get<2>(p.as_tuple())) {
				gcr.push_back(Gene::ComparisonResult::from_parsed_python(gene, compared, reference));
			}
			return ComparisonResult(compared, reference, std::move(gcr), std::get<3>(p.as_tuple()));
		}
	};

	int distance_to(const Genome& reference, const std::vector<Gene::ComparisonResult>& genes_compared) const {
		auto sort_genes = [](const std::vector<std::unique_ptr<Gene>>& genes) {
			std::vector<const Gene*> output;
			for(const std::unique_ptr<Gene>& gene: genes) {
				output.push_back(gene.get());
			}
			sort(output.begin(), output.end(), [](const Gene* g1, const Gene* g2){return g1->begin < g2->begin;});
			return output;
		};

		std::array<std::vector<const Gene*>, 2> sorted_genes{sort_genes(genes), sort_genes(reference.genes)};

		for(int j = 0; j < 2; ++j) {
			for(size_t i = 1; i < sorted_genes[j].size(); ++i) {
				if(sorted_genes[j][i]->begin <= sorted_genes[j][i - 1]->end) {
					return -2;
				}
			}
		}

		if(genes.size() != reference.genes.size()) {
			return -4;
		}

		for(size_t i = 0; i < genes.size(); ++i) {
			if(genes[i]->name != reference.genes[i]->name) {
				return -3;
			}
		}

		std::vector<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>> fragments;
		size_t ends[2] = {0, 0};

		for(size_t i = 0; i < genes.size(); ++i) {
			fragments.emplace_back(std::make_pair(ends[0], sorted_genes[0][i]->begin), std::make_pair(ends[1], sorted_genes[1][i]->begin));
			ends[0] = sorted_genes[0][i]->end;
			ends[1] = sorted_genes[1][i]->end;
		}

		fragments.emplace_back(std::make_pair(ends[0], data.size()), std::make_pair(ends[1], reference.data.size()));

		int score = 0;
		for(const auto&[f, s] : fragments) {
			score += lev_edit_distance_int<int>(data.data() + f.first, f.second - f.first, reference.data.data() + s.first, s.second - s.first);
		}

		for(const Gene::ComparisonResult& cr : genes_compared) {
			score += cr.distance;
		}

		return score;
	}

	static std::unique_ptr<Genome> load_reference(std::istream& ref_genes_file, std::istream& ref_genome_file) {
		std::unique_ptr<Genome> genome;
		{
			std::string header;
			std::string data;
			size_t total_length = 0;
			std::vector<std::string> data_lines;
			std::string line;
			while(std::getline(ref_genome_file, line)) {
				rstrip(line);
				if(!line.empty() && line[0] == '>') {
					header = std::move(line);
					continue;
				}
				total_length += line.size();
				data_lines.emplace_back(std::move(line));
			}
			data = concatenate_strings(data_lines, total_length);
			genome = std::unique_ptr<Genome>(new Genome(std::move(header), std::move(data), std::vector<std::unique_ptr<Gene>>()));
		}

		{
			std::string gene_name;
			size_t gene_begin;
			size_t gene_end;
			while(ref_genes_file >> gene_name >> gene_begin >> gene_end) {
				--gene_begin;
				if(gene_end <= gene_begin || gene_end >= genome->data.size()) {
					throw std::runtime_error("Invalid boundaries in line starting with " + gene_name);
				}
				std::string protein = gene_to_protein(std::string_view(genome->data.data() + gene_begin, gene_end - gene_begin));
				genome->genes.emplace_back(new Gene(std::move(gene_name), *genome, gene_begin, gene_end, std::move(protein)));
				gene_name.clear();
			}
			if(!gene_name.empty()) {
				throw std::runtime_error("Invalid line starting with " + gene_name);
			}
		}

		return genome;
	}

	static std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>> make_genome(std::string header, std::string data, const Genome& reference, size_t tolerance = 128) {
		std::unique_ptr<Genome> genome = std::unique_ptr<Genome>(new Genome(std::move(header), std::move(data), std::vector<std::unique_ptr<Gene>>()));
		std::vector<Gene::ComparisonResult> genes_compared;
		for(const std::unique_ptr<Gene>& ref_gene_ptr : reference.genes) {
			const Gene& ref_gene = *ref_gene_ptr;

			size_t boundary_left = ref_gene.begin > tolerance ? ref_gene.begin - tolerance : 0ULL;
			if(boundary_left >= genome->data.size()) {
				genes_compared.emplace_back(nullptr, ref_gene_ptr.get());
				continue;
			}
			size_t boundary_right = std::min(genome->data.size(), ref_gene.end + tolerance);

			std::string_view ref_prefix(reference.data.data() + ref_gene.begin, PRE_SUF_LENGTH);
			std::string_view ref_suffix(reference.data.data() + ref_gene.end - PRE_SUF_LENGTH, PRE_SUF_LENGTH);
			std::string_view seq(genome->data.data() + boundary_left, boundary_right - boundary_left);
			size_t prefix = std::string::npos, suffix = std::string::npos;
			int prefix_distance = INT_MAX, suffix_distance = INT_MAX;
			{
				for(size_t i = 0; i + PRE_SUF_LENGTH - 1 < seq.size(); ++i) {
					if(seq[i] == 'A' && seq[i + 1] == 'T' && seq[i + 2] == 'G') {
						int distance = lev_edit_distance_int<int>(seq.data() + i, PRE_SUF_LENGTH, ref_prefix.data(), PRE_SUF_LENGTH, prefix_distance);
						if(distance < prefix_distance) {
							prefix = i;
							prefix_distance = distance;
						}
					}
				}
			}
			{
				for(size_t i = PRE_SUF_LENGTH - 3; i + 2 < seq.size(); ++i) {
					if(seq[i] == 'T' && ((seq[i + 1] == 'A' && seq[i + 2] == 'G') || (seq[i + 1] == 'A' && seq[i + 2] == 'A') || (seq[i + 1] == 'G' && seq[i + 2] == 'A'))) {
						int distance = lev_edit_distance_int<int>(seq.data() + i + 3 - PRE_SUF_LENGTH, PRE_SUF_LENGTH, ref_suffix.data(), PRE_SUF_LENGTH, suffix_distance);
						if(distance < suffix_distance) {
							suffix = i + 3;
							suffix_distance = distance;
						}
					}
				}
			}
			if(prefix == std::string::npos || suffix == std::string::npos) {
				genes_compared.emplace_back(nullptr, ref_gene_ptr.get());
				continue;
			}
			if(prefix > suffix) {
				genes_compared.emplace_back(nullptr, ref_gene_ptr.get());
				continue;
			}
			std::string_view found_gene(seq.data() + prefix, suffix - prefix);

			genome->genes.emplace_back(new Gene(
				ref_gene.name,
				*genome,
				boundary_left + prefix,
				boundary_left + suffix,
				gene_to_protein(found_gene)
			));
			genes_compared.push_back(genome->genes.back()->compare_to(ref_gene));
		}
		std::unique_ptr<Genome::ComparisonResult> cr = std::make_unique<Genome::ComparisonResult>(*genome, reference, std::move(genes_compared), genome->distance_to(reference, genes_compared));
		return std::make_pair(std::move(genome), std::move(cr));
	}

	static std::vector<std::pair<std::string, std::string>> load_genomes_raw(std::istream& input, std::ostream& info) {
		std::vector<std::pair<std::string, std::string>> genomes_raw; //header, data
		std::vector<std::string> data;
		size_t genome_length = 0;
		std::string header;
		std::string line;
		while(std::getline(input, line)) {
			rstrip(line);
			if(!line.empty() && line[0] == '>') {
				if(!header.empty() && genome_length > MIN_GENOME_SIZE) {
					genomes_raw.emplace_back(std::move(header), concatenate_strings(data, genome_length));
				}
				genome_length = 0;
				data.clear();
				header = std::move(line);
			} else {
				genome_length += line.size();
				for(char& c : line) {
					c = toupper(c);
				}
				data.emplace_back(std::move(line));
			}
		}
		if(!header.empty() && genome_length > MIN_GENOME_SIZE) {
			genomes_raw.emplace_back(std::move(header), concatenate_strings(data, genome_length));
		}
		info << "We have in total " << genomes_raw.size() << " genomes collected from the fasta file" << std::endl;
		return genomes_raw;
	}

	static std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> load_genomes(std::istream& input, const Genome& reference, std::ostream& info) {
		std::vector<std::pair<std::string, std::string>> genomes_raw = load_genomes_raw(input, info); //header, data
		std::vector<std::pair<Genome*, Genome::ComparisonResult*>> genomes_p_r(genomes_raw.size(), std::make_pair<Genome*, Genome::ComparisonResult*>(nullptr, nullptr));
		try {
			transform_with_progress(genomes_raw.begin(), genomes_raw.end(), genomes_p_r.begin(), [&](std::pair<std::string, std::string>& genome_raw){
				auto genome = make_genome(std::move(genome_raw.first), std::move(genome_raw.second), reference);
				return std::make_pair(genome.first.release(), genome.second.release());
			}, info);
		} catch (...) {
			for(const auto&[g, cr] : genomes_p_r) {
				delete g;
				delete cr;
			}
			throw;
		}

		std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes;
		genomes.reserve(genomes_p_r.size());
		for(const auto&[g, cr] : genomes_p_r) {
			genomes.emplace_back(g, cr);
		}

		size_t errors = 0;
		for(const auto&[genome, cr] : genomes) {
			if(cr->distance < 0) {
				++errors;
			}
		}
		info << errors << " distance errors" << std::endl;
		info << "In total we have processed " << genomes.size() << " genomes" << std::endl;
		return genomes;
	}

	Genome::ComparisonResult compare_to(const Genome& reference) const {
		std::unordered_map<std::string, Gene::ComparisonResult> genes_compared;
		for(const auto& gene: genes) {
			bool found = false;
			for(const auto& gene_reference: reference.genes) {
				if(gene_reference->name == gene->name) {
					genes_compared.emplace(gene->name, gene->compare_to(*gene_reference));
					found = true;
					break;
				}
			}
			if(!found) {
				genes_compared.emplace(gene->name, Gene::ComparisonResult(gene.get(), nullptr));
			}
		}
		for(const auto& gene_reference: reference.genes) {
			if(genes_compared.find(gene_reference->name) == genes_compared.end()) {
				genes_compared.emplace(gene_reference->name, Gene::ComparisonResult(nullptr, gene_reference.get()));
			}
		}
		std::vector<Gene::ComparisonResult> genes_compared_vec;
		genes_compared_vec.reserve(genes_compared.size());
		for(auto it = genes_compared.begin(); it != genes_compared.end(); ++it) {
			genes_compared_vec.push_back(std::move(it->second));
		}
		int distance = distance_to(reference, genes_compared_vec);
		return Genome::ComparisonResult(*this, reference, std::move(genes_compared_vec), distance);
	}

	static std::unique_ptr<Genome> from_parsed_python(const python_print_type& p, std::string data) {
		auto genome = std::unique_ptr<Genome>(new Genome(std::get<0>(p.as_tuple()), std::move(data), {}));
		for(const Gene::python_print_type& gene : std::get<1>(p.as_tuple())) {
			genome->genes.push_back(Gene::from_parsed_python(gene, *genome));
		}
		return genome;
	}

private:
	Genome() = delete;
	Genome(const Genome&) = delete;
	Genome(Genome&&) = delete;
	Genome(std::string header, std::string data, std::vector<std::unique_ptr<Gene>> genes) : header(std::move(header)), data(std::move(data)), genes(std::move(genes)) {}
};

std::unique_ptr<Gene> Gene::from_parsed_python(const Gene::python_print_type& p, const Genome& genome) {
	auto t = p.as_tuple();
	if(std::get<2>(t) < std::get<1>(t) || std::get<2>(t) > genome.data.size()) {
		throw std::runtime_error("Invalid gene.");
	}
	return std::unique_ptr<Gene>(new Gene(std::get<0>(t), genome, std::get<1>(t), std::get<2>(t), gene_to_protein(std::string_view(genome.data.data() + std::get<1>(t), std::get<2>(t) - std::get<1>(t)))));
}

Gene::ComparisonResult Gene::ComparisonResult::from_parsed_python(const Gene::ComparisonResult::python_print_type& p, const Genome& r_compared, const Genome& r_reference) {
	const Gene* compared = nullptr;
	const Gene* reference = nullptr;
	for(const std::unique_ptr<Gene>& g : r_compared.genes) {
		if(g->name == std::get<0>(p.as_tuple())) {
			compared = g.get();
			break;
		}
	}
	for(const std::unique_ptr<Gene>& g : r_reference.genes) {
		if(g->name == std::get<0>(p.as_tuple())) {
			reference = g.get();
			break;
		}
	}
	return Gene::ComparisonResult(compared, reference, std::get<1>(p.as_tuple()), std::get<2>(p.as_tuple()), std::get<3>(p.as_tuple()), std::get<4>(p.as_tuple()));
}

size_t Gene::invalid_nucleotides() const {
	size_t n = 0;
	for(size_t i = begin; i < end; ++i) {
		char c = genome.data[i];
		if(c != 'A' && c != 'C' && c != 'T' && c != 'G') {
			++n;
		}
	}
	return n;
}

Gene::ComparisonResult Gene::compare_to(const Gene& reference) const {
	return ComparisonResult(
		this,
		&reference,
		lev_edit_distance_int<int>(genome.data.data() + begin, end - begin, reference.genome.data.data() + reference.begin, reference.end - reference.begin),
		lev_edit_distance_int<int>(protein.data(), protein.size(), reference.protein.data(), reference.protein.size()),
		get_differences(reference.genome.data.data() + reference.begin, genome.data.data() + begin, std::min(end - begin, reference.end - reference.begin)),
		get_differences(reference.protein.data(), protein.data(), std::min(protein.size(), reference.protein.size()))
	);
}

void print_as_python(std::ostream& o, const Genome& genome) {
	print_as_python_dict(o, "header"s, genome.header, "genes"s, genome.genes);
}

void print_as_python(std::ostream& o, const Genome::ComparisonResult& cr) {
	print_as_python_dict(o, "compared"s, cr.compared.header, "reference"s, cr.reference.header, "genes"s, cr.genes, "distance"s, cr.distance);
}

std::ifstream open_file_i(const std::string& filename) {
	std::ifstream f(filename);
	if(!f.is_open()) {
		throw std::runtime_error("Failed to open file " + filename);
	}
	return f;
}

std::ofstream open_file_o(const std::string& filename) {
	std::ofstream f(filename);
	if(!f.is_open()) {
		throw std::runtime_error("Failed to open file " + filename);
	}
	return f;
}

void calculate_genomes(const Genome& reference, const std::string& genomes_filename, const std::string& output_filename) {
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes;
	{
		auto input = open_file_i(genomes_filename);
		genomes = Genome::load_genomes(input, reference, std::cout);
	}
	{
		auto output = open_file_o(output_filename);
		output << genomes.size() << '\n';
		for(const auto& e : genomes) {
			print_as_python(output, e);
			output << '\n';
		}
	}
}

void extra_comparisons(const Genome& reference, const std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>>& genomes, const std::string& extra_comparisons_filename, const std::string& output_filename) {
	std::vector<std::pair<std::string, std::string>> extra_comparison_tasks;
	{
		auto input = open_file_i(extra_comparisons_filename);
		size_t tasks;
		input >> tasks;
		if(input.fail() || input.bad()) {
			std::cout << "Failed to read from " << extra_comparisons_filename << std::endl;
			return;
		}
		{
			std::string a;
			std::getline(input, a);
		}
		for(size_t i = 0; i < tasks; ++i) {
			std::string a, b;
			std::getline(input, a);
			std::getline(input, b);
			rstrip(a);
			rstrip(b);
			if(!input) {
				std::cout << "Failed to read from " << extra_comparisons_filename << std::endl;
				return;
			}
			extra_comparison_tasks.emplace_back(std::move(a), std::move(b));
		}
	}
	std::vector<std::shared_ptr<Genome::ComparisonResult>> extra_comparison_results(extra_comparison_tasks.size(), nullptr);
	{
		std::unordered_map<std::string, const Genome*> genomes_map;
		genomes_map.emplace(reference.header, &reference);
		for(const auto&[genome, cr] : genomes) {
			if(!genomes_map.emplace(genome->header, genome.get()).second) {
				throw std::runtime_error("Header" + genome->header + " specified twice");
			}
		}
		transform_with_progress(extra_comparison_tasks.begin(), extra_comparison_tasks.end(), extra_comparison_results.begin(), [&genomes_map](const std::pair<std::string, std::string>& t){
			auto r = genomes_map.find(t.first);
			auto c = genomes_map.find(t.second);
			if(r == genomes_map.end()) {
				throw std::runtime_error("Header " + t.first + " not found");
			}
			if(c == genomes_map.end()) {
				throw std::runtime_error("Header " + t.second + " not found");
			}
			return std::shared_ptr<Genome::ComparisonResult>(new Genome::ComparisonResult(c->second->compare_to(*r->second)));
		}, std::cout);
	}
	size_t errors = 0;
	for(const auto& cr : extra_comparison_results) {
		if(cr->distance < 0) {
			++errors;
		}
	}
	std::cout << errors << " distance errors in " << extra_comparison_results.size() << " comparisons." << std::endl;
	{
		auto output = open_file_o(output_filename);
		if(!output.is_open()) {
			std::cout << "Failed to open file " << output_filename << std::endl;
			return;
		}
		output << extra_comparison_results.size() << '\n';
		for(const auto& e : extra_comparison_results) {
			print_as_python(output, e);
			output << '\n';
		}
	}
}

std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> load_genomes_from_files(const Genome& reference, const std::string& genomes_filename, const std::string& input_filename) {
	std::unordered_map<std::string, std::string> genomes_raw_map; //header -> data
	{
		std::vector<std::pair<std::string, std::string>> genomes_raw; //header, data
		auto input = open_file_i(genomes_filename);
		genomes_raw = Genome::load_genomes_raw(input, std::cout);
		for(const auto&[h, v] : genomes_raw) {
			genomes_raw_map.emplace(h, std::move(v));
		}
	}
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> output;
	{
		using DataT = std::pair<Genome::python_print_type, Genome::ComparisonResult::python_print_type>;
		auto input = open_file_i(input_filename);
		size_t size;
		input >> size;
		if(input.bad() || input.fail()) {
			throw std::runtime_error("Error reading from file " + input_filename);
		}
		std::vector<DataT> genomes_extra_data;
		genomes_extra_data.reserve(size);
		for(size_t i = 0; i < size; ++i) {
			genomes_extra_data.push_back(python_loader<DataT>::load(input));
		}
		for(const auto& g : genomes_extra_data) {
			auto genome = Genome::from_parsed_python(g.first, std::move(genomes_raw_map.at(std::get<0>(g.first.as_tuple()))));
			auto comparison = std::make_unique<Genome::ComparisonResult>(Genome::ComparisonResult::from_parsed_python(g.second, *genome, reference));
			output.emplace_back(std::move(genome), std::move(comparison));
		}
	}
	return output;
}

int main(int argc, char* argv[]){
	using namespace filenames;
	if(argc < 2 || (argv[1] != "distances"s && argv[1] != "comparisons"s)) {
		std::cout << "USAGE:\nfindgenemutations <distances|comparisons>\n";
		return -1;
	}
	std::unique_ptr<Genome> reference;
	{
		auto ref_info_file = open_file_i(data_path(REFERENCE_GENES));
		auto ref_genome_file = open_file_i(data_path(REFERENCE_GENOME));
		reference = Genome::load_reference(ref_info_file, ref_genome_file);
	}
	if(argv[1] == "distances"s) {
		calculate_genomes(*reference, data_path(FILTERED_GENOMES), data_path(CALCULATED_GENOMES_DATA));
	} else {
		std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes = load_genomes_from_files(*reference, data_path(FILTERED_GENOMES), data_path(CALCULATED_GENOMES_DATA));
		std::cout << "Genomes loaded" << std::endl;
		extra_comparisons(*reference, genomes, data_path(EXTRA_COMPARISONS_REQUESTS), data_path(EXTRA_COMPARISONS_RESULTS));
	}
}
