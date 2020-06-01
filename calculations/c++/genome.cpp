#include "genome.h"
#include "util.h"

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

using namespace std::string_literals;

void print_as_python(std::ostream& o, const Gene& gene) {
	print_as_python_dict(o, "name"s, gene.name, "begin"s, gene.begin, "end"s, gene.end, "invalid_nucleotides"s, gene.invalid_nucleotides());
}

void print_as_python(std::ostream& o, const Gene::ComparisonResult& cr) {
	print_as_python_dict(o, "name"s, cr.compared.name, "mutations"s, cr.mutations, "protein_mutations"s, cr.protein_mutations);
}


Gene::ComparisonResult::ComparisonResult(const Gene& compared, const Gene& reference, std::vector<EditOperation> mutations, std::vector<EditOperation> protein_mutations) : compared(compared), reference(reference), mutations(std::move(mutations)), protein_mutations(std::move(protein_mutations)) {}

Gene::Gene(std::string name, const Genome& genome, size_t begin, size_t end, std::string protein) : name(std::move(name)), genome(genome), begin(begin), end(end), protein(std::move(protein)) {}

Genome::ComparisonResult::ComparisonResult(const Genome& compared, const Genome& reference, std::vector<Gene::ComparisonResult> genes, int distance) : compared(compared), reference(reference), genes(std::move(genes)), distance(distance) {}
Genome::ComparisonResult Genome::ComparisonResult::from_parsed_python(const python_print_type& p, const Genome& compared, const Genome& reference) {
	if(compared.header != std::get<0>(p.as_tuple()) || reference.header != std::get<1>(p.as_tuple())) {
		throw std::runtime_error("Invalid parsed comparison.");
	}
	std::vector<Gene::ComparisonResult> gcr;
	for(const Gene::ComparisonResult::python_print_type& gene : std::get<2>(p.as_tuple())) {
		gcr.push_back(Gene::ComparisonResult::from_parsed_python(gene, compared, reference));
	}
	return ComparisonResult(compared, reference, std::move(gcr), std::get<3>(p.as_tuple()));
}

int Genome::distance_to(const Genome& reference, const std::vector<Gene::ComparisonResult>& genes_compared) const {
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
		try {
			score += edit_distance_int(data.data() + f.first, f.second - f.first, reference.data.data() + s.first, s.second - s.first)->size();
		} catch (const std::runtime_error&) {
			return -5;
		}
	}

	for(const Gene::ComparisonResult& cr : genes_compared) {
		score += cr.mutations.size();
	}

	return score;
}

std::unique_ptr<Genome> Genome::load_reference(std::istream& ref_genes_file, std::istream& ref_genome_file) {
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

std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>> Genome::make_genome(std::string header, std::string data, const Genome& reference, size_t tolerance) {
	std::unique_ptr<Genome> genome = std::unique_ptr<Genome>(new Genome(std::move(header), std::move(data), std::vector<std::unique_ptr<Gene>>()));
	std::vector<Gene::ComparisonResult> genes_compared;
	for(const std::unique_ptr<Gene>& ref_gene_ptr : reference.genes) {
		const Gene& ref_gene = *ref_gene_ptr;

		size_t boundary_left = ref_gene.begin > tolerance ? ref_gene.begin - tolerance : 0ULL;
		if(boundary_left >= genome->data.size()) {
			continue;
		}
		size_t boundary_right = std::min(genome->data.size(), ref_gene.end + tolerance);

		std::string_view ref_prefix(reference.data.data() + ref_gene.begin, PRE_SUF_LENGTH);
		std::string_view ref_suffix(reference.data.data() + ref_gene.end - PRE_SUF_LENGTH, PRE_SUF_LENGTH);
		std::string_view seq(genome->data.data() + boundary_left, boundary_right - boundary_left);
		size_t prefix = std::string::npos, suffix = std::string::npos;
		int prefix_distance = std::numeric_limits<int>::max(), suffix_distance = std::numeric_limits<int>::max();
		{
			for(size_t i = 0; i + PRE_SUF_LENGTH - 1 < seq.size(); ++i) {
				if(seq[i] == 'A' && seq[i + 1] == 'T' && seq[i + 2] == 'G') {
					auto d = edit_distance_int(seq.data() + i, PRE_SUF_LENGTH, ref_prefix.data(), PRE_SUF_LENGTH, prefix_distance);
					if(d && (int) d->size() < prefix_distance) {
						prefix = i;
						prefix_distance = d->size();
					}
				}
			}
		}
		{
			for(size_t i = PRE_SUF_LENGTH - 3; i + 2 < seq.size(); ++i) {
				if(seq[i] == 'T' && ((seq[i + 1] == 'A' && seq[i + 2] == 'G') || (seq[i + 1] == 'A' && seq[i + 2] == 'A') || (seq[i + 1] == 'G' && seq[i + 2] == 'A'))) {
					auto d = edit_distance_int(seq.data() + i + 3 - PRE_SUF_LENGTH, PRE_SUF_LENGTH, ref_suffix.data(), PRE_SUF_LENGTH, suffix_distance);
					if(d && (int) d->size() < suffix_distance) {
						suffix = i + 3;
						suffix_distance = d->size();
					}
				}
			}
		}
		if(prefix == std::string::npos || suffix == std::string::npos) {
			continue;
		}
		if(prefix > suffix) {
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
		try {
			genes_compared.push_back(genome->genes.back()->compare_to(ref_gene));
		} catch (const std::runtime_error&) {
			genome->genes.pop_back();
		}
	}
	std::unique_ptr<Genome::ComparisonResult> cr = std::make_unique<Genome::ComparisonResult>(*genome, reference, std::move(genes_compared), genome->distance_to(reference, genes_compared));
	return std::make_pair(std::move(genome), std::move(cr));
}

std::vector<std::pair<std::string, std::string>> Genome::load_genomes_raw(std::istream& input, std::ostream& info) {
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

std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> Genome::calculate_genomes(std::istream& input, const Genome& reference, std::ostream& info) {
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

	size_t distance_errors = 0, incomplete_genomes = 0;
	for(const auto&[genome, cr] : genomes) {
		if(cr->distance < 0) {
			++distance_errors;
		}
		if(genome->genes.size() < reference.genes.size()) {
			++incomplete_genomes;
		}
	}
	info << distance_errors << " distance errors" << std::endl;
	info << incomplete_genomes << " incomplete genomes" << std::endl;
	info << "In total we have processed " << genomes.size() << " genomes" << std::endl;
	return genomes;
}

Genome::ComparisonResult Genome::compare_to(const Genome& reference) const {
	std::vector<Gene::ComparisonResult> genes_compared;
	for(const auto& gene: genes) {
		for(const auto& gene_reference: reference.genes) {
			if(gene_reference->name == gene->name) {
				try {
					genes_compared.emplace_back(gene->compare_to(*gene_reference));
				} catch (const std::runtime_error&) {}
				break;
			}
		}
	}
	int distance = distance_to(reference, genes_compared);
	return Genome::ComparisonResult(*this, reference, std::move(genes_compared), distance);
}

std::unique_ptr<Genome> Genome::from_parsed_python(const python_print_type& p, std::string data) {
	auto genome = std::unique_ptr<Genome>(new Genome(std::get<0>(p.as_tuple()), std::move(data), {}));
	for(const Gene::python_print_type& gene : std::get<1>(p.as_tuple())) {
		genome->genes.push_back(Gene::from_parsed_python(gene, *genome));
	}
	return genome;
}

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
	return Gene::ComparisonResult(*compared, *reference, std::get<1>(p.as_tuple()), std::get<2>(p.as_tuple()));
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
		*this,
		reference,
		*edit_distance_int(reference.genome.data.data() + reference.begin, reference.end - reference.begin, genome.data.data() + begin, end - begin),
		*edit_distance_int(reference.protein.data(), reference.protein.size(), protein.data(), protein.size())
	);
}

std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> Genome::load_from_files(const Genome& reference, std::istream& genomes_file, std::istream& input_file, std::ostream& info) {
	std::vector<std::pair<std::string, std::string>> genomes_raw = load_genomes_raw(genomes_file, info); //header, data
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> output;
	info << "Parsing distances...";
	info << std::endl;
	{
		using DataT = std::pair<Genome::python_print_type, Genome::ComparisonResult::python_print_type>;
		size_t size;
		input_file >> size;
		if(input_file.bad() || input_file.fail()) {
			throw std::runtime_error("Error reading from input_file in Genome::load_from_files");
		}
		std::vector<DataT> genomes_extra_data;
		genomes_extra_data.reserve(size);
		output.reserve(size);
		for(size_t i = 0; i < size; ++i) {
			genomes_extra_data.push_back(python_loader<DataT>::load(input_file));
		}
		if(size != genomes_raw.size()) {
			throw std::runtime_error("Size mismatch, recalculate!");
		}
		for(size_t i = 0; i < size; ++i) {
			auto genome = Genome::from_parsed_python(genomes_extra_data[i].first, std::move(genomes_raw[i].second));
			if(genome->header != genomes_raw[i].first) {
				throw std::runtime_error("Header mismatch at position " + std::to_string(i));
			}
			auto comparison = std::make_unique<Genome::ComparisonResult>(Genome::ComparisonResult::from_parsed_python(genomes_extra_data[i].second, *genome, reference));
			output.emplace_back(std::move(genome), std::move(comparison));
		}
	}
	info << "Parsing done.";
	info << std::endl;
	return output;
}

void print_as_python(std::ostream& o, const Genome& genome) {
	print_as_python_dict(o, "header"s, genome.header, "genes"s, genome.genes);
}

void print_as_python(std::ostream& o, const Genome::ComparisonResult& cr) {
	print_as_python_dict(o, "compared"s, cr.compared.header, "reference"s, cr.reference.header, "genes"s, cr.genes, "distance"s, cr.distance);
}
