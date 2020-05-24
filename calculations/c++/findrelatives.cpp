#include "genome.h"
#include "paths.h"
#include "util.h"
#include "edit_distance.h"

using namespace std::string_literals;
using namespace filenames;

#ifndef DEBUG
#define DEBUG 0
#endif

bool similar_genomes(const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& a, const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& b, size_t max_distance = 1) {
	if(a.second->distance < 0 || b.second->distance < 0) {
		return false;
	}
	if((size_t) std::abs(a.second->distance - b.second->distance) > max_distance) {
		return false;
	}
	try {
		return (bool) edit_distance_int(a.first->data.data(), a.first->data.size(), b.first->data.data(), b.first->data.size(), max_distance);
	} catch (...) {
		return false;
	}
}

void find_similar_genomes(const std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>>& genomes, const std::string& output_filename) {
	std::vector<std::vector<const Genome*>> results(genomes.size());
	transform_with_progress(genomes.begin(), genomes.end(), results.begin(), [&genomes](const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& genome) {
		size_t index = &genome - genomes.data();
		std::vector<const Genome*> result;
		for(size_t i = index + 1; i < genomes.size(); ++i) {
			if(similar_genomes(genome, genomes[i])) {
				result.push_back(genomes[i].first.get());
			}
		}
		return result;
	}, std::cout);
	size_t total_results = 0;
	for(const auto& v : results) {
		total_results += v.size();
	}
	std::ofstream output = open_file_o(output_filename);
	output << total_results << '\n';
	for(size_t i = 0; i < genomes.size(); ++i) {
		for(const Genome* o: results[i]) {
			output << genomes[i].first->header << '\n' << o->header << '\n';
		}
	}
}

std::string_view match_sequence(const std::string& string, const std::string& seq, size_t width_difference, size_t begin = 0, size_t end = std::numeric_limits<size_t>::max()) {
	std::string_view record;
	size_t record_score = std::numeric_limits<size_t>::max();
	for(size_t size = seq.size() < width_difference ? 0 : seq.size() - width_difference; size < seq.size() + width_difference + 1; ++size) {
		for(size_t i = begin; i < std::min(string.size() - size, end); ++i) {
			auto r = edit_distance_int(string.data() + i, size, seq.data(), seq.size(), record_score);
			if(r && r->size() < record_score) {
				record = std::string_view(string.data() + i, size);
				record_score = r->size();
			}
		}
	}
	return record;
}

void find_similar_genomes_by_sequence(const std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>>& genomes, const std::string& seq, size_t begin, size_t end, const std::string& sequences_output_filename, const std::string& comparisons_output_filename, size_t expected_error = 1) {
	std::vector<std::string_view> sequences(genomes.size());
	transform_with_progress(genomes.begin(), genomes.end(), sequences.begin(), [&seq, expected_error, begin, end](const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& genome){
		return match_sequence(genome.first->data, seq, expected_error, begin, end);
	}, std::cout);
	{
		std::ofstream output = open_file_o(sequences_output_filename);
		for(size_t i = 0; i < genomes.size(); ++i) {
			print_as_python_dict(output, "header", genomes[i].first->header, "begin", sequences[i].begin() - genomes[i].first->data.data(), "end", sequences[i].end() - genomes[i].first->data.data());
			output << "\n";
		}
	}
	std::vector<std::vector<std::pair<size_t, std::vector<EditOperation>>>> results(genomes.size());
	transform_with_progress(sequences.begin(), sequences.end(), results.begin(), [&sequences, expected_error](const std::string_view& sequence) {
		size_t index = &sequence - sequences.data();
		std::vector<std::pair<size_t, std::vector<EditOperation>>> result;
		for(size_t i = index + 1; i < sequences.size(); ++i) {
			auto r = edit_distance_int(sequence.data(), sequence.size(), sequences[i].data(), sequences[i].size(), expected_error);
			if(r && r->size() == expected_error) {
				result.emplace_back(i, std::move(*r));
			}
		}
		return result;
	}, std::cout);
	std::ofstream output = open_file_o(comparisons_output_filename);
	for(size_t i = 0; i < genomes.size(); ++i) {
		for(const auto& [o, mutations]: results[i]) {
			print_as_python_dict(output, "compared", genomes[i].first->header, "reference", genomes[o].first->header, "mutations", mutations);
			output << "\n";
		}
	}
}

std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> load_genomes_from_files(const Genome& reference, const std::string& genomes_filename, const std::string& input_filename) {
	std::ifstream genomes_file = open_file_i(genomes_filename);
	std::ifstream input_file = open_file_i(input_filename);
	return Genome::load_from_files(reference, genomes_file, input_file, std::cout);
}

int main(int argc, char* argv[]){
	std::unique_ptr<Genome> reference;
	{
		auto ref_info_file = open_file_i(data_path(REFERENCE_GENES));
		auto ref_genome_file = open_file_i(data_path(REFERENCE_GENOME));
		reference = Genome::load_reference(ref_info_file, ref_genome_file);
	}
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes = load_genomes_from_files(*reference, data_path(FILTERED_GENOMES), data_path(CALCULATED_GENOMES_DATA));
	std::cout << "Genomes loaded" << std::endl;
	switch(argc) {
	case 1:
		find_similar_genomes(genomes, data_path(EXTRA_COMPARISONS_REQUESTS));
		break;
	case 4:
		find_similar_genomes_by_sequence(genomes, argv[1], std::stoll(argv[2]), std::stoll(argv[3]), data_path("relatives_sequences.txt"), data_path("relatives_comparisons.txt"));
		break;
	default:
		std::cout << "Wrong usage.\n";
		return -1;
	}
}
