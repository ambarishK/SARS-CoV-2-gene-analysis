#include "genome.h"
#include "paths.h"
#include "util.h"
#include "edit_distance.h"

#include <iostream>

using namespace std::string_literals;
using namespace filenames;

#ifndef DEBUG
#define DEBUG 0
#endif

bool similar_genomes(const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& a, const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& b, size_t max_distance = 1) {
	if(std::abs(a.second->distance - b.second->distance) > max_distance) {
		return false;
	}
	return (bool) edit_distance_int(a.first->header.data(), a.first->header.size(), b.first->header.data(), b.first->header.size(), max_distance);
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

std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> load_genomes_from_files(const Genome& reference, const std::string& genomes_filename, const std::string& input_filename) {
	std::ifstream genomes_file = open_file_i(genomes_filename);
	std::ifstream input_file = open_file_i(input_filename);
	return Genome::load_from_files(reference, genomes_file, input_file, std::cout);
}

int main(){
	std::unique_ptr<Genome> reference;
	{
		auto ref_info_file = open_file_i(data_path(REFERENCE_GENES));
		auto ref_genome_file = open_file_i(data_path(REFERENCE_GENOME));
		reference = Genome::load_reference(ref_info_file, ref_genome_file);
	}
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes = load_genomes_from_files(*reference, data_path(FILTERED_GENOMES), data_path(CALCULATED_GENOMES_DATA));
	std::cout << "Genomes loaded" << std::endl;
	find_similar_genomes(genomes, data_path(EXTRA_COMPARISONS_REQUESTS));
}
