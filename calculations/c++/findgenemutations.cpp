#ifndef DEBUG
#define DEBUG 0
#endif

#include <iostream>
#include <vector>
#include <string>
#include <string_view>
#include <exception>
#include <unordered_map>
#include <fstream>
#include <cctype>

#include "edit_distance.h"
#include "python_printer.h"
#include "string_tools.h"
#include "paths.h"
#include "genome.h"
#include "util.h"

using namespace std::string_literals;

void calculate_genomes(const Genome& reference, const std::string& genomes_filename, const std::string& output_filename) {
	std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> genomes;
	{
		auto input = open_file_i(genomes_filename);
		genomes = Genome::calculate_genomes(input, reference, std::cout);
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
	std::ifstream genomes_file = open_file_i(genomes_filename);
	std::ifstream input_file = open_file_i(input_filename);
	return Genome::load_from_files(reference, genomes_file, input_file, std::cout);
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
