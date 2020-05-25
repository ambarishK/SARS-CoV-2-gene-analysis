#ifndef DEBUG
#define DEBUG 0
#endif

#include "genome.h"
#include "paths.h"
#include "util.h"
#include "edit_distance.h"

#include <map>
#include <set>

using namespace std::string_literals;
using namespace filenames;

std::optional<EditOperation> similar_genomes(const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& a, const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>& b) {
	if(a.second->distance < 0 || b.second->distance < 0) {
		return std::nullopt;
	}
	if((size_t) std::abs(a.second->distance - b.second->distance) > 1) {
		return std::nullopt;
	}
	const char* a_str = a.first->data.data();
	size_t a_len = a.first->data.size();
	const char* b_str = b.first->data.data();
	size_t b_len = b.first->data.size();
	while(a_len && b_len && a_str[0] == b_str[0]) {
		++a_str;
		++b_str;
		--a_len;
		--b_len;
	}
	while(a_len && b_len && a_str[a_len - 1] == b_str[b_len - 1]) {
		--a_len;
		--b_len;
	}
	if(a_len == 1 && b_len == 0) {
		return {EditOperation(a_str - a.first->data.data(), a_str[0], EditOperation::Type::DELETE)};
	}
	if(a_len == 0 && b_len == 1) {
		return {EditOperation(a_str - a.first->data.data(), b_str[0], EditOperation::Type::INSERT)};
	}
	if(a_len == 1 && b_len == 1) {
		return {EditOperation(a_str - a.first->data.data(), b_str[0], EditOperation::Type::SUBSTITUTE)};
	}
	return {};
}

void find_similar_genomes(const std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>>& genomes, const std::string& groups_output_filename, const std::string& comparisons_output_filename) {
	std::cout << "Grouping genomes by distance to reference..." << std::endl;
	std::map<decltype(Genome::ComparisonResult::distance), std::vector<const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>*>> by_edit_distance_to_reference;

	for(const auto& g : genomes) {
		auto distance = g.second->distance;
		if(distance >= 0) {
			by_edit_distance_to_reference[distance].push_back(&g);
		}
	}
	std::cout << "Finished grouping genomes, " << by_edit_distance_to_reference.size() << " distances." << std::endl;

	std::map<decltype(Genome::ComparisonResult::distance), std::vector<std::vector<const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>*>>> by_edit_distance_to_reference_equal;

	size_t group_count = 0;
	std::cout << "Grouping same genomes..." << std::endl;
	for(const auto&[distance, gens] : by_edit_distance_to_reference) {
		std::vector<std::vector<const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>*>> partition;
		for(const auto genome : gens) {
			bool match = false;
			for(auto& p : partition) {
				if(genome->first->data == (**p.begin()).first->data) {
					p.push_back(genome);
					match = true;
					break;
				}
			}
			if(!match) {
				partition.push_back({genome});
				++group_count;
			}
		}
		by_edit_distance_to_reference_equal.emplace(distance, move(partition));
		by_edit_distance_to_reference.at(distance).clear();
	}
	by_edit_distance_to_reference.clear();
	std::cout << "Finished grouping same genomes, " << group_count << " groups." << std::endl;

	std::cout << "Comparing groups..." << std::endl;
	std::vector<std::tuple<std::pair<decltype(Genome::ComparisonResult::distance), decltype(Genome::ComparisonResult::distance)>, std::pair<size_t, size_t>, EditOperation>> results;

	{
		std::mutex results_mutex;
		std::vector<decltype(Genome::ComparisonResult::distance)> distances;
		distances.reserve(by_edit_distance_to_reference_equal.size());
		for(const auto&[distance, gens] : by_edit_distance_to_reference_equal) {
			distances.push_back(distance);
		}
		for_each_with_progress(distances.begin(), distances.end(), [&by_edit_distance_to_reference_equal, &results_mutex, &results](decltype(Genome::ComparisonResult::distance) distance){
			auto other_gens = by_edit_distance_to_reference_equal.find(distance + 1);
			if(other_gens == by_edit_distance_to_reference_equal.end()) {
				return;
			}
			const auto& gens = by_edit_distance_to_reference_equal.at(distance);
			for(size_t a = 0; a < gens.size(); ++a) {
				const std::vector<const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>*>& a_gens = gens[a];
				for(size_t b = 0; b < other_gens->second.size(); ++b) {
					const std::vector<const std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>*>& b_gens = other_gens->second[b];
					auto r = similar_genomes(**a_gens.begin(), **b_gens.begin());
					if(r) {
						std::lock_guard<std::mutex> lock(results_mutex);
						results.emplace_back(std::make_pair(distance, distance + 1), std::make_pair(a, b), std::move((*r)));
					}
				}
			}
		}, std::cout);
	}
	std::cout << "Finished comparing groups, there are " << results.size() << " matches." << std::endl;


	{
		std::cout << "Printing groups..." << std::endl;
		std::ofstream output = open_file_o(groups_output_filename);
		for(const auto&[distance, gens] : by_edit_distance_to_reference_equal) {
			std::vector<std::vector<std::string>> groups;
			groups.reserve(gens.size());
			for(const auto& eq : gens) {
				std::vector<std::string> names;
				names.reserve(eq.size());
				for(const auto& genome : eq) {
					names.push_back(genome->first->header);
				}
				groups.push_back(move(names));
			}
			print_as_python_dict(output, "distance", distance, "groups", groups);
			output << '\n';
		}
		std::cout << "Finished printing groups." << std::endl;
	}

	{
		std::cout << "Printing comparisons..." << std::endl;
		std::ofstream output = open_file_o(comparisons_output_filename);
		for(const auto& [distances, group_ids, edit_operation] : results) {
			print_as_python_dict(output,
				"distance_compared", std::get<1>(distances),
				"distance_reference", std::get<0>(distances),
				"group_id_compared", std::get<1>(group_ids),
				"group_id_reference", std::get<0>(group_ids),
				"mutation", edit_operation);
			output << '\n';
		}
		std::cout << "Finished printing comparisons." << std::endl;
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
		find_similar_genomes(genomes, data_path("relatives_groups.txt"), data_path("relatives_comparisons.txt"));
		break;
	case 4:
		find_similar_genomes_by_sequence(genomes, argv[1], std::stoll(argv[2]), std::stoll(argv[3]), data_path("relatives_seq_sequences.txt"), data_path("relatives_seq_comparisons.txt"));
		break;
	default:
		std::cout << "Wrong usage.\n";
		return -1;
	}
}
