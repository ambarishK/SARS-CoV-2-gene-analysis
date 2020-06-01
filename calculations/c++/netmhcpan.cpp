#include "paths.h"
#include "edit_distance.h"
#include "genome.h"
#include "python_printer.h"
#include "util.h"

#include <cstdlib>
#include <iostream>
#include <string_view>
#include <sstream>
#include <map>

using namespace filenames;

const std::string NETMHCPAN_LOCATION = data_path("netMHCpan");

class PeptideNotFoundException : public std::exception {};

size_t find_peptide_in_genome(const std::string& gene_protein, const std::string& peptide_data, size_t peptide_location, size_t search_radius = 1000, size_t min_peptide_len = 8) {
	size_t record_distance, record_pos = peptide_location;
	if(peptide_data.size() < min_peptide_len) {
		throw std::runtime_error("Too short peptide.");
	}
	if(peptide_location + search_radius < peptide_location || peptide_location + peptide_data.size() < peptide_location || peptide_location + peptide_data.size() + search_radius < peptide_location + peptide_data.size()) {
		throw std::runtime_error("Overflow.");
	}
	if(peptide_location + peptide_data.size() < gene_protein.size()) {
		record_distance = edit_distance_int(gene_protein.data() + peptide_location, peptide_data.size(), peptide_data.data(), peptide_data.size())->size();
	} else {
		record_distance = std::numeric_limits<size_t>::max();
	}
	if(record_distance == 0) {
		return peptide_location;
	}
	for(size_t i = 1; i < search_radius; ++i) {
		size_t pos = peptide_location + i;
		if(pos + peptide_data.size() < gene_protein.size()) {
			auto distance = edit_distance_int(gene_protein.data() + pos, peptide_data.size(), peptide_data.data(), peptide_data.size(), record_distance - 1);
			if(distance) {
				record_pos = pos;
				record_distance = distance->size();
				if(record_distance == 0) {
					return pos;
				}
			}
		}
		pos = peptide_location - i;
		if(pos + peptide_data.size() < gene_protein.size() && pos < peptide_location) {
			auto distance = edit_distance_int(gene_protein.data() + pos, peptide_data.size(), peptide_data.data(), peptide_data.size(), record_distance - 1);
			if(distance) {
				record_pos = pos;
				record_distance = distance->size();
				if(record_distance == 0) {
					return pos;
				}
			}
		}
	}
	if(record_distance == std::numeric_limits<size_t>::max()) {
		throw PeptideNotFoundException();
	}
	return record_pos;
}

std::vector<std::tuple<std::string, std::string, size_t>> load_peptides(std::istream& input) {
	std::string peptide, gene;
	size_t position;
	std::vector<std::tuple<std::string, std::string, size_t>> ret;
	while(input >> peptide >> gene >> position) {
		ret.emplace_back(std::move(peptide), std::move(gene), position);
	}
	return ret;
}

std::vector<std::tuple<const std::string*, std::shared_ptr<std::string>, size_t>> find_peptides_in_genome(const std::vector<std::tuple<std::string, std::string, size_t>>& peptides, const Genome& genome) {
	std::vector<std::tuple<const std::string*, std::shared_ptr<std::string>, size_t>> result; // {reference_peptide, gene_protein, peptide_in_genome_location}
	for(const std::unique_ptr<Gene>& gene : genome.genes) {
		if(gene->invalid_nucleotides() > 0) {
			continue;
		}
		std::shared_ptr<std::string> protein;
		for(const auto& [peptide_data, peptide_gene, peptide_location] : peptides) {
			if(peptide_gene == gene->name) {
				if(!protein) {
					protein = std::make_unique<std::string>(gene_to_protein(std::string_view(genome.data.data() + gene->begin, gene->end - gene->begin)));
				}
				try {
					result.emplace_back(&peptide_data, protein, find_peptide_in_genome(*protein, peptide_data, peptide_location));
				} catch (const PeptideNotFoundException&) {}
			}
		}
	}
	return result;
}

struct NetMHCPanResultLine {
	int pos;
	std::string core;
	int of;
	int gp;
	int gl;
	int ip;
	int il;
	std::string icore;
	double score_el;
	double bindlevel;

	static constexpr const char* csv_header = "pos,core,of,gp,gl,ip,il,icore,score_el,bindlevel";

	static std::vector<NetMHCPanResultLine> parse(std::istream& input) {
		for(int i = 0; i < 49; ++i) {
			input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		std::string garbage;
		std::vector<NetMHCPanResultLine> result;
		while(true) {
			NetMHCPanResultLine l;
			input >> l.pos >> garbage >> garbage >> l.core >> l.of >> l.gp >> l.gl >> l.ip >> l.il >> l.icore >> garbage >> l.score_el >> l.bindlevel;
			if(input) {
				result.push_back(std::move(l));
			} else {
				break;
			}
		}
		return result;
	}
private:
	NetMHCPanResultLine() = default;
};

std::ostream& operator<<(std::ostream& o, const NetMHCPanResultLine& n) {
	return o << n.pos << ',' << n.core << ',' << n.of << ',' << n.gp << ',' << n.gl << ',' << n.ip << ',' << n.il << ',' << n.icore << ',' << n.score_el << ',' << n.bindlevel;
}

int main(int argc, char* argv[]){
	if(argc < 3) {
		std::cerr << "USAGE:\nnetmhcpan MAX_GENOMES HLAs...\n";
		return 1;
	}
	long long max_genomes = std::stoll(argv[1]);
	std::unique_ptr<Genome> reference;
	{
		std::ifstream genes = open_file_i(data_path(REFERENCE_GENES));
		std::ifstream genome = open_file_i(data_path(REFERENCE_GENOME));
		reference = Genome::load_reference(genes, genome);
	}
	std::vector<std::tuple<std::string, std::string, size_t>> peptides;
	{
		std::ifstream peptides_file = open_file_i(data_path(PEPTIDES_INPUT));
		peptides = load_peptides(peptides_file);
	}
	size_t good = 0, skipped = 0;
	std::map<std::string, std::vector<std::pair<std::string, std::string>>> to_netmhcpan; // peptide_in_genome -> {{header, reference_peptide}}
	{
		std::ifstream genomes_file = open_file_i(data_path(CALCULATED_GENOMES_DATA));
		std::ifstream genomes_raw_file = open_file_i(data_path(FILTERED_GENOMES));
		size_t size;
		genomes_file >> size;
		std::vector<std::pair<std::string, std::string>> genomes_raw = Genome::load_genomes_raw(genomes_raw_file, std::cerr);
		if(size != genomes_raw.size()) {
			throw std::runtime_error("Size mismatch, recalculate!");
		}
		if(max_genomes > 0) {
			std::cerr << "Genomes limit is set at " << max_genomes << std::endl;
			size = std::min(size, (size_t) max_genomes);
		}
		for(size_t i = 0; i < size; ++i) {
			genomes_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			genomes_file.get();
			std::unique_ptr<Genome> genome = Genome::from_parsed_python(python_loader<Genome::python_print_type>::load(genomes_file), std::move(genomes_raw[i].second));
			if(genome->header != genomes_raw[i].first) {
				throw std::runtime_error("Header mismatch at position " + std::to_string(i));
			}
			std::vector<std::tuple<const std::string*, std::shared_ptr<std::string>, size_t>> peptides_in_genome = find_peptides_in_genome(peptides, *genome);
			good += peptides_in_genome.size();
			skipped += peptides.size() - peptides_in_genome.size();
			for(const auto&[reference_peptide, genome_protein, peptide_in_genome_location] : peptides_in_genome) {
				to_netmhcpan[genome_protein->substr(peptide_in_genome_location, reference_peptide->size())].emplace_back(genomes_raw[i].first, *reference_peptide);
			}
		}
	}
	if(to_netmhcpan.empty()) {
		throw std::runtime_error("Empty input or all skipped!");
	}
	std::cerr << good << " good entries, " << skipped << " skipped, " << to_netmhcpan.size() << " different peptides to analyze." << std::endl;
	TempFile<std::ostream> netmhc_input;
	for(const auto&[peptide_in_genome, _] : to_netmhcpan) {
		netmhc_input << peptide_in_genome << '\n';
	}
	netmhc_input.flush();
	std::cout << "header,mhc,reference_peptide,matched_peptide," << NetMHCPanResultLine::csv_header << '\n';
	for(int i = 2; i < argc; ++i) {
		TempFile<std::istream> netmhc_output;
		const char* mhc = argv[i];
		{
			std::string command = NETMHCPAN_LOCATION + " -a " + mhc + " -p -f " + netmhc_input.get_filename() + " > " + netmhc_output.get_filename();
			auto netmhcpan_return_code = system(command.data());
			if(netmhcpan_return_code) {
				throw std::runtime_error("Invalid return code from netMHCpan " + std::to_string(netmhcpan_return_code));
			}
		}
		std::vector<NetMHCPanResultLine> parsed = NetMHCPanResultLine::parse(netmhc_output);
		if(parsed.size() != to_netmhcpan.size()) {
			throw std::runtime_error("Invalid number of records from netMHCpan " + std::to_string(parsed.size()) + " instead of " + std::to_string(to_netmhcpan.size()));
		}
		auto parsed_result_iter = parsed.begin();
		for(const auto&[peptide_in_genome, v] : to_netmhcpan) {
			std::string common_part = [&](){std::stringstream ss; ss << *parsed_result_iter; return ss.str();}();
			for(const auto&[header, reference_peptide] : v) {
				std::cout << header << ',' << mhc << ',' << reference_peptide << ',' << peptide_in_genome << ',' << common_part << '\n';
			}
			++parsed_result_iter;
		}
	}
}
