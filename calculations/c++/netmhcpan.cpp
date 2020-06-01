#include "paths.h"
#include "edit_distance.h"
#include "genome.h"
#include "python_printer.h"
#include "util.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>

using namespace filenames;

const std::string NETMHCPAN_LOCATION = data_path("netMHCpan");

std::vector<std::string> get_mhc_names() {
	std::string mhc_pseudo_filename;
	std::string allelenames_filename;
	{
		TempFile<std::iostream> tmp_file;
		auto netmhcpan_return_code = system((NETMHCPAN_LOCATION + " > " + tmp_file.get_filename()).data());
		if(netmhcpan_return_code) {
			throw std::runtime_error("Invalid return code from netMHCpan " + std::to_string(netmhcpan_return_code));
		}
		for(int i = 0; i < 11; ++i) {
			tmp_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		std::getline(tmp_file, mhc_pseudo_filename);
		for(int i = 0; i < 17; ++i) {
			tmp_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		std::getline(tmp_file, allelenames_filename);
		mhc_pseudo_filename = mhc_pseudo_filename.substr(24, mhc_pseudo_filename.size() - 24 - 31);
		allelenames_filename = allelenames_filename.substr(23, allelenames_filename.size() - 23 - 34);
	}
	std::set<std::string> mhc_pseudo_entries;
	{
		std::ifstream mhc_pseudo = open_file_i(mhc_pseudo_filename);
		std::string entry;
		while(true) {
			mhc_pseudo >> entry;
			if(mhc_pseudo) {
				mhc_pseudo_entries.insert(std::move(entry));
				mhc_pseudo.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			} else {
				break;
			}
		}
	}
	std::vector<std::string> result;
	std::ifstream allelenames = open_file_i(allelenames_filename);
	std::string entry;
	while(true) {
		allelenames >> entry;
		if(allelenames) {
			if(mhc_pseudo_entries.count(entry)) {
				result.push_back(std::move(entry));
			}
			allelenames.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} else {
			return result;
		}
	}
}

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
	double rank_el;

	static std::vector<NetMHCPanResultLine> parse(std::istream& input) {
		for(int i = 0; i < 49; ++i) {
			input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		std::string garbage;
		std::vector<NetMHCPanResultLine> result;
		while(true) {
			NetMHCPanResultLine l;
			input >> l.pos >> garbage >> garbage >> l.core >> l.of >> l.gp >> l.gl >> l.ip >> l.il >> l.icore >> garbage >> l.score_el >> l.rank_el;
			input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			if(input) {
				result.push_back(std::move(l));
			} else {
				return result;;
			}
		}
	}
private:
	NetMHCPanResultLine() = default;
	NetMHCPanResultLine(int pos, std::string core, int of, int gp, int gl, int ip, int il, std::string icore, double score_el, double rank_el) : pos(pos), core(std::move(core)), of(of), gp(gp), gl(gl), ip(ip), il(il), icore(std::move(icore)), score_el(score_el), rank_el(rank_el) {}
	friend struct python_loader<NetMHCPanResultLine>;
};

void print_as_python(std::ostream& o, const NetMHCPanResultLine& l) {
	print_as_python_dict(o, "pos", l.pos, "core", l.core, "of", l.of, "gp", l.gp, "gl", l.gl, "ip", l.ip, "il", l.il, "icore", l.icore, "score_el", l.score_el, "rank_el", l.rank_el);
}

template<>
struct python_loader<NetMHCPanResultLine> {
	static NetMHCPanResultLine load(std::istream& input) {
		return std::apply([](auto&&... args){return NetMHCPanResultLine(std::forward<decltype(args)>(args)...);}, decltype(python_loader_struct<int, std::string, int, int, int, int, int, std::string, double, double>::with_names_values("pos"_sl, "core"_sl, "of"_sl, "gp"_sl, "gl"_sl, "ip"_sl, "il"_sl, "icore"_sl, "score_el"_sl, "rank_el"_sl))::load(input).as_tuple());
	}
};

std::tuple<size_t, size_t, std::map<std::string, std::vector<std::pair<size_t, std::string>>>> get_data_to_netmhcpan(std::istream& genomes_file, std::istream& genomes_raw_file, ssize_t max_genomes, const std::vector<std::tuple<std::string, std::string, size_t>>& peptides, std::ostream& info) {
	std::map<std::string, std::vector<std::pair<size_t, std::string>>> to_netmhcpan; // peptide_in_genome -> {{id, reference_peptide}}
	size_t size, good = 0, skipped = 0;
	genomes_file >> size;
	std::vector<std::pair<std::string, std::string>> genomes_raw = Genome::load_genomes_raw(genomes_raw_file, info);
	if(size != genomes_raw.size()) {
		throw std::runtime_error("Size mismatch, recalculate!");
	}
	if(max_genomes > 0) {
		info << "Genomes limit is set at " << max_genomes << std::endl;
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
			to_netmhcpan[genome_protein->substr(peptide_in_genome_location, reference_peptide->size())].emplace_back(i, *reference_peptide);
		}
	}
	info << good << " good entries, " << skipped << " skipped, " << to_netmhcpan.size() << " different peptides to analyze." << std::endl;
	return {good, skipped, to_netmhcpan};
}

void print_groups(std::ostream& out, const std::map<std::string, std::vector<std::pair<size_t, std::string>>>& to_netmhcpan) {
	out << to_netmhcpan.size() << "\n";
	for(const auto&[peptide_in_genome, v] : to_netmhcpan) {
		std::unordered_map<std::string, std::vector<size_t>> different_peptides;
		for(const auto&[id, reference_peptide] : v) {
			different_peptides[reference_peptide].push_back(id);
		}
		out << peptide_in_genome << '\n';
		print_as_python(out, different_peptides);
		out << '\n';
	}
}

int main(int argc, char* argv[]){
	if(argc < 3) {
		std::cerr << "USAGE:\nnetmhcpan MAX_GENOMES MAX_MHCs\n";
		return 1;
	}
	long long max_genomes = std::stoll(argv[1]);
	long long max_mhcs = std::stoll(argv[2]);
	std::vector<std::string> mhcs = get_mhc_names();
	std::cerr << "Total MHCs: " << mhcs.size() << std::endl;
	if(max_mhcs > 0 && (size_t) max_mhcs < mhcs.size()) {
		std::cerr << "Choosing " << max_mhcs << std::endl;
		std::vector<std::string> mhcs_sample;
		mhcs_sample.reserve(max_mhcs);
		std::sample(mhcs.begin(), mhcs.end(), std::back_inserter(mhcs_sample), max_mhcs, std::mt19937_64{/*std::random_device{}()*/});
		mhcs = std::move(mhcs_sample);
	}

	std::map<std::string, std::vector<std::pair<size_t, std::string>>> to_netmhcpan; // peptide_in_genome -> {{id, reference_peptide}}
	{
		std::ifstream genomes_file = open_file_i(data_path(CALCULATED_GENOMES_DATA));
		std::ifstream genomes_raw_file = open_file_i(data_path(FILTERED_GENOMES));
		std::ifstream peptides_file = open_file_i(data_path(PEPTIDES_INPUT));
		std::vector<std::tuple<std::string, std::string, size_t>> peptides = load_peptides(peptides_file);
		auto t = get_data_to_netmhcpan(genomes_file, genomes_raw_file, max_genomes, peptides, std::cerr);
		to_netmhcpan = std::move(std::get<2>(t));
	}
	if(to_netmhcpan.empty()) {
		throw std::runtime_error("Empty input or all skipped!");
	}

	std::cerr << "Printing netMHCpan input..." << std::endl;
	TempFile<std::ostream> netmhc_input;
	for(const auto&[peptide_in_genome, _] : to_netmhcpan) {
		netmhc_input << peptide_in_genome << '\n';
	}
	netmhc_input.flush();
	std::cerr << "netMHCpan input printed." << std::endl;

	{
		std::cerr << "Printing groups..." << std::endl;
		print_groups(std::cout, to_netmhcpan);
		std::cerr << "Groups printed." << std::endl;
	}
	std::cout << argc - 2  << '\n';
	std::mutex output_mutex;
	for_each_with_progress(mhcs.begin(), mhcs.end(), [&output_mutex, expected_size{to_netmhcpan.size()}, &netmhc_input_filename = netmhc_input.get_filename()](const std::string& mhc){
		TempFile<std::istream> netmhc_output;
		{
			std::string command = NETMHCPAN_LOCATION + " -a " + mhc + " -p -f " + netmhc_input_filename + " > " + netmhc_output.get_filename();
			auto netmhcpan_return_code = system(command.data());
			if(netmhcpan_return_code) {
				throw std::runtime_error("Invalid return code from netMHCpan " + std::to_string(netmhcpan_return_code));
			}
		}
		std::vector<NetMHCPanResultLine> parsed = NetMHCPanResultLine::parse(netmhc_output);
		if(parsed.size() != expected_size) {
			throw std::runtime_error("Invalid number of records from netMHCpan " + std::to_string(parsed.size()) + " instead of " + std::to_string(expected_size));
		}
		{
			const std::lock_guard<std::mutex> lock(output_mutex);
			std::cout << mhc << '\n';
			for(const NetMHCPanResultLine& result_line: parsed) {
				print_as_python(std::cout, result_line);
				std::cout << '\n';
			}
		}
	}, std::cerr);
}
