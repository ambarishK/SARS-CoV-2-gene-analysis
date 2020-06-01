#pragma once

#include "python_printer.h"
#include "edit_distance.h"

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
		using python_print_type = decltype(python_loader_struct<std::string, std::vector<EditOperation>, std::vector<EditOperation>>::with_names_values("name"_sl, "mutations"_sl, "protein_mutations"_sl));
		const Gene& compared;
		const Gene& reference;
		std::vector<EditOperation> mutations;
		std::vector<EditOperation> protein_mutations;
		ComparisonResult() = delete;
		ComparisonResult(const Gene& compared, const Gene& reference, std::vector<EditOperation> mutations, std::vector<EditOperation> protein_mutations);

		static ComparisonResult from_parsed_python(const python_print_type& p, const Genome& r_compared, const Genome& r_reference);
	};

	ComparisonResult compare_to(const Gene& reference) const;

	static std::unique_ptr<Gene> from_parsed_python(const python_print_type& p, const Genome& genome);

private:
	friend struct Genome;
	Gene() = delete;
	Gene(const Gene&) = delete;
	Gene(Gene&&) = delete;
	Gene(std::string name, const Genome& genome, size_t begin, size_t end, std::string protein);
};

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
		ComparisonResult(const Genome& compared, const Genome& reference, std::vector<Gene::ComparisonResult> genes, int distance);
		static ComparisonResult from_parsed_python(const python_print_type& p, const Genome& compared, const Genome& reference);
	};

	static std::unique_ptr<Genome> from_parsed_python(const python_print_type& p, std::string data);

	int distance_to(const Genome& reference, const std::vector<Gene::ComparisonResult>& genes_compared) const;

	static std::unique_ptr<Genome> load_reference(std::istream& ref_genes_file, std::istream& ref_genome_file);

	static std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>> make_genome(std::string header, std::string data, const Genome& reference, size_t tolerance = 128);

	static std::vector<std::pair<std::string, std::string>> load_genomes_raw(std::istream& input, std::ostream& info);

	static std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> calculate_genomes(std::istream& input, const Genome& reference, std::ostream& info); /* this function has to preserve the order */

	static std::vector<std::pair<std::unique_ptr<Genome>, std::unique_ptr<Genome::ComparisonResult>>> load_from_files(const Genome& reference, std::istream& genomes_file, std::istream& input_file, std::ostream& info);

	Genome::ComparisonResult compare_to(const Genome& reference) const;

private:
	Genome() = delete;
	Genome(const Genome&) = delete;
	Genome(Genome&&) = delete;
	Genome(std::string header, std::string data, std::vector<std::unique_ptr<Gene>> genes) : header(std::move(header)), data(std::move(data)), genes(std::move(genes)) {}
};

void print_as_python(std::ostream& o, const Gene& gene);

void print_as_python(std::ostream& o, const Gene::ComparisonResult& cr);

void print_as_python(std::ostream& o, const Genome& genome);

void print_as_python(std::ostream& o, const Genome::ComparisonResult& cr);

std::string gene_to_protein(std::string_view gene);
