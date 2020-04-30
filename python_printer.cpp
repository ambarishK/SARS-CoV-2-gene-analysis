#include "python_printer.h"

#include <iomanip>

void print_as_python(std::ostream& o, const std::string& string) {
	o << std::quoted(string, '\'');
}

void print_as_python(std::ostream& o, char c) {
	print_as_python(o, std::string(1, c));
}

void print_as_python_tuple_(std::ostream&) {}

void print_as_python_dict_(std::ostream&) {}

std::string python_loader<std::string>::load(std::istream& input) {
	std::string str;
	input >> std::quoted(str, '\'');
	if(input.bad() || input.fail()) {
		throw std::runtime_error("Incorrect input.");
	}
	return str;
}

char python_loader<char>::load(std::istream& input) {
	std::string str;
	input >> std::quoted(str, '\'');
	if(input.bad() || input.fail() || str.size() != 1) {
		throw std::runtime_error("Incorrect input.");
	}
	return str[0];
}
