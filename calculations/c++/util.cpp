#include "util.h"

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
