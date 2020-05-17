#include "string_tools.h"

void rstrip(std::string& string) {
	for(size_t i = string.size() - 1; i != std::string::npos; --i) {
		if(string[i] == ' ' || string[i] == '\t' || string[i] == '\n' || string[i] == '\r') {
			string.resize(i);
		} else {
			break;
		}
	}
}

std::string concatenate_strings(const std::vector<std::string>& strings, size_t reserve_size) {
	std::string result;
	result.reserve(reserve_size);
	for(const std::string& string : strings) {
		result += string;
	}
	return result;
}
