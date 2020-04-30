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

std::vector<std::tuple<size_t, char, char>> get_differences(const char* str_a, const char* str_b, size_t len) {
	std::vector<std::tuple<size_t, char, char>> differences;
	for(size_t i = 0; i < len; ++i) {
		if(str_a[i] != str_b[i]) {
			differences.emplace_back(i, str_a[i], str_b[i]);
		}
	}
	return differences;
}
