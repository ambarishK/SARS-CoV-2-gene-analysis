#include "paths.h"

std::string filenames::data_path(const char* file, int sub_level) {
	std::string result;
	for(int i = 0; i < sub_level; ++i) {
		result += "../";
	}
	return result + "../../data/" + file;
}
