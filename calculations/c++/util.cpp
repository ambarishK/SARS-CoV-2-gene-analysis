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

template<typename T>
TempFile<T>::TempFile() : T(nullptr), filename("/tmp/mytemp.XXXXXX"), filebuf() {
	int fd = mkstemp(filename.data());
	if(fd == 0) {
		throw std::runtime_error("Failed to create a temporary file.");
	}
	std::ios_base::openmode openmode;
	if constexpr (std::is_same_v<T, std::istream>) {
		openmode = std::ios_base::in;
	} else if constexpr (std::is_same_v<T, std::ostream>) {
		openmode = std::ios_base::out;
	} else if constexpr (std::is_same_v<T, std::iostream>) {
		openmode = std::ios_base::in | std::ios_base::out;
	}
	filebuf = __gnu_cxx::stdio_filebuf<char>(fd, openmode);
	this->rdbuf(&filebuf);
}
template<typename T>
const std::string& TempFile<T>::get_filename() const {
	return filename;
}
template<typename T>
TempFile<T>::~TempFile() {
	this->exceptions(std::ios_base::goodbit);
	if constexpr(std::is_same_v<T, std::ostream> || std::is_same_v<T, std::iostream>) {
		this->flush();
	}
	this->rdbuf(nullptr);
}

template class TempFile<std::istream>;
template class TempFile<std::ostream>;
template class TempFile<std::iostream>;
