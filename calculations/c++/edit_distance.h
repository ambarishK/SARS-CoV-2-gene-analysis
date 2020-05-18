#pragma once

#include <vector>
#include <variant>
#include <limits>
#include <tuple>
#include <optional>
#include <string>

#include "python_printer.h"

struct EditOperation {
	enum class Type : char {
		DELETE = 'd',
		INSERT = 'i',
		SUBSTITUTE = 's'
	};
	uint16_t position;
	char arg;
	Type type;
	EditOperation() = delete;
	EditOperation(uint16_t position, char arg, Type type);
};

void print_as_python(std::ostream& o, const EditOperation::Type& eo);

template<>
struct python_loader<EditOperation::Type> {
	static EditOperation::Type load(std::istream& input);
};

void print_as_python(std::ostream& o, const EditOperation& eo);

template<>
struct python_loader<EditOperation> {
	static EditOperation load(std::istream& input);
};

std::optional<std::vector<EditOperation>> edit_distance_int(const char* source, uint16_t source_len, const char* target, uint16_t target_len, uint16_t limit = std::numeric_limits<uint16_t>::max());

std::string rebuild(const std::vector<EditOperation>& eos, const char *string1, uint16_t len1);
