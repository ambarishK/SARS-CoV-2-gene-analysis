#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstring>
#include <algorithm>
#include <vector>
#include <variant>
#include <limits>
#include <tuple>
#include <string>
#include <optional>

namespace {
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

EditOperation::EditOperation(uint16_t position, char arg, Type type) : position(position), arg(arg), type(type) {}

static constexpr const uint16_t MAX_PIECE_LEN = 1000;
static constexpr const uint16_t SKIP_LEN = 50;
static constexpr const uint16_t SKIP_TRIES = 10;

template<typename T>
struct TwiceSizeInt;

template<>
struct TwiceSizeInt<uint16_t> {
	using Type = uint32_t;
};

template<typename T>
struct SimpleIntsPair {
	T first;
	T second;
	SimpleIntsPair() = delete;
	SimpleIntsPair(T first, T second) : first(first), second(second) {}
	bool operator==(const SimpleIntsPair& other) const {
		return reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(*this) == reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(other);
	}
	bool operator!=(const SimpleIntsPair& other) const {
		return reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(*this) != reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(other);
	}
	SimpleIntsPair(const SimpleIntsPair& other) {
		reinterpret_cast<typename TwiceSizeInt<T>::Type&>(*this) = reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(other);
	}
	SimpleIntsPair& operator=(const SimpleIntsPair& other) {
		reinterpret_cast<typename TwiceSizeInt<T>::Type&>(*this) = reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(other);
		return *this;
	}
};
}

static bool ignored_nucleotide_pair(char a, char b) {
	static const std::array<std::array<unsigned char, 8>, 256> IGNORED_PAIRS = [](){
		std::array<std::array<unsigned char, 8>, 256> ignored_pairs{};
		for(size_t i = 0; i < 256; ++i) {
			ignored_pairs[i][0] = i;
		}
		ignored_pairs['A'] = {'V', 'R', 'W', 'M', 'H', 'N', 'A', 'D'};
		ignored_pairs['B'] = {'G', 'T', 'B', 'C', 'U'};
		ignored_pairs['C'] = {'V', 'Y', 'M', 'H', 'N', 'B', 'C', 'S'};
		ignored_pairs['D'] = {'G', 'T', 'D', 'A', 'U'};
		ignored_pairs['G'] = {'V', 'G', 'R', 'K', 'N', 'D', 'B', 'S'};
		ignored_pairs['H'] = {'T', 'H', 'A', 'C', 'U'};
		ignored_pairs['K'] = {'G', 'T', 'K', 'U'};
		ignored_pairs['M'] = {'A', 'C', 'M'};
		ignored_pairs['N'] = {'G', 'T', 'N', 'A', 'C', 'U'};
		ignored_pairs['R'] = {'A', 'R', 'G'};
		ignored_pairs['S'] = {'C', 'G', 'S'};
		ignored_pairs['T'] = {'W', 'T', 'Y', 'H', 'K', 'N', 'D', 'B'};
		ignored_pairs['U'] = {'W', 'Y', 'H', 'K', 'N', 'D', 'B', 'U'};
		ignored_pairs['V'] = {'V', 'G', 'C', 'A'};
		ignored_pairs['W'] = {'A', 'W', 'T', 'U'};
		ignored_pairs['Y'] = {'C', 'U', 'T', 'Y'};
		return ignored_pairs;
	}();
	for(size_t i = 0; i < 8; ++i) {
		if(IGNORED_PAIRS[a][i] == b) {
			return true;
		}
	}
	return false;
}

template<typename SIZE_TYPE, typename CHAR_TYPE>
static std::optional<std::vector<EditOperation>> _edit_distance(const CHAR_TYPE* s1, SIZE_TYPE s1_size, const CHAR_TYPE* s2, SIZE_TYPE s2_size, const SIZE_TYPE limit, bool ignore_uncertain_nucleotides) {
	static constexpr const SIZE_TYPE SAME = 0;
	static constexpr const SIZE_TYPE BEGIN = -1;
	static_assert(static_cast<char>(EditOperation::Type::INSERT) > 0);
	static_assert(static_cast<char>(EditOperation::Type::DELETE) > 0);
	static_assert(static_cast<char>(EditOperation::Type::SUBSTITUTE) > 0);
	std::vector<std::vector<SimpleIntsPair<SIZE_TYPE>>> matrix(s1_size + 1, std::vector<SimpleIntsPair<SIZE_TYPE>>(s2_size + 1, SimpleIntsPair<SIZE_TYPE>(0, 0)));

	matrix[0][0] = SimpleIntsPair<SIZE_TYPE>(0, BEGIN);

	for(SIZE_TYPE i = 1; i <= s2_size; ++i) {
		matrix[0][i] = SimpleIntsPair<SIZE_TYPE>(i, static_cast<SIZE_TYPE>(EditOperation::Type::INSERT));
	}

	for(SIZE_TYPE i = 1; i <= s1_size; ++i) {
		matrix[i][0] = SimpleIntsPair<SIZE_TYPE>(i, static_cast<SIZE_TYPE>(EditOperation::Type::DELETE));
		SIZE_TYPE m = i;
		for(SIZE_TYPE j = 1; j <= s2_size; ++j) {
			SIZE_TYPE ins = matrix[i][j-1].first + 1;
			SIZE_TYPE del = matrix[i-1][j].first + 1;
			SIZE_TYPE substitution_cost = ignore_uncertain_nucleotides ? !ignored_nucleotide_pair(s1[i - 1], s2[j - 1]) : (s1[i - 1] == s2[j - 1] ? 0 : 1);
			SIZE_TYPE sub = matrix[i - 1][j - 1].first + substitution_cost;
			if(ins < del) {
				if(sub < ins) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(sub, substitution_cost ? static_cast<SIZE_TYPE>(EditOperation::Type::SUBSTITUTE) : SAME);
				} else {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(ins, static_cast<SIZE_TYPE>(EditOperation::Type::INSERT));
				}
			} else {
				if(sub < del) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(sub, substitution_cost ? static_cast<SIZE_TYPE>(EditOperation::Type::SUBSTITUTE) : SAME);
				} else {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(del, static_cast<SIZE_TYPE>(EditOperation::Type::DELETE));
				}
			}
			if(matrix[i][j].first < m) {
				m = matrix[i][j].first;
			}
		}
		if(m > limit) {
			return std::nullopt;
		}
	}

	std::vector<EditOperation> result;
	SimpleIntsPair position(s1_size, s2_size);
	result.reserve(matrix[position.first][position.second].first);
	bool loop = true;
	while(loop) {
		SimpleIntsPair r = matrix[position.first][position.second];
		switch(r.second) {
		case BEGIN:
			loop = false;
			break;
		case SAME:
			--position.first;
			--position.second;
			break;
		case static_cast<SIZE_TYPE>(EditOperation::Type::INSERT):
			--position.second;
			result.emplace_back(position.first, s2[position.second], EditOperation::Type::INSERT);
			break;
		case static_cast<SIZE_TYPE>(EditOperation::Type::DELETE):
			--position.first;
			result.emplace_back(position.first, s1[position.first], EditOperation::Type::DELETE);
			break;
		case static_cast<SIZE_TYPE>(EditOperation::Type::SUBSTITUTE):
			--position.first;
			--position.second;
			result.emplace_back(position.first, s2[position.second], EditOperation::Type::SUBSTITUTE);
			break;
		default:
			throw std::logic_error("ASDF!");
		}
	}
	std::reverse(result.begin(), result.end());

	return result;
}

template<typename SIZE_TYPE, typename CHAR_TYPE>
static std::pair<SIZE_TYPE, SIZE_TYPE> find_breaks_once(const CHAR_TYPE* string1, SIZE_TYPE len1, const CHAR_TYPE* string2, SIZE_TYPE len2) {
	std::pair<SIZE_TYPE, SIZE_TYPE> breaks = {len1 / 2, len2 / 2 - SKIP_LEN / 2};
	for(SIZE_TYPE i = 0; i < SKIP_TRIES; ++i) {
		for(SIZE_TYPE j = 0; j < SKIP_LEN; ++j) {
			if(memcmp(string1 + breaks.first, string2 + breaks.second, SKIP_LEN * sizeof(CHAR_TYPE)) == 0) {
				return breaks;
			}
			++breaks.second;
		}
		breaks.first += SKIP_LEN;
	}
	throw std::runtime_error("Can't find a break.");
}

template<typename SIZE_TYPE, typename CHAR_TYPE>
static void find_breaks_do(const CHAR_TYPE* string1, SIZE_TYPE len1, const CHAR_TYPE* string2, SIZE_TYPE len2, std::vector<std::pair<std::pair<const CHAR_TYPE*, SIZE_TYPE>, std::pair<const CHAR_TYPE*, SIZE_TYPE>>>& breaks_vec) {
	if(len1 > MAX_PIECE_LEN && len2 > MAX_PIECE_LEN) {
		std::pair<SIZE_TYPE, SIZE_TYPE> breaks = find_breaks_once(string1, len1, string2, len2);
		find_breaks_do(string1, breaks.first, string2, breaks.second, breaks_vec);
		find_breaks_do<SIZE_TYPE>(string1 + breaks.first, len1 - breaks.first, string2 + breaks.second, len2 - breaks.second, breaks_vec);
	} else {
		breaks_vec.emplace_back(std::make_pair(string1, len1), std::make_pair(string2, len2));
	}
}

template<typename SIZE_TYPE, typename CHAR_TYPE>
static std::vector<std::pair<std::pair<const CHAR_TYPE*, SIZE_TYPE>, std::pair<const CHAR_TYPE*, SIZE_TYPE>>> find_breaks(const CHAR_TYPE* string1, SIZE_TYPE len1, const CHAR_TYPE* string2, SIZE_TYPE len2) {
	std::vector<std::pair<std::pair<const CHAR_TYPE*, SIZE_TYPE>, std::pair<const CHAR_TYPE*, SIZE_TYPE>>> breaks;
	find_breaks_do(string1, len1, string2, len2, breaks);
	return breaks;
}

static std::optional<std::vector<EditOperation>> edit_distance_int(const char *string1, uint16_t len1, const char *string2, uint16_t len2, uint16_t limit, bool ignore_uncertain_nucleotides) {
	std::vector<int> s1(len1);
	std::vector<int> s2(len2);
	for(uint16_t i = 0; i < len1; ++i) {
		s1[i] = string1[i];
	}
	for(uint16_t i = 0; i < len2; ++i) {
		s2[i] = string2[i];
	}
	auto breaks = find_breaks(s1.data(), len1, s2.data(), len2);
	std::vector<EditOperation> result;
	uint16_t offset = 0;
	for(const auto&[f, s] : breaks) {
		auto r = _edit_distance(f.first, f.second, s.first, s.second, limit, ignore_uncertain_nucleotides);
		if(!r) {
			return std::nullopt;
		}
		for(const EditOperation& op : *r) {
			result.emplace_back(op.position + offset, op.arg, op.type);
		}
		offset += f.second;
		if(result.size() > limit) {
			return std::nullopt;
		}
	}
	return result;
}

static PyObject* py_get_mutations(PyObject*, PyObject *args) {
    const char* string1;
    const char* string2;
    int ignore_uncertain_nucleotides;
    if (!PyArg_ParseTuple(args, "ssp", &string1, &string2, &ignore_uncertain_nucleotides))
        return nullptr;
    std::optional<std::vector<EditOperation>> r;
    try {
	r = edit_distance_int(string1, strlen(string1), string2, strlen(string2), (uint16_t) -1U, ignore_uncertain_nucleotides);
    } catch (const std::exception& e) {
	PyErr_SetString(PyExc_RuntimeError, e.what());
	return nullptr;
    } catch (...) {
	PyErr_SetString(PyExc_RuntimeError, "Hell broke loose.");
	return nullptr;
    }
    if(!r) {
	Py_INCREF(Py_None);
	return Py_None;
    }
    PyObject* ret = PyList_New(0);
    if(!ret) {
	PyErr_SetString(PyExc_RuntimeError, "Failed to create a new list.");
	return nullptr;
    }
    for(const EditOperation& eo : *r) {
	const char* error = nullptr;
	bool err = true;
	PyObject *literal_position, *literal_arg, *literal_type, *position_value, *arg_value, *type_value;
	char buffer[2];
	buffer[1] = 0;
	PyObject* dict = PyDict_New();
	if(!dict) {
		error = "Failed to create a new dict.";
		goto err0;
	}
	literal_position = PyUnicode_FromString("position");
	if(!literal_position) {
		error = "Failed to create 'position' literal.";
		goto err1;
	}
	literal_arg = PyUnicode_FromString("arg");
	if(!literal_arg) {
		error = "Failed to create 'arg' literal.";
		goto err2;
	}
	literal_type = PyUnicode_FromString("type");
	if(!literal_type) {
		error = "Failed to create 'type' literal.";
		goto err3;
	}
	position_value = PyLong_FromLong(eo.position);
	if(!position_value) {
		error = "Failed to create position value.";
		goto err4;
	}
	buffer[0] = eo.arg;
	arg_value = PyUnicode_FromString(buffer);
	if(!arg_value) {
		error = "Failed to create arg value.";
		goto err5;
	}
	buffer[0] = static_cast<char>(eo.type);
	type_value = PyUnicode_FromString(buffer);
	if(!type_value) {
		error = "Failed to create type value.";
		goto err6;
	}
	if(PyDict_SetItem(dict, literal_position, position_value) || PyDict_SetItem(dict, literal_arg, arg_value) || PyDict_SetItem(dict, literal_type, type_value)) {
		error = "Failed to set dictionary.";
		goto err7;
	}
	err = false;
	err7:
	Py_DECREF(type_value);
	err6:
	Py_DECREF(arg_value);
	err5:
	Py_DECREF(position_value);
	err4:
	Py_DECREF(literal_type);
	err3:
	Py_DECREF(literal_arg);
	err2:
	Py_DECREF(literal_position);
	if(!err) {
		if(PyList_Append(ret, dict)) {
			Py_DECREF(dict);
			Py_DECREF(ret);
			return nullptr;
		} else {
			Py_DECREF(dict);
			continue;
		}
	}
	err1:
	Py_DECREF(dict);
	err0:
	Py_DECREF(ret);
	PyErr_SetString(PyExc_RuntimeError, error);
	return nullptr;
    }
    return ret;
}

static PyMethodDef methods[] = {
    {"mutations",  py_get_mutations, METH_VARARGS, "Get mutations between two strings."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "buttplug", /* name of module */
    nullptr, /* module documentation */
    -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods,
    nullptr,
    nullptr,
    0,
    nullptr
};

PyMODINIT_FUNC PyInit_buttplug() {
    return PyModule_Create(&module);
}
