#include "edit_distance.h"

#include <cstring>
#include <algorithm>
#include <map>
#include <set>

using namespace std::string_literals;

void print_as_python(std::ostream& o, const EditOperation& eo) {
	print_as_python_dict(o, "position"s, eo.position, "arg"s, eo.arg, "type"s, eo.type);
}

EditOperation python_loader<EditOperation>::load(std::istream& input) {
	auto p = decltype(python_loader_struct<uint16_t, char, EditOperation::Type>::with_names_values("position"_sl, "arg"_sl, "type"_sl))::load(input);
	auto& t = p.as_tuple();
	return EditOperation(std::get<0>(t), std::get<1>(t), std::get<2>(t));
}

void print_as_python(std::ostream& o, const EditOperation::Type& t) {
	print_as_python(o, static_cast<char>(t));
}

EditOperation::Type python_loader<EditOperation::Type>::load(std::istream& input) {
	return static_cast<EditOperation::Type>(python_loader<char>::load(input));
}

EditOperation::EditOperation(uint16_t position, char arg, Type type) : position(position), arg(arg), type(type) {}

namespace {
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
	size_t hash() const {
		return reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(*this);
	}
	SimpleIntsPair& operator=(const SimpleIntsPair& other) {
		reinterpret_cast<typename TwiceSizeInt<T>::Type&>(*this) = reinterpret_cast<const typename TwiceSizeInt<T>::Type&>(other);
		return *this;
	}
};
}

namespace std {
template<typename T>
struct hash<SimpleIntsPair<T>> {
	size_t operator()(const SimpleIntsPair<T>& k) const {
		return k.hash();
	}
};
}

template<typename SIZE_TYPE, typename CHAR_TYPE>
static std::optional<std::vector<EditOperation>> _edit_distance(const CHAR_TYPE* s1, SIZE_TYPE s1_size, const CHAR_TYPE* s2, SIZE_TYPE s2_size, const SIZE_TYPE limit) {
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
			SIZE_TYPE sub = matrix[i - 1][j - 1].first + (s1[i - 1] == s2[j - 1] ? 0 : 1);
			if(ins < del) {
				if(sub < ins) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(sub, s1[i - 1] == s2[j - 1] ? SAME : static_cast<SIZE_TYPE>(EditOperation::Type::SUBSTITUTE));
				} else {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(ins, static_cast<SIZE_TYPE>(EditOperation::Type::INSERT));
				}
			} else {
				if(sub < del) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(sub, s1[i - 1] == s2[j - 1] ? SAME : static_cast<SIZE_TYPE>(EditOperation::Type::SUBSTITUTE));
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


static void assert_throw(bool a) {
	if(!a) {
		throw std::logic_error("Assertion failed.");
	}
}

std::string rebuild(const std::vector<EditOperation>& eos, const char *string1, uint16_t len1) {
	std::string str;
	std::map<uint16_t, std::vector<EditOperation>> e;
	for(const EditOperation& eo : eos) {
		e[eo.position].push_back(eo);
	}
	std::set<uint16_t> deleted;
	for(uint16_t i = 0; i < len1; ++i) {
		if(e.count(i)) {
			for(const EditOperation& eo: e.at(i)) {
				switch(eo.type) {
				case EditOperation::Type::DELETE:
					assert_throw(string1[i] == eo.arg);
					assert_throw(deleted.emplace(i).second);
					break;
				case EditOperation::Type::SUBSTITUTE:
					assert_throw(deleted.emplace(i).second);
					// fall-through
				case EditOperation::Type::INSERT:
					str += eo.arg;
					break;
				default:
					assert_throw(false);
				}
			}
		}
		if(deleted.count(i) == 0) {
			str += string1[i];
		}
	}
	if(e.count(len1)) {
		for(const EditOperation& eo: e.at(len1)) {
			switch(eo.type) {
			case EditOperation::Type::INSERT:
				str += eo.arg;
				break;
			case EditOperation::Type::DELETE:
				// fall-through
			case EditOperation::Type::SUBSTITUTE:
				// fall-through
			default:
				assert_throw(false);
			}
		}
	}
	return str;
}

std::optional<std::vector<EditOperation>> edit_distance_int(const char *string1, uint16_t len1, const char *string2, uint16_t len2, uint16_t limit) {
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
		auto r = _edit_distance(f.first, f.second, s.first, s.second, limit);
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
