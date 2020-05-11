#include "edit_distance.h"

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
static std::optional<std::vector<std::variant<SIZE_TYPE, std::pair<SIZE_TYPE, CHAR_TYPE>>>> _edit_distance(const CHAR_TYPE* s1, SIZE_TYPE s1_size, const CHAR_TYPE* s2, SIZE_TYPE s2_size, const SIZE_TYPE limit) {
	std::vector<std::vector<SimpleIntsPair<SIZE_TYPE>>> matrix(s1_size + 1, std::vector<SimpleIntsPair<SIZE_TYPE>>(s2_size + 1, SimpleIntsPair<SIZE_TYPE>(0, 0)));
	static constexpr const SIZE_TYPE SAME = 3, DELETE = 2, INSERT = 1, BEGIN = 0;
	for(SIZE_TYPE i = 1; i <= s2_size; ++i) {
		matrix[0][i] = SimpleIntsPair<SIZE_TYPE>(i, INSERT);
	}

	for(SIZE_TYPE i = 1; i <= s1_size; ++i) {
		matrix[i][0] = SimpleIntsPair<SIZE_TYPE>(i, DELETE);
		SIZE_TYPE m = i;
		for(SIZE_TYPE j = 1; j <= s2_size; ++j) {
			SIZE_TYPE ins = matrix[i][j-1].first + 1;
			SIZE_TYPE del = matrix[i-1][j].first + 1;
			if(ins < del) {
				if(s1[i - 1] == s2[j - 1] && matrix[i - 1][j - 1].first < ins) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(matrix[i - 1][j - 1].first, SAME);
				} else {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(ins, INSERT);
				}
			} else {
				if(s1[i - 1] == s2[j - 1] && matrix[i - 1][j - 1].first < del) {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(matrix[i - 1][j - 1].first, SAME);
				} else {
					matrix[i][j] = SimpleIntsPair<SIZE_TYPE>(del, DELETE);
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

	std::vector<std::variant<SIZE_TYPE, std::pair<SIZE_TYPE, CHAR_TYPE>>> result;
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
		case INSERT:
			--position.second;
			result.emplace_back(std::make_pair(position.first, s2[position.second]));
			break;
		case DELETE:
			--position.first;
			result.emplace_back(position.first);
			break;
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

std::optional<std::vector<std::variant<uint16_t, std::pair<uint16_t, char>>>> edit_distance_int(const char *string1, uint16_t len1, const char *string2, uint16_t len2, uint16_t limit) {
	std::vector<int> s1(len1);
	std::vector<int> s2(len2);
	for(uint16_t i = 0; i < len1; ++i) {
		s1[i] = string1[i];
	}
	for(uint16_t i = 0; i < len2; ++i) {
		s2[i] = string2[i];
	}
	auto breaks = find_breaks(s1.data(), len1, s2.data(), len2);
	std::vector<std::variant<uint16_t, std::pair<uint16_t, char>>> result;
	uint16_t offset = 0;
	for(const auto&[f, s] : breaks) {
		auto r = _edit_distance(f.first, f.second, s.first, s.second, limit);
		if(!r) {
			return std::nullopt;
		}
		for(const auto& op : *r) {
			if(std::holds_alternative<uint16_t>(op)) {
				result.emplace_back(std::get<uint16_t>(op) + offset);
			} else {
				auto p = std::get<std::pair<uint16_t, int>>(op);
				result.emplace_back(std::make_pair(p.first + offset, (char) p.second));
			}
		}
		offset += f.second;
		if(result.size() > limit) {
			return std::nullopt;
		}
	}
	return result;
}
