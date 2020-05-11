#pragma once
#include <string>
#include <tuple>
#include <vector>

void rstrip(std::string& string);

std::string concatenate_strings(const std::vector<std::string>& strings, size_t reserve_size);

template <char... Cs>
struct string_literal {
	static const char* str() {
		static const char value[sizeof...(Cs) + 1] = {Cs..., '\0'};
		return value;
	}
};

template <typename CharT, CharT... Cs>
constexpr string_literal<Cs...> operator ""_sl() {
	return {};
}
