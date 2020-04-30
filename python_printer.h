#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <type_traits>
#include <memory>
#include <unordered_map>
#include <utility>
#include <functional>

#include "string_tools.h"

template<typename... Ts>
void print_as_python(std::ostream& o, const std::tuple<Ts...>& tp);

void print_as_python(std::ostream& o, const std::string& string);

void print_as_python(std::ostream& o, char c);

template<typename T>
void print_as_python(std::ostream& o, const std::vector<T>& vec);

template<typename T, typename std::enable_if_t<(std::is_arithmetic_v<T> && !std::is_same_v<T, char>), int> = 0>
void print_as_python(std::ostream& o, const T& arg) {
	o << arg;
}

void print_as_python_tuple_(std::ostream&);

template<typename T, typename... Ts>
void print_as_python_tuple_(std::ostream& o, const T& arg, const Ts&... args) {
	print_as_python(o, arg);
	o << ", ";
	print_as_python_tuple_(o, args...);
}

template<typename... Ts>
void print_as_python(std::ostream& o, const std::tuple<Ts...>& tp) {
	o << '(';
	std::apply([&o](const auto&... args){print_as_python_tuple_(o, args...);}, tp);
	o << ')';
}

template<typename T>
void print_as_python(std::ostream& o, const std::unique_ptr<T>& ptr) {
	if(ptr) {
		print_as_python(o, *ptr);
	} else {
		o << "None";
	}
}

template<typename T>
void print_as_python(std::ostream& o, const std::shared_ptr<T>& ptr) {
	if(ptr) {
		print_as_python(o, *ptr);
	} else {
		o << "None";
	}
}

template<typename T>
void print_as_python(std::ostream& o, const std::vector<T>& vec) {
	o << '[';
	if(!vec.empty()) {
		print_as_python(o, vec[0]);
	}
	for(size_t i = 1; i < vec.size(); ++i) {
		o << ", ";
		print_as_python(o, vec[i]);
	}
	o << ']';
}

template<typename A, typename B>
void print_as_python(std::ostream& o, const std::pair<A, B>& p) {
	print_as_python(o, std::forward_as_tuple(p.first, p.second));
}

void print_as_python_dict_(std::ostream&);

template<typename K, typename V, typename... Ts>
void print_as_python_dict_(std::ostream& o, const K& key, const V& value, const Ts&... args) {
	o << ", ";
	print_as_python(o, key);
	o << ": ";
	print_as_python(o, value);
	print_as_python_dict_(o, args...);
}

template<typename K, typename V, typename... Ts>
void print_as_python_dict(std::ostream& o, const K& key, const V& value, const Ts&... args) {
	o << "{";
	print_as_python(o, key);
	o << ": ";
	print_as_python(o, value);
	print_as_python_dict_(o, args...);
	o << "}";
}

template<typename T>
struct python_loader;

template<typename K, typename V>
struct python_loader<std::unordered_map<K, V>> {
	static std::unordered_map<K, V> load(std::istream& input) {
		input >> std::ws;
		if(input.peek() != '{') {
			throw std::runtime_error("Incorrect input.");
		}
		input.get();
		std::unordered_map<K, V> map;
		while(true) {
			input >> std::ws;
			if(input.peek() == '}') {
				break;
			}
			K key = python_loader<K>::load(input);
			input >> std::ws;
			if(input.peek() != ':') {
				throw std::runtime_error("Incorrect input.");
			}
			input.get();
			V value = python_loader<K>::load(input);
			map.emplace(std::move(key, std::move(value)));
			input >> std::ws;
			if(input.peek() == ',') {
				input.get();
				continue;
			} else {
				break;
			}
		}
		input >> std::ws;
		if(input.peek() != '}') {
			throw std::runtime_error("Incorrect input.");
		}
		input.get();
		return map;
	}
};

template<>
struct python_loader<std::string> {
	static std::string load(std::istream& input);
};

template<>
struct python_loader<char> {
	static char load(std::istream& input);
};

template<typename T>
struct python_loader_int {
	static T load(std::istream& input) {
		T num;
		input >> num;
		if(input.fail() || input.bad()) {
			throw std::runtime_error("Incorrect input.");
		}
		return num;
	}
};

template<>
struct python_loader<int> : python_loader_int<int> {};

template<>
struct python_loader<unsigned> : python_loader_int<unsigned> {};

template<>
struct python_loader<long> : python_loader_int<long> {};

template<>
struct python_loader<unsigned long> : python_loader_int<unsigned long> {};

template<>
struct python_loader<long long> : python_loader_int<long long> {};

template<>
struct python_loader<unsigned long long> : python_loader_int<unsigned long long> {};

template<typename T>
struct python_loader<std::vector<T>> {
	static std::vector<T> load(std::istream& input) {
		input >> std::ws;
		if(input.peek() != '[') {
			throw std::runtime_error("Incorrect input.");
		}
		input.get();
		std::vector<T> vec;
		while(true) {
			input >> std::ws;
			if(input.peek() == ']') {
				break;
			}
			vec.push_back(python_loader<T>::load(input));
			input >> std::ws;
			if(input.peek() == ',') {
				input.get();
				continue;
			} else {
				break;
			}
		}
		input >> std::ws;
		if(input.peek() != ']') {
			throw std::runtime_error("Incorrect input.");
		}
		input.get();
		return vec;
	}
};


template<typename... Ts>
struct python_loader<std::tuple<Ts...>> {
	static std::tuple<Ts...> load(std::istream& input) {
		input >> std::ws;
		if(input.peek() != '(') {
			throw std::runtime_error("Incorrect input.");
		}
		input.get();
		return get_tuple_contents_dispatch<Ts...>(input);
	}
private:
	template<typename U, typename... Us>
	static std::tuple<U, Us...> get_tuple_contents(std::istream& input) {
		U element = python_loader<U>::load(input);
		input >> std::ws;
		if(input.peek() == ',') {
			input.get();
			return std::tuple_cat(std::make_tuple<U>(std::move(element)), get_tuple_contents_dispatch<Us...>(input));
		} else {
			if constexpr (sizeof...(Us) == 0) {
				return std::make_tuple<U>(std::move(element));
			} else {
				throw std::runtime_error("Incorrect input.");
			}
		}
	}

	template<typename... Us>
	static std::tuple<Us...> get_tuple_contents_dispatch(std::istream& input) {
		if constexpr(sizeof...(Us) > 0) {
			return get_tuple_contents<Us...>(input);
		} else {
			input >> std::ws;
			if(input.peek() != ')') {
				throw std::runtime_error("Incorrect input.");
			}
			input.get();
			return std::tuple<>();
		}
	}
};

template<typename A, typename B>
struct python_loader<std::pair<A, B>> {
	static std::pair<A, B> load(std::istream& input) {
		auto tuple = python_loader<std::tuple<A, B>>::load(input);
		return std::make_pair<A, B>(std::move(std::get<0>(tuple)), std::move(std::get<1>(tuple)));
	}
};

template<size_t I, typename... Ts>
struct get_type {
	using Type = std::remove_reference_t<decltype(std::get<I>(std::declval<std::tuple<Ts...>>()))>;
};

template<typename... Ts>
struct python_loader_struct {
	template<typename... Ns>
	class with_names {
		std::tuple<Ts...> tuple;
		static_assert(sizeof...(Ts) == sizeof...(Ns));

		template<size_t I>
		static void make_ids_emplace(std::unordered_map<std::string, size_t>& map) {
				map.emplace(get_type<I, Ns...>::Type::str(), I);
			if constexpr(I + 1 < sizeof...(Ns)) {
				make_ids_emplace<I + 1>(map);
			}
		}

		static std::unordered_map<std::string, size_t> make_ids() {
			std::unordered_map<std::string, size_t> map;
			if constexpr (sizeof...(Ns) > 0) {
				make_ids_emplace<0>(map);
			}
			return map;
		}

		template<size_t I = 0>
		static void parse(size_t num, std::tuple<std::unique_ptr<Ts>...>& tuple, std::istream& input) {
			if constexpr (I >= sizeof...(Ns)) {
				throw std::runtime_error("Incorrect input.");
			} else {
				if(num == I) {
					using Type = typename get_type<I, Ts...>::Type;
					std::get<I>(tuple) = std::unique_ptr<Type>(new Type(python_loader<Type>::load(input)));
				} else {
					parse<I + 1>(num, tuple, input);
				}
			}
		}

		template<size_t... I>
		static std::tuple<Ts...> unwrap_tuple_(std::tuple<std::unique_ptr<Ts>...>&& tuple, std::integer_sequence<size_t, I...>) {
			return std::tuple<Ts...>(std::move(*std::get<I>(tuple))...);
		}

		static std::tuple<Ts...> unwrap_tuple(std::tuple<std::unique_ptr<Ts>...>&& tuple) {
			return unwrap_tuple_(std::move(tuple), std::make_index_sequence<sizeof...(Ts)>());
		}
	public:
		std::tuple<Ts...>& as_tuple() {
			return tuple;
		}

		const std::tuple<Ts...>& as_tuple() const {
			return tuple;
		}
		with_names() = delete;
		with_names(std::tuple<Ts...>&& tuple) : tuple(std::move(tuple)) {}
		static with_names load(std::istream& input) {
			std::unordered_map<std::string, size_t> ids = make_ids();
			std::tuple<std::unique_ptr<Ts>...> parsed;
			input >> std::ws;
			if(input.peek() != '{') {
				throw std::runtime_error("Incorrect input.");
			}
			input.get();
			while(true) {
				input >> std::ws;
				if(input.peek() == '}') {
					break;
				}
				size_t id = ids.at(python_loader<std::string>::load(input));
				input >> std::ws;
				if(input.peek() != ':') {
					throw std::runtime_error("Incorrect input.");
				}
				input.get();
				parse(id, parsed, input);
				input >> std::ws;
				if(input.peek() == ',') {
					input.get();
					continue;
				} else {
					break;
				}
			}
			input >> std::ws;
			if(input.peek() != '}') {
				throw std::runtime_error("Incorrect input.");
			}
			input.get();
			return unwrap_tuple(std::move(parsed));
		}
	};
	template<typename... Ns>
	static with_names<Ns...> with_names_values(const Ns&...);
};

template<typename T>
struct python_loader {
	static T load(std::istream& input) {
		return T::load(input);
	}
};
