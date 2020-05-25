#pragma once
#include <execution>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <mutex>

template<typename It1, typename It2, typename F>
void transform_with_progress(It1 input_begin, It1 input_end, It2 output_begin, F&& function, std::ostream& output) {
	std::atomic<size_t> done = 0;
	std::mutex output_mutex;
	size_t total = std::distance(input_begin, input_end);
	std::transform(
	#if DEBUG
	std::execution::seq,
	#else
	std::execution::par_unseq,
	#endif
	input_begin, input_end, output_begin, [function(std::move(function)), &done, &output, total, &output_mutex](auto& arg){
		auto result = function(arg);
		int d = done.fetch_add(1);
		if(d * 100 / total != (d - 1) * 100 / total) {
			const std::lock_guard<std::mutex> lock(output_mutex);
			output << std::to_string(d * 100 / total) + "%\n";
			output.flush();
		}
		return result;
	});
}

template<typename It1, typename F>
void for_each_with_progress(It1 input_begin, It1 input_end, F&& function, std::ostream& output) {
	std::atomic<size_t> done = 0;
	std::mutex output_mutex;
	size_t total = std::distance(input_begin, input_end);
	std::for_each(
	#if DEBUG
	std::execution::seq,
	#else
	std::execution::par_unseq,
	#endif
	input_begin, input_end, [function(std::move(function)), &done, &output, total, &output_mutex](auto& arg){
		function(arg);
		int d = done.fetch_add(1);
		if(d * 100 / total != (d - 1) * 100 / total) {
			const std::lock_guard<std::mutex> lock(output_mutex);
			output << std::to_string(d * 100 / total) + "%\n";
			output.flush();
		}
	});
}

std::ifstream open_file_i(const std::string& filename);

std::ofstream open_file_o(const std::string& filename);
