if [ "$1" = "debug" ]; then
	E_FLAGS="-O0 -DDEBUG=1"
else
	E_FLAGS="-O3 -DDEBUG=0"
fi
g++ -g -std=c++17 edit_distance.cpp python_printer.cpp paths.cpp string_tools.cpp findgenemutations.cpp -o findgenemutations $E_FLAGS -Werror -Wall -Wextra -march=native -ltbb
