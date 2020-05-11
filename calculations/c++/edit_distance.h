#pragma once

#include <vector>
#include <variant>
#include <limits>
#include <tuple>
#include <optional>
#include <string>

std::optional<std::vector<std::variant<uint16_t, std::pair<uint16_t, char>>>> edit_distance_int(const char *string1, uint16_t len1, const char *string2, uint16_t len2, uint16_t limit = std::numeric_limits<uint16_t>::max());
