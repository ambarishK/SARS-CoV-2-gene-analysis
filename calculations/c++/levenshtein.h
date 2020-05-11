#pragma once
#include <vector>
#include <limits>
#include <cstring>
#include <cstdlib>

template<typename SIZE_TYPE, typename CHAR_TYPE>
// https://github.com/miohtama/python-Levenshtein/blob/master/Levenshtein.c
SIZE_TYPE lev_edit_distance(const CHAR_TYPE *string1, SIZE_TYPE len1, const CHAR_TYPE *string2, SIZE_TYPE len2) {
	SIZE_TYPE i;
	SIZE_TYPE *row;
	SIZE_TYPE *end;
	SIZE_TYPE half;

	while (len1 > 0 && len2 > 0 && *string1 == *string2) {
		len1--;
		len2--;
		string1++;
		string2++;
	}

	while (len1 > 0 && len2 > 0 && string1[len1-1] == string2[len2-1]) {
		len1--;
		len2--;
	}

	if (len1 == 0)
		return len2;
	if (len2 == 0)
		return len1;

	if (len1 > len2) {
		std::swap(len1, len2);
		std::swap(string1, string2);
	}
	if (len1 == 1) {
		return len2 - (memchr(string2, *string1, len2) != NULL);
	}
	len1++;
	len2++;
	half = len1 >> 1;

	row = (SIZE_TYPE*)malloc(len2*sizeof(SIZE_TYPE));
	if (!row)
		return (SIZE_TYPE)(-1);
	end = row + len2 - 1;
	for (i = 0; i < len2 - half; i++)
		row[i] = i;

	row[0] = len1 - half - 1;
	for (i = 1; i < len1; i++) {
		SIZE_TYPE *p;
		const CHAR_TYPE char1 = string1[i - 1];
		const CHAR_TYPE *char2p;
		SIZE_TYPE D, x;
		if (i >= len1 - half) {
			SIZE_TYPE offset = i - (len1 - half);
			SIZE_TYPE c3;

			char2p = string2 + offset;
			p = row + offset;
			c3 = *(p++) + (char1 != *(char2p++));
			x = *p;
			x++;
			D = x;
			if (x > c3)
				x = c3;
			*(p++) = x;
		}
		else {
			p = row + 1;
			char2p = string2;
			D = x = i;
		}
		if (i <= half + 1)
			end = row + len2 + i - half - 2;
		while (p <= end) {
			SIZE_TYPE c3 = --D + (char1 != *(char2p++));
			x++;
			if (x > c3)
				x = c3;
			D = *p;
			D++;
			if (x > D)
				x = D;
			*(p++) = x;
		}
		if (i <= half) {
			SIZE_TYPE c3 = --D + (char1 != *char2p);
			x++;
			if (x > c3)
				x = c3;
			*p = x;
		}
	}

	i = *end;
	free(row);
	return i;
}

//https://github.com/tdebatty/java-string-similarity/blob/master/src/main/java/info/debatty/java/stringsimilarity/Levenshtein.java
template<typename SIZE_TYPE, typename CHAR_TYPE>
SIZE_TYPE lev_edit_distance(const CHAR_TYPE* s1, SIZE_TYPE s1_size, const CHAR_TYPE* s2, SIZE_TYPE s2_size, const SIZE_TYPE limit) {
	while (s1_size > 0 && s2_size > 0 && *s1 == *s2) {
		s1_size--;
		s2_size--;
		s1++;
		s2++;
	}

	while (s1_size > 0 && s2_size > 0 && s1[s1_size-1] == s2[s2_size-1]) {
		s1_size--;
		s2_size--;
	}

	if (s1_size == 0) {
		return s2_size;
	}

	if (s2_size == 0) {
		return s1_size;
	}

	// create two work vectors of integer distances
	std::vector<SIZE_TYPE> v0(s2_size + 1);
	std::vector<SIZE_TYPE> v1(s2_size + 1);

	// initialize v0 (the previous row of distances)
	// this row is A[0][i]: edit distance for an empty s
	// the distance is just the number of characters to delete from t
	for (SIZE_TYPE i = 0; i < (SIZE_TYPE) v0.size(); ++i) {
		v0[i] = i;
	}

	for (SIZE_TYPE i = 0; i < s1_size; ++i) {
		// calculate v1 (current row distances) from the previous row v0
		// first element of v1 is A[i+1][0]
		//   edit distance is delete (i+1) chars from s to match empty t
		v1[0] = i + 1;

		SIZE_TYPE minv1 = v1[0];

		// use formula to fill in the rest of the row
		for (SIZE_TYPE j = 0; j < s2_size; ++j) {
			SIZE_TYPE cost = 1;
			if (s1[i] == s2[j]) {
				cost = 0;
			}
			v1[j + 1] = std::min(
					v1[j] + 1,			  // Cost of insertion
					std::min(
							v0[j + 1] + 1,  // Cost of remove
							v0[j] + cost)); // Cost of substitution

			minv1 = std::min(minv1, v1[j + 1]);
		}

		if (minv1 >= limit) {
			return limit;
		}

		// Flip references to current and previous row
		std::swap(v0, v1);

	}

	return v0[s2_size];
}

template<typename INT_TYPE = int>
INT_TYPE lev_edit_distance_int(const char *string1, INT_TYPE len1, const char *string2, INT_TYPE len2, INT_TYPE limit = std::numeric_limits<INT_TYPE>::max()) {
	std::vector<INT_TYPE> s1(len1);
	std::vector<INT_TYPE> s2(len2);
	for(INT_TYPE i = 0; i < len1; ++i) {
		s1[i] = string1[i];
	}
	for(INT_TYPE i = 0; i < len2; ++i) {
		s2[i] = string2[i];
	}
	if(limit == std::numeric_limits<INT_TYPE>::max()) {
		return lev_edit_distance(s1.data(), len1, s2.data(), len2);
	} else {
		return lev_edit_distance(s1.data(), len1, s2.data(), len2, limit);
	}
}
