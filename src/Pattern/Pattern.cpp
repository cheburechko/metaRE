/*
 * Pattern.cpp
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#include "Pattern.h"
#include <algorithm>
#include "PatternWrongSizeException.hpp"

using namespace std;

const map<char, const vector<unsigned> > Pattern::lettersToInt = {
		{'a', {0}},
		{'c', {1}},
		{'g', {2}},
		{'t', {3}},
		{'w', {0, 3}},
		{'s', {1, 2}},
		{'m', {0, 1}},
		{'k', {2, 3}},
		{'r', {0, 2}},
		{'y', {1, 3}},
		{'b', {1,2,3}},
		{'d', {0,2,3}},
		{'h', {0,1,3}},
		{'v', {0,1,2}},
		{'n', {0,1,2,3}}
};

Pattern::Pattern() :
		mask(0),
		maskSize(0),
		patternSize(0)
{}

Pattern::Pattern(std::string pattern) :
	Pattern() {
	init(pattern);
}

void Pattern::recursiveAdd(unsigned start, string end) {
	for (auto iter = lettersToInt.at(end[0]).begin();
			iter != lettersToInt.at(end[0]).end();
			iter++) {
		unsigned next = start * 4 + *iter;
		if (end.length() == 1) {
			markPattern(next);
		} else {
			recursiveAdd(next, end.substr(1));
		}
	}
}

void Pattern::add(std::string pattern) {
	if (maskSize == 0) {
		init(pattern);
	} else {
		if (pattern.length() != patternSize) {
			throw PatternWrongSizeException("Attempt to add kmer to pattern of different size");
		}
		transform(pattern.begin(), pattern.end(), pattern.begin(), ::tolower);
		recursiveAdd(0, pattern);
	}
}

bool Pattern::check(unsigned kmer) const {
	if ((kmer >> MASK_SHIFT) >= maskSize) {
		return false;
	}
	return mask[kmer >> MASK_SHIFT] & (((unsigned long) 1) << (kmer &  MASK_FILTER));
}

void Pattern::markPattern(unsigned kmer) {
	if ((kmer >> MASK_SHIFT) < maskSize) {
		mask[kmer >> MASK_SHIFT] |= ((unsigned long) 1) << (kmer &  MASK_FILTER);
	}
}

void Pattern::init(std::string pattern) {
	patternSize = pattern.length();
	maskSize = 1 << (2*patternSize - MASK_SHIFT);

	mask.resize(maskSize, 0);
	add(pattern);
}

std::vector<unsigned> Pattern::getKmers() const {
	std::vector<unsigned> result;
	for (unsigned kmer = 0; kmer < maskSize << MASK_SHIFT; kmer++) {
		if (check(kmer)) {
			result.push_back(kmer);
		}
	}
	return result;
}
