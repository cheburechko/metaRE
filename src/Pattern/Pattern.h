/*
 * Pattern.h
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#ifndef PATTERN_H_
#define PATTERN_H_

#include <string>
#include <map>
#include <vector>

/* TODO
 * Make this an interface with 2 implementations:
 * FixedSizePattern, VariableSizePattern
 * Add checks in FixedSizePattern on pattern length.
 */
/**
 * Each possible string has a bit assigned to it
 * E.g. TGT|CTC
 * Last 3 letters are used to find a bit in int64
 * First letters - to find position in the array of int64.
 */
class Pattern {
private:
	const unsigned MASK_SHIFT = 6, MASK_FILTER = 63;
	unsigned long maskSize, patternSize;
	std::vector<unsigned long> mask;
	static const std::map<char, const std::vector<unsigned> > lettersToInt;
	void recursiveAdd(unsigned start, std::string end);
	void markPattern(unsigned kmer);
	void init(std::string pattern);
public:
	Pattern();
	Pattern(std::string pattern);
	void add(std::string pattern);
	std::vector<unsigned> getKmers() const;
	bool check(unsigned kmer) const;
	virtual ~Pattern() {};
};

#endif /* PATTERN_H_ */
