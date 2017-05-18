/*
 * Utils.cpp
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#include "Utils.h"
#include "UtilsNucleotideException.hpp"
#include "UtilsKmerTooBigException.hpp"


Utils::Utils() {
	// TODO Auto-generated constructor stub

}

unsigned Utils::reverseComplement(unsigned kmer, unsigned k) {
	if (kmer >= getKmersTotal(k)) {
		throw UtilsKmerTooBigException("Kmer is too big");
	}
	std::vector<char> digits(k);
	for (unsigned i = 0; i < k; i++) {
		digits[i] = kmer % 4;
		kmer /= 4;
	}
	int result = 0;
	for (unsigned i = 0; i < k; i++) {
		result = result * 4 + reverseComplementInt(digits[i]);
	}
	return result;
}

unsigned Utils::getKmersTotal(unsigned k) {
	return 1 << (k*2);
}

std::unique_ptr<unsigned[]> Utils::allReverseComplements(unsigned k) {
	unsigned kmersTotal = getKmersTotal(k);
	auto rcs = std::unique_ptr<unsigned[]>(new unsigned[kmersTotal]);
	for (unsigned i = 0; i < kmersTotal; i++) {
		rcs[i] = reverseComplement(i, k);
	}
	return rcs;
}

std::unique_ptr<unsigned[]> Utils::getPlusMinusMapping(unsigned k) {
	unsigned kmersTotal = getKmersTotal(k);
	auto plusAndMinus = std::unique_ptr<unsigned[]>(new unsigned[kmersTotal]);
	std::unique_ptr<unsigned[]> rcs = allReverseComplements(k);
	for (unsigned i = 0; i < kmersTotal; i++) {
		plusAndMinus[i] = kmersTotal;
	}
	for (unsigned i = 0; i < kmersTotal; i++) {
		if (plusAndMinus[rcs[i]] != kmersTotal) {
			plusAndMinus[i] = rcs[i];
		} else {
			plusAndMinus[i] = i;
		}
	}
	return plusAndMinus;
}

Utils::~Utils() {
	// TODO Auto-generated destructor stub
}

std::string Utils::reverseComplement(std::string kmer, bool toupper) {
	std::string result(kmer);
	for (unsigned i = 0; i < kmer.length(); i++) {
		result[i] = reverseComplementChar(kmer[kmer.length()-1-i], toupper);
	}
	return result;
}

unsigned Utils::stringToInt(std::string kmer) {
	unsigned result = 0;
	for (unsigned i = 0; i < kmer.length(); i++) {
		result = result * 4 + charToInt(kmer[i]);
	}
	return result;
}

std::string Utils::intToString(unsigned kmer, unsigned k, bool toupper) {
	if (kmer >= getKmersTotal(k)) {
		throw UtilsKmerTooBigException("Kmer is too big");
	}
	std::string result(k, 'A');
	for (int i = k-1; i >= 0 ; i--) {
		result[i] = intToChar(kmer % 4, toupper);
		kmer /= 4;
	}
	return result;
}

std::vector<std::string> Utils::generateAllKmerLabels(unsigned k, bool toupper) {
	unsigned kmersTotal = getKmersTotal(k);
	std::vector<std::string> elementLabels;
	for (unsigned kmer = 0; kmer < kmersTotal; kmer++) {
		elementLabels.push_back(Utils::intToString(kmer, k, toupper));
	}
	return elementLabels;
}
