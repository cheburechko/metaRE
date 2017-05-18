/*
 * Utils/Utils.h
 *
 *  Created on: 21 авг. 2015 г.
 *      Author: thepunchy
 */

#ifndef Utils_H_
#define Utils_H_

#include <vector>
#include <string>
#include <memory>
#include "UtilsNucleotideException.hpp"
#include "UtilsKmerTooBigException.hpp"

class Utils {
private:
	Utils();
	virtual ~Utils();
public:
	static const unsigned MC_PLUS_MINUS = 1;

	inline static char reverseComplementChar(char nucleotide, bool toupper = false)  {
	    switch(nucleotide) {
	    case 'a': case 'A': return (toupper ? 'T' : 't');
	    case 'c': case 'C': return (toupper ? 'G' : 'g');
	    case 't': case 'T': return (toupper ? 'A' : 'a');
	    case 'g': case 'G': return (toupper ? 'C' : 'c');
	    case 'n': case 'N': return (toupper ? 'N' : 'n');
	    default: throw UtilsNucleotideException("wrong character");
	    }
	};
	inline static unsigned reverseComplementInt(unsigned nucleotide){
	    if (nucleotide > 3) {
	        throw UtilsNucleotideException("wrong unsigned");
	    }
	    return 3 - nucleotide;
	};
	inline static char intToChar(unsigned nucleotide, bool toupper = false)  {
	    switch(nucleotide) {
	    case 0: return (toupper ? 'A' : 'a');
	    case 1: return (toupper ? 'C' : 'c');
	    case 2: return (toupper ? 'G' : 'g');
	    case 3: return (toupper ? 'T' : 't');
	    default: throw UtilsNucleotideException("wrong unsigned ");
	    }
	};
	inline static unsigned charToInt(char nucleotide) {
	    switch(nucleotide) {
	    case 'a': case 'A': return 0;
	    case 'c': case 'C': return 1;
	    case 'g': case 'G': return 2;
	    case 't': case 'T': return 3;
	    default: return -1;
	    }
	};

	static unsigned reverseComplement(unsigned kmer, unsigned k);
	static std::string reverseComplement(std::string kmer, bool toupper = false);
	static unsigned stringToInt(std::string kmer);
	static std::string intToString(unsigned kmer, unsigned k, bool toupper = false);

	static unsigned getKmersTotal(unsigned k);
	/// Generates lookup table for kmers reverse complements
	static std::unique_ptr<unsigned[]> allReverseComplements(unsigned k);
	static std::unique_ptr<unsigned[]> getPlusMinusMapping(unsigned k);
	static std::vector<std::string> generateAllKmerLabels(unsigned k, bool toupper = false);
};

#endif /* Utils_H_ */
