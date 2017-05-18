#ifndef ENCODINGS_H_
#define ENCODINGS_H_

typedef uint64_t cell;
typedef char base;

const char COMPACT_TO_IUPAC[] = {1, 2, 4, 8};
const char IUPAC_TO_COMPACT[] = {-1, 0, 1, -1,
                                 2, -1, -1, -1,
                                 3, -1, -1, -1,
                                 -1, -1, -1, -1};
const char COMPACT_TO_CHAR[] = {'a', 'c', 'g', 't'};
const char IUPAC_TO_CHAR[] = {'-', 'a', 'c', 'm',
                              'g', 'r', 's', 'v',
                              't', 'w', 'y', 'h',
                              'k', 'd', 'b', 'n'};
const char IUPAC_RC[] = {0, 8, 4, 12,
                         2, 10, 6, 14,
                         1, 9, 5, 13,
                         3, 11, 7, 15};
const char COMPACT_RC[] = {3, 2, 1, 0};

class BadCharException : public std::exception {
private:
    std::string msg;
public:
    BadCharException(std::string encoding, char chr) {
        msg = "Bad nucleotide for " + encoding + ": " + chr;
    };
    virtual const char * what() const noexcept {
        return msg.c_str();
    };
};

inline char CHAR_TO_IUPAC(char c) {
    switch (c) {
    case '-': return 0;
    case 'a': case 'A': return 1;
    case 'c': case 'C': return 2;
    case 'm': case 'M': return 3;
    case 'g': case 'G': return 4;
    case 'r': case 'R': return 5;
    case 's': case 'S': return 6;
    case 'v': case 'V': return 7;
    case 't': case 'T': return 8;
    case 'w': case 'W': return 9;
    case 'y': case 'Y': return 10;
    case 'h': case 'H': return 11;
    case 'k': case 'K': return 12;
    case 'd': case 'D': return 13;
    case 'b': case 'B': return 14;
    case 'n': case 'N': return 15;
    }
    throw(BadCharException("IUPAC", c));
}

inline char CHAR_TO_COMPACT(char c) {
    switch(c) {
    case 'a': case 'A': return 0;
    case 'c': case 'C' :return 1;
    case 'g': case 'G' :return 2;
    case 't': case 'T': return 3;
    }
    throw(BadCharException("Compact", c));
}

const unsigned COMPACT_SIZE = 2;
const cell COMPACT_MASK = 3;
const unsigned COMPACTS_PER_CELL = sizeof(cell) * 8 / COMPACT_SIZE;
const unsigned IUPAC_SIZE = 4;
const cell IUPAC_MASK = 0xf;
const unsigned IUPAC_PER_CELL = sizeof(cell) * 8 / IUPAC_SIZE;

inline base getIUPAC(cell * buf, unsigned pos) {
    return (buf[pos / IUPAC_PER_CELL] >> (pos % IUPAC_PER_CELL * IUPAC_SIZE)) & IUPAC_MASK;
};

inline void putIUPAC(cell * buf, unsigned pos, cell val) {
    unsigned shift =  pos % IUPAC_PER_CELL * IUPAC_SIZE;
    buf[pos / IUPAC_PER_CELL] &= ~(IUPAC_MASK << shift);
    buf[pos / IUPAC_PER_CELL] |= (val & IUPAC_MASK) << shift;
};

inline base getCompact(cell * buf, unsigned pos) {
    return (buf[pos / COMPACTS_PER_CELL] >> (pos % COMPACTS_PER_CELL * COMPACT_SIZE)) & COMPACT_MASK;
};

inline void putCompact(cell * buf, unsigned pos, cell val) {
    unsigned shift = pos % COMPACTS_PER_CELL * COMPACT_SIZE;
    buf[pos / COMPACTS_PER_CELL] &= ~(COMPACT_MASK << shift);
    buf[pos / COMPACTS_PER_CELL] |= (val & COMPACT_MASK) << shift;
};

#endif
