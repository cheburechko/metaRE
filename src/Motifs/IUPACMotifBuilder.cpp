#include "IUPACMotifBuilder.h"
#include <iostream>

IUPACMotifBuilder::IUPACMotifBuilder(unsigned size) :
    length(size),
    bufSize((length - 1) / IUPAC_PER_CELL + 1),
    accumulated(0)
{
    buf = new cell[bufSize]();
    complement = new cell[bufSize]();
}

bool IUPACMotifBuilder::ready(unsigned size) const {
    return accumulated >= size && size <= length;
}

void IUPACMotifBuilder::skip() {
    accumulated = 0;
};

void IUPACMotifBuilder::put(unsigned nucleotide) {
    cell carry = nucleotide & IUPAC_MASK;
    for (unsigned i = 0; i < bufSize; i++) {
        cell nextCarry = (buf[i] >> ((IUPAC_PER_CELL-1) * IUPAC_SIZE)) & IUPAC_MASK;
        buf[i] = (buf[i] << IUPAC_SIZE) | carry;
        carry = nextCarry;
    }
    if (accumulated < length) {
        putIUPAC(complement, accumulated, IUPAC_RC[nucleotide]);
    } else {
        cell carry = ((cell)IUPAC_RC[nucleotide]) << ((length-1) % IUPAC_PER_CELL * IUPAC_SIZE);
        for(int i = bufSize-1; i >= 0; i--) {
            cell nextCarry = (complement[i] & IUPAC_MASK) << ((IUPAC_PER_CELL-1) * IUPAC_SIZE);
            complement[i] = (complement[i] >> IUPAC_SIZE) | carry;
            carry = nextCarry;
        }
    }
    accumulated = accumulated < length ? accumulated+1 : accumulated;
};

void IUPACMotifBuilder::putCompact(unsigned nucleotide) {
    put(COMPACT_TO_IUPAC[nucleotide]);
}

void IUPACMotifBuilder::clear() {
    accumulated = 0;
}

IUPACMotif * IUPACMotifBuilder::build(unsigned size) const {
    if (accumulated < size) {
        // TODO add proper exception
        throw std::exception();
    }
    return new IUPACMotif(buf, size);
}

void IUPACMotifBuilder::fillComplementBuf(unsigned size, cell compBuf[]) const {
    //http://stackoverflow.com/questions/2773890/efficient-bitshifting-an-array-of-int
    // Need to shift array by (length-patternLength)*IUPAC_SIZE bits right
    unsigned shift = (std::min(accumulated, length)-size)*IUPAC_SIZE;
    unsigned cellSize = (sizeof(cell) * 8);
    unsigned bitShift = shift % cellSize;
    unsigned arrayShift = shift / cellSize;
    cell mask = (-1L) ^ ((1L << bitShift) - 1L);
    for (int i = 0;  i < bufSize - arrayShift - 1;  ++i)
    {
        compBuf[i] = (complement[i+arrayShift] >> bitShift) |
                     ((complement[i+arrayShift+1] << (cellSize-bitShift)) & mask);
    }
    compBuf[bufSize - arrayShift - 1] = complement[bufSize - 1] >> bitShift;
}

IUPACMotif * IUPACMotifBuilder::buildComplement(unsigned size) const {
    if (accumulated < size) {
        // TODO add proper exception
        throw std::exception();
    }
    cell compBuf[bufSize];
    fillComplementBuf(size, compBuf);
    return new IUPACMotif(compBuf, size);
}

bool IUPACMotifBuilder::matches(IUPACMotif &pattern, bool rc) const {
    if (accumulated < pattern.getLength()) {
        return false;
    }
    if (rc) {
        cell compBuf[bufSize];
        fillComplementBuf(pattern.getLength(), compBuf);

        bool result = pattern.includes(buf, pattern.getLength()) || pattern.includes(compBuf, pattern.getLength());
        return result;
    }
    return pattern.includes(buf, pattern.getLength());
}

IUPACMotifBuilder::~IUPACMotifBuilder() {
    delete buf;
    delete complement;
};
