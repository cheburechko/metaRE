#include "CompactMotifBuilder.h"
#include <iostream>
#include <cstring>

CompactMotifBuilder::CompactMotifBuilder(unsigned size) :
    length(size),
    bufSize((length - 1) / COMPACTS_PER_CELL + 1),
    accumulated(0)
{
    buf = new cell[bufSize]();
    complement = new cell[bufSize]();
    pos = 0;
}

bool CompactMotifBuilder::ready(unsigned size) const {
    return accumulated >= size && size <= length;
}

void CompactMotifBuilder::skip() {
    accumulated = 0;
    pos++;
};

void CompactMotifBuilder::put(unsigned nucleotide) {
    cell carry = nucleotide & COMPACT_MASK;
    for (unsigned i = 0; i < bufSize; i++) {
        cell nextCarry = (buf[i] >> ((COMPACTS_PER_CELL-1) * COMPACT_SIZE)) & COMPACT_MASK;
        buf[i] = (buf[i] << COMPACT_SIZE) | carry;
        carry = nextCarry;
    }
    if (accumulated < length) {
        putCompact(complement, accumulated, COMPACT_RC[nucleotide]);
    } else {
        cell carry = ((cell)COMPACT_RC[nucleotide]) << ((length-1) % COMPACTS_PER_CELL * COMPACT_SIZE);
        for(int i = bufSize-1; i >= 0; i--) {
            cell nextCarry = (complement[i] & COMPACT_MASK) << ((COMPACTS_PER_CELL-1) * COMPACT_SIZE);
            complement[i] = (complement[i] >> COMPACT_SIZE) | carry;
            carry = nextCarry;
        }
    }
    accumulated = accumulated < length ? accumulated+1 : accumulated;
    pos++;
};

void CompactMotifBuilder::putIUPAC(unsigned nucleotide) {
    if (IUPAC_TO_COMPACT[nucleotide] == -1) {
        throw std::exception();
    }
    put(IUPAC_TO_COMPACT[nucleotide]);
}

void CompactMotifBuilder::clear() {
    accumulated = 0;
    pos = 0;
}

CompactMotif * CompactMotifBuilder::build(unsigned size) const {
    if (accumulated < size) {
        // TODO add proper exception
        throw std::exception();
    }
    return new CompactMotif(buf, size);
}

void CompactMotifBuilder::fillComplementBuf(unsigned size, cell compBuf[]) const {
    //http://stackoverflow.com/questions/2773890/efficient-bitshifting-an-array-of-int
    // Need to shift array by (length-patternLength)*COMPACT_SIZE bits right
    unsigned shift = (std::min(accumulated, length)-size)*COMPACT_SIZE;
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

CompactMotif * CompactMotifBuilder::buildComplement(unsigned size) const {
    if (accumulated < size) {
        // TODO add proper exception
        throw std::exception();
    }
    cell compBuf[bufSize];
    fillComplementBuf(size, compBuf);
    return new CompactMotif(compBuf, size);
}

void CompactMotifBuilder::write(cell * buf, unsigned size) const {
    unsigned lastCell = (size-1) / COMPACTS_PER_CELL;
    memcpy(buf, this->buf, (lastCell+1)*sizeof(cell));
    buf[lastCell] &= (1 << (size % COMPACTS_PER_CELL * COMPACT_SIZE)) - 1;
}

CompactMotifBuilder::~CompactMotifBuilder() {
    delete buf;
    delete complement;
};
