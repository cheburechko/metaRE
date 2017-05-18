#include "IUPACMotif.h"
#include <boost/algorithm/string.hpp>

IUPACMotif::IUPACMotif(std::string sequence) :
    length(sequence.size()),
    sequence(sequence),
    bufSize((length - 1) / IUPAC_PER_CELL + 1)
{
    boost::algorithm::to_lower(this -> sequence);
    buf = new cell[bufSize]();
    for (unsigned i = 0; i < length; i++) {
        putIUPAC(buf, i, CHAR_TO_IUPAC(sequence[sequence.length()-1-i]));
    }
    checkDegenerate();
}

void IUPACMotif::checkDegenerate() {
    //https://en.wikipedia.org/wiki/Hamming_weight
    const uint64_t m1  = 0x5555555555555555; //binary: 0101...
    const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
    degenerate = false;
    cell mask = 0x1111111111111111;

    for (unsigned i = 0; i < bufSize && !degenerate; i++) {
        cell x = buf[i];
        x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
        x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits

        if (i == bufSize-1 && length % IUPAC_PER_CELL != 0) {
            cell maskMask = ((cell)1 << ((length-i*IUPAC_PER_CELL) * IUPAC_SIZE)) - (cell)1;
            mask &= maskMask;
        }

        degenerate |= x != mask;
    }
}


IUPACMotif::IUPACMotif(cell buf[], unsigned length) :
    length(length),
    bufSize((length - 1) / IUPAC_PER_CELL + 1),
    sequence(length, 'a')
{
    this -> buf = (cell*)memcpy(new cell[bufSize](), buf, bufSize*sizeof(cell));
    for (unsigned i = 0; i < length; i++) {
        sequence[sequence.length()-1-i] = IUPAC_TO_CHAR[getIUPAC(buf, i)];
    }
    checkDegenerate();
}

IUPACMotif::IUPACMotif(const IUPACMotif& motif) :
    length(motif.length),
    bufSize(motif.bufSize),
    sequence(motif.sequence)
{
    buf = new cell[bufSize];
    std::memcpy(buf, motif.buf, bufSize*sizeof(cell));
    checkDegenerate();
}

void IUPACMotif::encodeIUPAC(cell buf[]) const {
    std::memcpy(buf, this->buf, bufSize*sizeof(cell));
}

cell * IUPACMotif::encodeIUPAC() const {
    return (cell *) std::memcpy(new cell[bufSize](), this->buf, bufSize*sizeof(cell));
}

void IUPACMotif::encodeCompact(cell buf[]) const {
    if (isDegenerate()) {
        throw DegenerateToCompactException(sequence);
    }
    for (unsigned i = 0; i < length; i++) {
        putCompact(buf, i, IUPAC_TO_COMPACT[getIUPAC(this->buf, i)]);
    }
}
cell * IUPACMotif::encodeCompact() const {
    cell * result = new cell[(length - 1) / COMPACTS_PER_CELL + 1]();
    encodeCompact(result);
    return result;
}

unsigned IUPACMotif::getLength() const {
    return length;
}

std::string IUPACMotif::getString() const {
    return sequence;
}

Motif * IUPACMotif::getReverseComplement() const {
    auto motif = new IUPACMotif();

    motif->length = length;
    motif->buf = new cell[bufSize]();
    motif->bufSize = bufSize;
    motif->sequence = std::string(length, 'a');

    for (unsigned i = 0; i < length; i++) {
        base val = IUPAC_RC[getIUPAC(buf, i)];
        unsigned pos = length-i-1;
        putIUPAC(motif->buf, pos, val);
        motif->sequence[i] = IUPAC_TO_CHAR[val];
    }
    motif->checkDegenerate();
    return motif;
};


bool IUPACMotif::isReverseComplementOf(const Motif& motif) const {
    return (*this) == (*std::unique_ptr<Motif>(motif.getReverseComplement()));
}

bool IUPACMotif::isDegenerate() const {
    return degenerate;
}
bool IUPACMotif::includes(const Motif& motif) const {
    cell * otherBuf = motif.encodeIUPAC();
    bool result = includes(otherBuf, motif.getLength());
    delete otherBuf;
    return result;
}

bool IUPACMotif::includes(cell otherBuf[], unsigned length) const {
    if (this->length != length) {
        return false;
    }

    bool result = true;
    cell mask = 0x1111111111111111;
    cell presentMask = 0xffffffffffffffff;
    for (unsigned i = 0; i < bufSize && result; i++) {
        if (i == bufSize-1 && length % IUPAC_PER_CELL != 0) {
            presentMask = (1L << ((length % IUPAC_PER_CELL) * IUPAC_SIZE)) - 1L;
            mask &= presentMask;
        }

        bool cond1 = (otherBuf[i] & ~buf[i] & presentMask) == 0;
        cell X = otherBuf[i] & buf[i];
        cell Y = (X | X >> 1 | X >> 2 | X >> 3);

        bool cond2 = (Y & mask) == mask;
        result &= cond1 && cond2;
    }
    return result;
}

bool IUPACMotif::operator==(const Motif& motif) const {
    if (motif.getLength() != length) {
        return false;
    }

    cell buf[bufSize];
    motif.encodeIUPAC(buf);
    return memcmp(buf, this->buf, bufSize*sizeof(cell)) == 0;
}

bool IUPACMotif::operator> (const Motif& motif) const {
    if (length > motif.getLength()) {
        return true;
    } else if (length < motif.getLength()) {
        return false;
    } else if (motif.isDegenerate() && !isDegenerate()) {
        return false;
    } else if (isDegenerate() && !motif.isDegenerate()) {
        return true;
    }
    cell buf[bufSize];
    motif.encodeIUPAC(buf);
    return memcmp(this->buf, buf, bufSize*sizeof(cell)) > 0;
}

bool IUPACMotif::operator>=(const Motif& motif) const {
    if (length > motif.getLength()) {
        return true;
    } else if (length < motif.getLength()) {
        return false;
    } else if (motif.isDegenerate() && !isDegenerate()) {
        return false;
    } else if (isDegenerate() && !motif.isDegenerate()) {
        return true;
    }
    cell buf[bufSize];
    motif.encodeIUPAC(buf);
    return memcmp(this->buf, buf, bufSize*sizeof(cell)) >= 0;
}

bool IUPACMotif::operator< (const Motif& motif) const {
    return !((*this) >= motif);
};

bool IUPACMotif::operator<=(const Motif& motif) const {
    return !((*this) > motif);
};

IUPACMotif::~IUPACMotif() {
    delete buf;
}
