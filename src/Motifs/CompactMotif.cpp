#include "CompactMotif.h"

#include <boost/algorithm/string.hpp>

#include <cstring>
#include <memory>
#include <iostream>

CompactMotif::CompactMotif(std::string sequence) :
    length(sequence.size()),
    sequence(sequence),
    bufSize((length - 1) / COMPACTS_PER_CELL + 1)
{
    boost::algorithm::to_lower(this -> sequence);
    buf = new cell[bufSize]();
    for (unsigned i = 0; i < length; i++) {
        putCompact(buf, i, CHAR_TO_COMPACT(sequence[sequence.length()-1-i]));
    }
}

CompactMotif::CompactMotif(cell buf[], unsigned length) :
    length(length),
    bufSize((length - 1) / COMPACTS_PER_CELL + 1),
    sequence(length, 'a')
{
    this -> buf = (cell*)memcpy(new cell[bufSize](), buf, bufSize*sizeof(cell));
    for (unsigned i = 0; i < length; i++) {
        sequence[sequence.length()-1-i] = COMPACT_TO_CHAR[getCompact(buf, i)];
    }
}

CompactMotif::CompactMotif(const CompactMotif& motif) :
    length(motif.length),
    bufSize(motif.bufSize),
    sequence(motif.sequence)
{
    buf = new cell[bufSize];
    std::memcpy(buf, motif.buf, bufSize*sizeof(cell));
}

CompactMotif::~CompactMotif() {
    delete buf;
};

void CompactMotif::encodeIUPAC(cell buf[]) const {
    for (unsigned i = 0; i < length; i++) {
        putIUPAC(buf, i, COMPACT_TO_IUPAC[getCompact(this->buf, i)]);
    }
}

cell * CompactMotif::encodeIUPAC() const {
    cell * result = new cell[(length - 1) / IUPAC_PER_CELL + 1]();
    encodeIUPAC(result);
    return result;
}

void CompactMotif::encodeCompact(cell buf[]) const {
    std::memcpy(buf, this->buf, bufSize*sizeof(cell));
}

cell * CompactMotif::encodeCompact() const {
    return (cell *) std::memcpy(new cell[bufSize](), this->buf, bufSize*sizeof(cell));
}

unsigned CompactMotif::getLength() const {
    return length;
}

std::string CompactMotif::getString() const {
    return sequence;
}

Motif * CompactMotif::getReverseComplement() const {
    auto motif = new CompactMotif();

    motif->length = length;
    motif->buf = new cell[bufSize]();
    motif->bufSize = bufSize;
    motif->sequence = std::string(length, 'a');

    for (unsigned i = 0; i < length; i++) {
        base val = 3-getCompact(buf, i);
        unsigned pos = length-i-1;
        putCompact(motif->buf, pos, val);
        motif->sequence[i] = COMPACT_TO_CHAR[val];
    }
    return motif;
};

bool CompactMotif::isReverseComplementOf(const Motif& motif) const {
    return (*this) == (*std::unique_ptr<Motif>(motif.getReverseComplement()));
}

bool CompactMotif::isDegenerate() const {
    return false;
};

bool CompactMotif::includes(const Motif& motif) const {
    return (*this) == motif;
};

bool CompactMotif::operator==(const Motif& motif) const {
    if (motif.isDegenerate()) {
        return false;
    }
    if (motif.getLength() != length) {
        return false;
    }

    auto buf = std::unique_ptr<cell>(motif.encodeCompact());
    return memcmp(buf.get(), this->buf, bufSize*sizeof(cell)) == 0;
};

bool CompactMotif::operator> (const Motif& motif) const {
    if (length > motif.getLength()) {
        return true;
    } else if (length < motif.getLength()) {
        return false;
    } else if (motif.isDegenerate()) {
        return false;
    }
    cell buf[bufSize];
    motif.encodeCompact(buf);
    return memcmp(this->buf, buf, bufSize*sizeof(cell)) > 0;
};
bool CompactMotif::operator>=(const Motif& motif) const {
    if (length > motif.getLength()) {
        return true;
    } else if (length < motif.getLength()) {
        return false;
    } else if (motif.isDegenerate()) {
        return false;
    }
    cell buf[bufSize];
    motif.encodeCompact(buf);
    return memcmp(this->buf, buf, bufSize*sizeof(cell)) >= 0;
};
bool CompactMotif::operator< (const Motif& motif) const {
    return !((*this) >= motif);
};
bool CompactMotif::operator<=(const Motif& motif) const {
    return !((*this) > motif);
};
