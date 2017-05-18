#ifndef MOTIF_H_
#define MOTIF_H_

#include <string>
#include <exception>
#include "encodings.h"

class Motif {
public:
    virtual void encodeIUPAC(cell buf[]) const = 0;
    virtual cell * encodeIUPAC() const = 0;

    virtual void encodeCompact(cell buf[]) const = 0;
    virtual cell * encodeCompact() const = 0;

    virtual unsigned getLength() const = 0;
    virtual std::string getString() const = 0;
    virtual Motif * getReverseComplement() const = 0;

    virtual bool isReverseComplementOf(const Motif& motif) const = 0;
    virtual bool isDegenerate() const = 0;
    virtual bool includes(const Motif& motif) const = 0;

    virtual bool operator==(const Motif& motif) const = 0;
    virtual bool operator> (const Motif& motif) const = 0;
    virtual bool operator>=(const Motif& motif) const = 0;
    virtual bool operator< (const Motif& motif) const = 0;
    virtual bool operator<=(const Motif& motif) const = 0;

    virtual ~Motif() {};
};

#endif
