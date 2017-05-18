#ifndef COMPACT_MOTIF_H_
#define COMPACT_MOTIF_H_

#include "Motif.h"

class CompactMotif: public Motif {
private:
    cell * buf;
    unsigned length, bufSize;
    std::string sequence;

    CompactMotif() {};

public:
    CompactMotif(std::string sequence);
    CompactMotif(cell buf[], unsigned length);
    CompactMotif(const CompactMotif& motif);

    virtual ~CompactMotif();

    virtual void encodeIUPAC(cell buf[]) const;
    virtual cell * encodeIUPAC() const;

    virtual void encodeCompact(cell buf[]) const;
    virtual cell * encodeCompact() const;

    virtual unsigned getLength() const;
    virtual std::string getString() const;
    virtual Motif * getReverseComplement() const;

    virtual bool isReverseComplementOf(const Motif& motif) const;
    virtual bool isDegenerate() const;
    virtual bool includes(const Motif& motif) const;

    virtual bool operator==(const Motif& motif) const;
    virtual bool operator> (const Motif& motif) const;
    virtual bool operator>=(const Motif& motif) const;
    virtual bool operator< (const Motif& motif) const;
    virtual bool operator<=(const Motif& motif) const;
};

#endif
