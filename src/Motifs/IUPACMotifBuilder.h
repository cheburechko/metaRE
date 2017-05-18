#ifndef IUPACMOTIFBUILDER_H_
#define IUPACMOTIFBUILDER_H_

#include "IUPACMotif.h"

class IUPACMotifBuilder {
private:
    cell * buf, * complement;
    unsigned length, bufSize, accumulated;

    void fillComplementBuf(unsigned size, cell target[]) const;
public:
    IUPACMotifBuilder(unsigned size);

    /// Shows if motif of given size is ready for output
    bool ready(unsigned size) const;
    bool ready() const { return ready(length); };

    void skip();
    void put(unsigned nucleotide);
    void putCompact(unsigned nucleotide);
    void clear();

    IUPACMotif * build(unsigned size) const;
    IUPACMotif * build() const { return build(length); };

    IUPACMotif * buildComplement(unsigned size) const;
    IUPACMotif * buildComplement() const { return buildComplement(length); };

    bool matches(IUPACMotif &pattern, bool rc = false) const;

    virtual ~IUPACMotifBuilder();
};

#endif
