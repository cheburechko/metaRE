#ifndef COMPACTMOTIFBUILDER_H_
#define COMPACTMOTIFBUILDER_H_

#include "CompactMotif.h"

class CompactMotifBuilder {
private:
    cell * buf, * complement;
    unsigned length, bufSize, accumulated, pos;

    void fillComplementBuf(unsigned size, cell target[]) const;
public:
    CompactMotifBuilder(unsigned size);

    /// Shows if motif of given size is ready for output
    bool ready(unsigned size) const;
    bool ready() const { return ready(length); };

    void skip();
    void put(unsigned nucleotide);
    void putIUPAC(unsigned nucleotide);
    void clear();

    CompactMotif * build(unsigned size) const;
    CompactMotif * build() const { return build(length); };

    CompactMotif * buildComplement(unsigned size) const;
    CompactMotif * buildComplement() const { return buildComplement(length); };

    void write(cell * buf, unsigned size) const;
    void write(cell * buf) const { write(buf, length); };

    unsigned getLength() const { return length; };
    unsigned getPos() const { return pos; };

    virtual ~CompactMotifBuilder();
};

#endif /* COMPACTMOTIFBUILDER_H_ */
