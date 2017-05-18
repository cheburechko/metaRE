#ifndef IUPACMOTIF_H_
#define IUPACMOTIF_H_

#include "Motif.h"
#include <string>
#include <exception>
#include "encodings.h"

class IUPACMotif : public Motif{
private:
    cell * buf;
    unsigned length, bufSize;
    std::string sequence;
    bool degenerate;

    void checkDegenerate();
    IUPACMotif() {};
public:
    IUPACMotif(std::string sequence);
    IUPACMotif(cell buf[], unsigned length);
    IUPACMotif(const IUPACMotif& motif);

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
    virtual bool includes(cell buf[], unsigned length) const;

    virtual bool operator==(const Motif& motif) const;
    virtual bool operator> (const Motif& motif) const;
    virtual bool operator>=(const Motif& motif) const;
    virtual bool operator< (const Motif& motif) const;
    virtual bool operator<=(const Motif& motif) const;

    virtual ~IUPACMotif();
};

class DegenerateToCompactException : public std::exception {
private:
    std::string msg;
public:
    DegenerateToCompactException(std::string motif) {
        msg = "Motif can not be encoded as compact: " + motif;
    };
    virtual const char * what() const noexcept {
        return msg.c_str();
    };
};

#endif
