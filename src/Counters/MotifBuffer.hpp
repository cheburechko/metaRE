#ifndef COUNTERS_MOTIFBUFFER_H_
#define COUNTERS_MOTIFBUFFER_H_

#include <boost/circular_buffer.hpp>
#include "../Motifs/encodings.h"

template<typename T>
class MotifBuffer {
private:
    int pos;
    const int center, bufferSize;

    boost::circular_buffer<unsigned> buffer;
    boost::circular_buffer<bool> validity;

public:
    MotifBuffer(unsigned size);
    unsigned size() const { return bufferSize; };
    void reset();

    T put(T motif);
    T skip();
    void invalidate(unsigned pos);

    bool centerAvailable() const;
    T getCenter() const;
    unsigned getCenterPos() const;

    bool isValid(unsigned pos) const;// { return true; };
    T get(unsigned pos) const;

    int posFromCenter(unsigned pos) const;
    int absPosition(unsigned pos) const;
};

template<typename T> MotifBuffer<T>::MotifBuffer(unsigned size) :
    buffer(size),
    validity(size),
    bufferSize(size),
    center(size/2)
{
    reset();
}

template<typename T>
void MotifBuffer<T>::reset() {
    buffer.clear();
    validity.clear();
    for (int i = 0; i < bufferSize; i++) {
        skip();
    }
    pos = -bufferSize;
}

template<typename T>
T MotifBuffer<T>::put(T motif) {
    pos++;
    T value = buffer[0];
    buffer.push_back(motif);
    validity.push_back(true);
    return value;
}

template<typename T>
T MotifBuffer<T>::skip() {
    pos++;
    T value = buffer[0];
    buffer.push_back(T());
    validity.push_back(false);
    return value;
}

template<typename T>
void MotifBuffer<T>::invalidate(unsigned pos) {
    validity[pos] = false;
}

template<typename T>
bool MotifBuffer<T>::centerAvailable() const {
    return validity[center];
}

template<typename T>
T MotifBuffer<T>::getCenter() const {
    return buffer[center];
}

template<typename T>
unsigned MotifBuffer<T>::getCenterPos() const {
    return center;
}


template<typename T>
bool MotifBuffer<T>::isValid(unsigned pos) const {
    return validity[pos];
}

template<typename T>
T MotifBuffer<T>::get(unsigned pos) const {
    return buffer[pos];
}

template<typename T>
int MotifBuffer<T>::posFromCenter(unsigned pos) const {
    return ((int)pos)-center;
}

template<typename T>
int MotifBuffer<T>::absPosition(unsigned pos) const {
    return ((int)pos) + this->pos;
}

#endif /*  COUNTERS_MOTIFBUFFER_H_ */
