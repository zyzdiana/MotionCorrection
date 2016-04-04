#ifndef FFTWBuffer_h
#define FFTWBuffer_h

#include <cstddef>

template <typename T>
class FFTWBuffer {
  public:
    typedef T value_type;

  public:
    FFTWBuffer(size_t numElements);
    ~FFTWBuffer();

    T& at(size_t index) {
      #ifdef DEBUG
      if(index >= numElements) {
        throw std::out_of_range ("FFTWBuffer::at() : index is out of range"); 
      }
      #endif
      return buffer[index]; 
    }
    
    const T& at(size_t index) const {
      #ifdef DEBUG
      if(index >= numElements) {
        throw std::out_of_range ("FFTWBuffer::at() : index is out of range"); 
      }
      #endif
      return buffer[index]; 
    }

    template <typename X>
    friend class FFTOp;

  protected:
    T *buffer;
    #ifdef DEBUG
    const size_t numElements;
    #endif
};

#endif
