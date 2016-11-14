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

    size_t size() {
      return this->numElements; 
    }

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
    const size_t numElements;
};

#endif
