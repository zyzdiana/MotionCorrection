#ifndef BinaryFile_h
#define BinaryFile_h

#include <vector>
#include <string>
#include <iostream>
#include <fcntl.h>


template <typename T>
class BinaryFile{
  public:
    static int read(std::vector<T> *buffer, std::string filePath) {
        int inputFile = open(filePath.c_str(), O_RDONLY);
    
        if(-1 == inputFile) {
          std::cerr << "Could not open file: " << filePath << std::endl;
          exit(1);
        }

        int bytesRead =
            ::read(inputFile, &(buffer->at(0)), sizeof(T) * buffer->size());

        close(inputFile);

        return bytesRead;
    }
    
    static int write(std::vector<T> *buffer, std::string filePath) {
        // create the file with user-only permissions
        int outputFile = open(filePath.c_str(), O_WRONLY | O_CREAT, 0600);
    
        if(-1 == outputFile) {
          std::cerr << "Could not open file: " << filePath << std::endl;
          exit(1);
        }

        int bytesWritten = 
            ::write(outputFile, &(buffer->at(0)), sizeof(T) * buffer->size());
    
        close(outputFile);

        return bytesWritten;
    }
};

#endif
