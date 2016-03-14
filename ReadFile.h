#ifndef ReadFile_h
#define ReadFile_h

#include "BinaryFile.h"
//#include "Interpolation_3D.h"
//#include "Volume.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fcntl.h>
#include <complex>
#include <fstream>
#include <stdio.h>

using namespace std;

template <typename T>
class ReadFile{
  public:
    typedef vector<T> vector1D;
    typedef complex<T> dataT;

    static void read_volume(vector1D *volume, string path, int num_slice){
        vector<dataT> slice(num_slice*num_slice);
        int bytesRead;

        for (int i = 0; i < num_slice; ++i){
            std::stringstream ss;
            ss << i;
            string slice_path = path+"_slice_" + ss.str() +".dat";
            bytesRead = BinaryFile<dataT>::read(&slice, slice_path);
            for(int j = 0;j < num_slice; ++j){
                for (int k = 0;k < num_slice; ++k){
                    volume->at(32*32*i+32*j+k) = sqrt(norm(slice[32*j+k]));
                }
            }
        }
    }

    static vector1D read_reference_volume(string path, int size){
      vector1D volume(size);
      int bytesRead = BinaryFile<T>::read(&volume, path);
      return volume;
  }

};
#endif
