#ifndef ReadFile_h
#define ReadFile_h

#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <string>
#include <iostream>
#include <fcntl.h>
#include <complex>

using namespace std;

template <typename T>
class ReadFile{
  public:
    typedef vector<vector<vector<vector<T> > > > vector4D;
    typedef vector<vector<vector<T> > > vector3D;
    typedef vector<vector<T> > vector2D;
    typedef vector<T> vector1D;
    typedef complex<T> dataT;

    static vector3D read_volume(string path, int num_slice){
        vector3D volume(num_slice, vector2D(num_slice, vector1D(num_slice)));
        vector<dataT> slice(num_slice*num_slice);
        int bytesRead;

        for (int i = 0; i < num_slice; ++i){
            string slice_path = path+"_slice_" + to_string(i) +".dat";
            bytesRead = BinaryFile<dataT>::read(&slice, slice_path);

            for(int j = 0;j < num_slice; ++j){
                for (int k = 0;k < num_slice; ++k){
                    volume[i][j][k] = sqrt(norm(slice[32*j+k]));
                }
            }
        }
        return volume;
    }

    static vector3D read_reference_volume(string path, int num_slice){
      vector3D volume(num_slice, vector2D(num_slice, vector1D(num_slice)));
      vector<float> buffer(num_slice*num_slice*num_slice);
      int bytesRead;
      bytesRead = BinaryFile<float>::read(&buffer, path);

      for (int i = 0; i < num_slice; ++i){
          for(int j = 0;j < num_slice; ++j){
              for (int k = 0;k < num_slice; ++k){
                  volume[i][j][k] = buffer[32*32*i+32*j+k];
              }
          }
      }
      return volume;
  }
};

#endif
