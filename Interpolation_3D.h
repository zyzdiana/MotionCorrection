#ifndef Interpolation_3D_h
#define Interpolation_3D_h

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <complex>
#include <vector>
#include <string>

using namespace std;

// MAX_DERIVS_SIZE is used to allocate a 2D array to store derivatives
static const int MAX_DERIVS_SIZE = 70*70*70;
static const float PI = 3.141592653589793238463;
static const int X_inv[64][64] ={
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
    {-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
    { 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
    {-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
    {18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
    { 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
    { 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
    {-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
    { 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
};

template <typename T>
class Interpolation_3D{
  public:
    typedef vector<T> vector1D;
    typedef vector<T*> pointer1D; //vector of pointers for storing the coefficients.

    static int clip(int input, int min, int max){
        if(input > max)
            return max;
        else if (input < min)
            return min;
        else
            return input;
    }


    // Function for computing the dot product of the inverse matrix with a vector Y
    static T* X_inv_dot(const vector1D *Y, T coeffs[64]){
        int shape = Y->size();

        for (int i = 0; i < shape; ++i){
            float tmp = 0;
            for (int j=0; j<shape; ++j){
               tmp += X_inv[i][j]*Y->at(j);
            }
            coeffs[i] = tmp;
        }
        return coeffs;
    }

    // Function for computing the dot product of two vectors Y.dot(A)
    static float dot(const vector1D *Y, const T* A){
        int shape = Y->size();
        float result = 0; 
        for (int i = 0; i < shape; ++i){
            result += Y->at(i) * (*(A+i));
        }
        return result;
    }

    static void get_target_Y(vector1D *Y, float z, float y, float x){
        Y->at(0)= 1.;
        Y->at(1)= x;
        Y->at(2)= x*x;
        Y->at(3)= x*x*x;
        Y->at(4)= y;
        Y->at(5)= x*y;
        Y->at(6)= x*x*y;
        Y->at(7)= x*x*x*y;
        Y->at(8)= y*y;
        Y->at(9)= x*y*y;
        Y->at(10) = x*x*y*y;
        Y->at(11) = x*x*x*y*y;
        Y->at(12) = y*y*y;
        Y->at(13) = x*y*y*y;
        Y->at(14) = x*x*y*y*y;
        Y->at(15) = x*x*x*y*y*y;
        

        Y->at(16) = z;
        Y->at(17) = x*z;
        Y->at(18) = x*x*z;
        Y->at(19) = x*x*x*z;
        Y->at(20) = y*z;
        Y->at(21) = x*y*z;
        Y->at(22) = x*x*y*z;
        Y->at(23) = x*x*x*y*z;
        Y->at(24) = y*y*z;
        Y->at(25) = x*y*y*z;
        Y->at(26) = x*x*y*y*z;
        Y->at(27) = x*x*x*y*y*z;
        Y->at(28) = y*y*y*z;
        Y->at(29) = x*y*y*y*z;
        Y->at(30) = x*x*y*y*y*z;
        Y->at(31) = x*x*x*y*y*y*z;

        Y->at(32) = z*z;
        Y->at(33) = x*z*z;
        Y->at(34) = x*x*z*z;
        Y->at(35) = x*x*x*z*z;
        Y->at(36) = y*z*z;
        Y->at(37) = x*y*z*z;
        Y->at(38) = x*x*y*z*z;
        Y->at(39) = x*x*x*y*z*z;
        Y->at(40) = y*y*z*z;
        Y->at(41) = x*y*y*z*z;
        Y->at(42) = x*x*y*y*z*z;
        Y->at(43) = x*x*x*y*y*z*z;
        Y->at(44) = y*y*y*z*z;
        Y->at(45) = x*y*y*y*z*z;
        Y->at(46) = x*x*y*y*y*z*z;
        Y->at(47) = x*x*x*y*y*y*z*z;

        Y->at(48) = z*z*z;
        Y->at(49) = x*z*z*z;
        Y->at(50) = x*x*z*z*z;
        Y->at(51) = x*x*x*z*z*z;
        Y->at(52) = y*z*z*z;
        Y->at(53) = x*y*z*z*z;
        Y->at(54) = x*x*y*z*z*z;
        Y->at(55) = x*x*x*y*z*z*z;
        Y->at(56) = y*y*z*z*z;
        Y->at(57) = x*y*y*z*z*z;
        Y->at(58) = x*x*y*y*z*z*z;
        Y->at(59) = x*x*x*y*y*z*z*z;
        Y->at(60) = y*y*y*z*z*z;
        Y->at(61) = x*y*y*y*z*z*z;
        Y->at(62) = x*x*y*y*y*z*z*z;
        Y->at(63) = x*x*x*y*y*y*z*z*z;
    }

    static void tricubic_derivatives(const vector1D *volume, pointer1D *derivatives){
        // get shape of the volume and derivatives
        int shape = cbrt(volume->size());
        int derives_shape = cbrt(derivatives->size());

        // allocate this big static array to store all the coefficients
        static T arrayCoeffs[MAX_DERIVS_SIZE][64];

        // a vector used to temporary store the polynomials of x, y, and z 
        vector1D Y(64);

        for(int i=-15; i<shape+15; ++i){
            for(int j=-15; j<shape+15; ++j){
                for(int k=-15; k<shape+15; ++k){
                    int x1 = i;
                    int y1 = j;
                    int z1 = k;
                    
                    int x0 = x1 - 1;
                    int x2 = x1 + 1;
                    int x3 = x2 + 1;
                    int y0 = y1 - 1;
                    int y2 = y1 + 1;
                    int y3 = y2 + 1;
                    int z0 = z1 - 1;
                    int z2 = z1 + 1;
                    int z3 = z2 + 1;

                    //Wrap Around
                    x0 = (x0 + shape) % shape;
                    x1 = (x1 + shape) % shape;
                    x2 = (x2 + shape) % shape;
                    x3 = (x3 + shape) % shape;

                    y0 = (y0 + shape) % shape;
                    y1 = (y1 + shape) % shape;
                    y2 = (y2 + shape) % shape;
                    y3 = (y3 + shape) % shape;

                    z0 = (z0 + shape) % shape;
                    z1 = (z1 + shape) % shape;
                    z2 = (z2 + shape) % shape;
                    z3 = (z3 + shape) % shape;                    

                    //values of f(x,y,z) at each corner.
                    Y[0] = volume->at(shape*shape*y1 + shape*x1 + z1);
                    Y[1] = volume->at(shape*shape*y1 + shape*x1 + z2);
                    Y[2] = volume->at(shape*shape*y1 + shape*x2 + z1);
                    Y[3] = volume->at(shape*shape*y1 + shape*x2 + z2);
                    Y[4] = volume->at(shape*shape*y2 + shape*x1 + z1);
                    Y[5] = volume->at(shape*shape*y2 + shape*x1 + z2);
                    Y[6] = volume->at(shape*shape*y2 + shape*x2 + z1);
                    Y[7] = volume->at(shape*shape*y2 + shape*x2 + z2);

                    //values of df/dx
                    Y[8] = ((volume->at(shape*shape*y1 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x1 + z0))/2.);
                    Y[9] = ((volume->at(shape*shape*y1 + shape*x1 + z3)-volume->at(shape*shape*y1 + shape*x1 + z1))/2.);
                    Y[10] = ((volume->at(shape*shape*y1 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x2 + z0))/2.);
                    Y[11] = ((volume->at(shape*shape*y1 + shape*x2 + z3)-volume->at(shape*shape*y1 + shape*x2 + z1))/2.);
                    Y[12] = ((volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y2 + shape*x1 + z0))/2.);
                    Y[13] = ((volume->at(shape*shape*y2 + shape*x1 + z3)-volume->at(shape*shape*y2 + shape*x1 + z1))/2.);
                    Y[14] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x2 + z0))/2.);
                    Y[15] = ((volume->at(shape*shape*y2 + shape*x2 + z3)-volume->at(shape*shape*y2 + shape*x2 + z1))/2.);

                    //values of df/dy
                    Y[16] = ((volume->at(shape*shape*y1 + shape*x2 + z1)-volume->at(shape*shape*y1 + shape*x0 + z1))/2.);
                    Y[17] = ((volume->at(shape*shape*y1 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x0 + z2))/2.);
                    Y[18] = ((volume->at(shape*shape*y1 + shape*x3 + z1)-volume->at(shape*shape*y1 + shape*x1 + z1))/2.);
                    Y[19] = ((volume->at(shape*shape*y1 + shape*x3 + z2)-volume->at(shape*shape*y1 + shape*x1 + z2))/2.);
                    Y[20] = ((volume->at(shape*shape*y2 + shape*x2 + z1)-volume->at(shape*shape*y2 + shape*x0 + z1))/2.);
                    Y[21] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x0 + z2))/2.);
                    Y[22] = ((volume->at(shape*shape*y2 + shape*x3 + z1)-volume->at(shape*shape*y2 + shape*x1 + z1))/2.);
                    Y[23] = ((volume->at(shape*shape*y2 + shape*x3 + z2)-volume->at(shape*shape*y2 + shape*x1 + z2))/2.);

                    //values of df/dz
                    Y[24] = ((volume->at(shape*shape*y2 + shape*x1 + z1)-volume->at(shape*shape*y0 + shape*x1 + z1))/2.);
                    Y[25] = ((volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y0 + shape*x1 + z2))/2.);
                    Y[26] = ((volume->at(shape*shape*y2 + shape*x2 + z1)-volume->at(shape*shape*y0 + shape*x2 + z1))/2.);
                    Y[27] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y0 + shape*x2 + z2))/2.);
                    Y[28] = ((volume->at(shape*shape*y3 + shape*x1 + z1)-volume->at(shape*shape*y1 + shape*x1 + z1))/2.);
                    Y[29] = ((volume->at(shape*shape*y3 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x1 + z2))/2.);
                    Y[30] = ((volume->at(shape*shape*y3 + shape*x2 + z1)-volume->at(shape*shape*y1 + shape*x2 + z1))/2.);
                    Y[31] = ((volume->at(shape*shape*y3 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x2 + z2))/2.);

                    //values of d2f/dxdy
                    Y[32] = ((volume->at(shape*shape*y1 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x0 + z2)-volume->at(shape*shape*y1 + shape*x2 + z0)+volume->at(shape*shape*y1 + shape*x0 + z0))/4.);
                    Y[33] = ((volume->at(shape*shape*y1 + shape*x2 + z3)-volume->at(shape*shape*y1 + shape*x0 + z3)-volume->at(shape*shape*y1 + shape*x2 + z1)+volume->at(shape*shape*y1 + shape*x0 + z1))/4.);
                    Y[34] = ((volume->at(shape*shape*y1 + shape*x3 + z2)-volume->at(shape*shape*y1 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x3 + z0)+volume->at(shape*shape*y1 + shape*x1 + z0))/4.);
                    Y[35] = ((volume->at(shape*shape*y1 + shape*x3 + z3)-volume->at(shape*shape*y1 + shape*x1 + z3)-volume->at(shape*shape*y1 + shape*x3 + z1)+volume->at(shape*shape*y1 + shape*x1 + z1))/4.);
                    Y[36] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x0 + z2)-volume->at(shape*shape*y2 + shape*x2 + z0)+volume->at(shape*shape*y2 + shape*x0 + z0))/4.);
                    Y[37] = ((volume->at(shape*shape*y2 + shape*x2 + z3)-volume->at(shape*shape*y2 + shape*x0 + z3)-volume->at(shape*shape*y2 + shape*x2 + z1)+volume->at(shape*shape*y2 + shape*x0 + z1))/4.);
                    Y[38] = ((volume->at(shape*shape*y2 + shape*x3 + z2)-volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y2 + shape*x3 + z0)+volume->at(shape*shape*y2 + shape*x1 + z0))/4.);
                    Y[39] = ((volume->at(shape*shape*y2 + shape*x3 + z3)-volume->at(shape*shape*y2 + shape*x1 + z3)-volume->at(shape*shape*y2 + shape*x3 + z1)+volume->at(shape*shape*y2 + shape*x1 + z1))/4.);

                    //values of d2f/dxdz
                    Y[40] = ((volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y0 + shape*x1 + z2)-volume->at(shape*shape*y2 + shape*x1 + z0)+volume->at(shape*shape*y0 + shape*x1 + z0))/4.);
                    Y[41] = ((volume->at(shape*shape*y2 + shape*x1 + z3)-volume->at(shape*shape*y0 + shape*x1 + z3)-volume->at(shape*shape*y2 + shape*x1 + z1)+volume->at(shape*shape*y0 + shape*x1 + z1))/4.);
                    Y[42] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y0 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x2 + z0)+volume->at(shape*shape*y0 + shape*x2 + z0))/4.);
                    Y[43] = ((volume->at(shape*shape*y2 + shape*x2 + z3)-volume->at(shape*shape*y0 + shape*x2 + z3)-volume->at(shape*shape*y2 + shape*x2 + z1)+volume->at(shape*shape*y0 + shape*x2 + z1))/4.);
                    Y[44] = ((volume->at(shape*shape*y3 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x1 + z2)-volume->at(shape*shape*y3 + shape*x1 + z0)+volume->at(shape*shape*y1 + shape*x1 + z0))/4.);
                    Y[45] = ((volume->at(shape*shape*y3 + shape*x1 + z3)-volume->at(shape*shape*y1 + shape*x1 + z3)-volume->at(shape*shape*y3 + shape*x1 + z1)+volume->at(shape*shape*y1 + shape*x1 + z1))/4.);
                    Y[46] = ((volume->at(shape*shape*y3 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x2 + z2)-volume->at(shape*shape*y3 + shape*x2 + z0)+volume->at(shape*shape*y1 + shape*x2 + z0))/4.);
                    Y[47] = ((volume->at(shape*shape*y3 + shape*x2 + z3)-volume->at(shape*shape*y1 + shape*x2 + z3)-volume->at(shape*shape*y3 + shape*x2 + z1)+volume->at(shape*shape*y1 + shape*x2 + z1))/4.);

                    //values of d2f/dydz
                    Y[48] = ((volume->at(shape*shape*y2 + shape*x2 + z1)-volume->at(shape*shape*y2 + shape*x0 + z1)-volume->at(shape*shape*y0 + shape*x2 + z1)+volume->at(shape*shape*y0 + shape*x0 + z1))/4.);
                    Y[49] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x0 + z2)-volume->at(shape*shape*y0 + shape*x2 + z2)+volume->at(shape*shape*y0 + shape*x0 + z2))/4.);
                    Y[50] = ((volume->at(shape*shape*y2 + shape*x3 + z1)-volume->at(shape*shape*y2 + shape*x1 + z1)-volume->at(shape*shape*y0 + shape*x3 + z1)+volume->at(shape*shape*y0 + shape*x1 + z1))/4.);
                    Y[51] = ((volume->at(shape*shape*y2 + shape*x3 + z2)-volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y0 + shape*x3 + z2)+volume->at(shape*shape*y0 + shape*x1 + z2))/4.);
                    Y[52] = ((volume->at(shape*shape*y3 + shape*x2 + z1)-volume->at(shape*shape*y3 + shape*x0 + z1)-volume->at(shape*shape*y1 + shape*x2 + z1)+volume->at(shape*shape*y1 + shape*x0 + z1))/4.);
                    Y[53] = ((volume->at(shape*shape*y3 + shape*x2 + z2)-volume->at(shape*shape*y3 + shape*x0 + z2)-volume->at(shape*shape*y1 + shape*x2 + z2)+volume->at(shape*shape*y1 + shape*x0 + z2))/4.);
                    Y[54] = ((volume->at(shape*shape*y3 + shape*x3 + z1)-volume->at(shape*shape*y3 + shape*x1 + z1)-volume->at(shape*shape*y1 + shape*x3 + z1)+volume->at(shape*shape*y1 + shape*x1 + z1))/4.);
                    Y[55] = ((volume->at(shape*shape*y3 + shape*x3 + z2)-volume->at(shape*shape*y3 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x3 + z2)+volume->at(shape*shape*y1 + shape*x1 + z2))/4.);

                    //values of d3f/dxdydz
                    Y[56] = ((volume->at(shape*shape*y2 + shape*x2 + z2)-volume->at(shape*shape*y2 + shape*x0 + z2)-volume->at(shape*shape*y2 + shape*x2 + z0)+volume->at(shape*shape*y2 + shape*x0 + z0))
                              -(volume->at(shape*shape*y0 + shape*x2 + z2)-volume->at(shape*shape*y0 + shape*x0 + z2)-volume->at(shape*shape*y0 + shape*x2 + z0)+volume->at(shape*shape*y0 + shape*x0 + z0)))/8.;
                    Y[57] = ((volume->at(shape*shape*y2 + shape*x2 + z3)-volume->at(shape*shape*y2 + shape*x0 + z3)-volume->at(shape*shape*y2 + shape*x2 + z1)+volume->at(shape*shape*y2 + shape*x0 + z1))
                              -(volume->at(shape*shape*y0 + shape*x2 + z3)-volume->at(shape*shape*y0 + shape*x0 + z3)-volume->at(shape*shape*y0 + shape*x2 + z1)+volume->at(shape*shape*y0 + shape*x0 + z1)))/8.;
                    Y[58] = ((volume->at(shape*shape*y2 + shape*x3 + z2)-volume->at(shape*shape*y2 + shape*x1 + z2)-volume->at(shape*shape*y2 + shape*x3 + z0)+volume->at(shape*shape*y2 + shape*x1 + z0))
                              -(volume->at(shape*shape*y0 + shape*x3 + z2)-volume->at(shape*shape*y0 + shape*x1 + z2)-volume->at(shape*shape*y0 + shape*x3 + z0)+volume->at(shape*shape*y0 + shape*x1 + z0)))/8.;
                    Y[59] = ((volume->at(shape*shape*y2 + shape*x3 + z3)-volume->at(shape*shape*y2 + shape*x1 + z3)-volume->at(shape*shape*y2 + shape*x3 + z1)+volume->at(shape*shape*y2 + shape*x1 + z1))
                              -(volume->at(shape*shape*y0 + shape*x3 + z3)-volume->at(shape*shape*y0 + shape*x1 + z3)-volume->at(shape*shape*y0 + shape*x3 + z1)+volume->at(shape*shape*y0 + shape*x1 + z1)))/8.;

                    Y[60] = ((volume->at(shape*shape*y3 + shape*x2 + z2)-volume->at(shape*shape*y3 + shape*x0 + z2)-volume->at(shape*shape*y3 + shape*x2 + z0)+volume->at(shape*shape*y3 + shape*x0 + z0))
                              -(volume->at(shape*shape*y1 + shape*x2 + z2)-volume->at(shape*shape*y1 + shape*x0 + z2)-volume->at(shape*shape*y1 + shape*x2 + z0)+volume->at(shape*shape*y1 + shape*x0 + z0)))/8.;
                    Y[61] = ((volume->at(shape*shape*y3 + shape*x2 + z3)-volume->at(shape*shape*y3 + shape*x0 + z3)-volume->at(shape*shape*y3 + shape*x2 + z1)+volume->at(shape*shape*y3 + shape*x0 + z1))
                              -(volume->at(shape*shape*y1 + shape*x2 + z3)-volume->at(shape*shape*y1 + shape*x0 + z3)-volume->at(shape*shape*y1 + shape*x2 + z1)+volume->at(shape*shape*y1 + shape*x0 + z1)))/8.;
                    Y[62] = ((volume->at(shape*shape*y3 + shape*x3 + z2)-volume->at(shape*shape*y3 + shape*x1 + z2)-volume->at(shape*shape*y3 + shape*x3 + z0)+volume->at(shape*shape*y3 + shape*x1 + z0))
                              -(volume->at(shape*shape*y1 + shape*x3 + z2)-volume->at(shape*shape*y1 + shape*x1 + z2)-volume->at(shape*shape*y1 + shape*x3 + z0)+volume->at(shape*shape*y1 + shape*x1 + z0)))/8.;
                    Y[63] = ((volume->at(shape*shape*y3 + shape*x3 + z3)-volume->at(shape*shape*y3 + shape*x1 + z3)-volume->at(shape*shape*y3 + shape*x3 + z1)+volume->at(shape*shape*y3 + shape*x1 + z1))
                              -(volume->at(shape*shape*y1 + shape*x3 + z3)-volume->at(shape*shape*y1 + shape*x1 + z3)-volume->at(shape*shape*y1 + shape*x3 + z1)+volume->at(shape*shape*y1 + shape*x1 + z1)))/8.;

                    // compute the index to the derivatives vector from loop indicies i, j and k
                    int idx = derives_shape*derives_shape*(j+15) + derives_shape*(i+15) + (k+15);

                    // store the derivatives
                    derivatives->at(idx) = X_inv_dot(&Y, arrayCoeffs[idx]);
                }
            }
        }
    }

    static float tricubic_interp(const pointer1D *derivatives, const float x, const float y, const float z){
        int derives_shape = cbrt(derivatives->size());
        int x1 = int(floor(x));
        int y1 = int(floor(y));
        int z1 = int(floor(z));

        vector1D target_Y(64);
        get_target_Y(&target_Y, y-y1, x-x1, z-z1);

        return dot(&target_Y, derivatives->at(derives_shape*derives_shape*(y1+15) + derives_shape*(x1+15) + (z1+15)) );
    }

    static float to_radian(float angle){
        return angle*PI/180.0;
    }

    static void rotate_coords_3d(vector1D *result, float x, float y, float z, float theta, float wx, float wy, float wz, float ox, float oy, float oz){
        theta = to_radian(theta);
        // make sure w is a unit vetor:
        float tmp = wx*wx + wy*wy + wz*wz;
        if ((tmp != 1) && (tmp != 0)){
            float norm = sqrt(tmp);
            wx = wx/norm;
            wy = wy/norm;
            wz = wz/norm;
        }

        float s = sin(theta);
        float c = cos(theta);
        x = x - ox;
        y = y - oy;
        z = z - oz;

        result->at(0) = c*x+s*(wy*z-wz*y)+(1-c)*(wx*x+wy*y+wz*z)*wx + ox;
        result->at(1) = c*y+s*(wz*x-wx*z)+(1-c)*(wx*x+wy*y+wz*z)*wy + oy;
        result->at(2) = c*z+s*(wx*y-wy*x)+(1-c)*(wx*x+wy*y+wz*z)*wz + oz;
    }
/*
    static vector1D rotate_volume_trilinear(const vector1D *volume, float theta, float wx, float wy, float wz){
       float shape = float( cbrt(volume->size()) );
       vector1D dest(shape*shape*shape);
        //find center of the volume
        float ox = (shape-1)/2;
        float oy = (shape-1)/2;
        float oz = (shape-1)/2;

        // interpolate every point
        vector1D dest_coords(3);
        for(int i = 0; i < shape; ++i){
            for(int j = 0; j < shape; ++j){
                for(int k = 0; k < shape; ++k){
                    rotate_coords_3d(&dest_coords, i, j, k, theta, wx, wy, wz, ox, oy, oz);
                    dest[shape*shape*j+shape*i+k] = trilinear_interp(volume,dest_coords[0],dest_coords[1],dest_coords[2]);
                }
            }
        }
        return dest;
    }
*/
    static vector1D rotate_volume_tricubic(const pointer1D *derivs, float theta, float wx, float wy, float wz){
        float shape = cbrt(derivs->size()) - 30;
        vector1D dest(shape*shape*shape);

        // int derivs_shape = shape+30;
        // pointer1D derivs(derivs_shape*derivs_shape*derivs_shape);
        // tricubic_derivatives(volume, &derivs);

        //find center of the volume
        float ox = (float(shape)-1)/2;
        float oy = (float(shape)-1)/2;
        float oz = (float(shape)-1)/2;

        // interpolate every point
        vector1D dest_coords(3);
        for(int i = 0; i < shape; ++i){
            for(int j = 0; j < shape; ++j){
                for(int k = 0; k < shape; ++k){
                    rotate_coords_3d(&dest_coords, i, j, k, theta, wx, wy, wz, ox, oy, oz);
                    dest[shape*shape*j+shape*i+k] = tricubic_interp(derivs,dest_coords[0],dest_coords[1],dest_coords[2]);
                }
            }
        }
        return dest;
    }
};
#endif
