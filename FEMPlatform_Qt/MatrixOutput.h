#pragma once

#include <iostream>


bool printdoubleMatrix(const char* fname, int rowsize, int colsize, double** mat);

bool printintMatrix(const char* fname, int rowsize, int colsize, int** mat);

bool printintVector(const char* fname, int size, int* vec);

bool printdoubleVector(const char* fname, int size, double* vec);
