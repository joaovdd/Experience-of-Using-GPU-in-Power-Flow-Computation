#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

#ifndef FLUMEN_GPU
    // #include "externo/cuComplex.h"
    #include <cuComplex.h>
#else
    #include <cuComplex.h>
#endif


#include "estruturas.h"
#include "leitura/config.h"

#include "global.h"

void showVec(bool* show, int dim, int prec);

void showVec(float_type* show, int dim, int prec);

void showVec(int* show, int dim, int prec);

void showVec(int* show, int dim, int prec);

void showVecR(complex_type* show, int dim, int prec);

void showVecI(complex_type* show, int dim, int prec);

void showVec(complex_type* show, int dim, int prec);

void showMatI(complex_type* show, int dim, int prec);

void showMatI(complex_type* show, int dim, int prec);

void showMat(const complex_type* show, const int dim);

void showMat(const float_type* show, const int dim);

void showMat(const float* show, const int dim);

void showMat(const bool* show, const int dim);

void showMatRI(const complex_type* show, const int dim);

void printAll(sistemaType& sistema, barraType &barra, ramoType &ramo);