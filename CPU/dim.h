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

void showVec(bool* show, unsigned int dim, unsigned int prec);

void showVec(float_type* show, unsigned int dim, unsigned int prec);

void showVec(unsigned int* show, unsigned int dim, unsigned int prec);

void showVec(int* show, int dim, unsigned int prec);

void showVecR(complex_type* show, unsigned int dim, unsigned int prec);

void showVecI(complex_type* show, unsigned int dim, unsigned int prec);

void showVec(complex_type* show, unsigned int dim, unsigned int prec);

void showMatI(complex_type* show, unsigned int dim, unsigned int prec);

void showMatI(complex_type* show, unsigned int dim, unsigned int prec);

void showMat(const complex_type* show, const unsigned int dim);

void showMat(const float_type* show, const unsigned int dim);

void showMat(const float* show, const unsigned int dim);

void showMat(const bool* show, const unsigned int dim);

void showMatRI(const complex_type* show, const unsigned int dim);

void printAll(sistemaType& sistema, barraType &barra, ramoType &ramo);