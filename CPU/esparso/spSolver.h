#pragma once
#include <iostream>

#define EIGEN_USE_MKL_ALL

#include "../externo/Eigen/PardisoSupport"

#include "../idx.h"
#include "../estruturas.h"
#ifndef FLUMEN_GPU
    #include "../externo/Eigen/Sparse"
#endif
#include "../global.h"

void spEigenSolverNaive(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

bool spSolve(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);