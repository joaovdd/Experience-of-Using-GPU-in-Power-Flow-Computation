#pragma once
#define EIGEN_USE_MKL_ALL

#include "../global.h"
#include "../idx.h"
#include "../estruturas.h"
#ifndef FLUMEN_GPU
    #include "../externo/Eigen/Dense"
#endif

void dnEigenSolver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);