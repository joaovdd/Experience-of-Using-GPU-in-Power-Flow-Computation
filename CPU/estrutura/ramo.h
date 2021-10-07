#pragma once

#ifndef FLUMEN_GPU

    #include <cuComplex.h>
#else
    #include <cuComplex.h>
#endif

#include "../global.h"
#include "sistema.h"

struct ramoType {
	complex_type* z; 
	float_type* bsh; 
	complex_type* tap; 

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor> eigen_phi;

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	int* de; 
	int* para; 
};

void initBranch(sistemaType &sistema, ramoType &ramo);
void finBranch(ramoType &ramo);

void InitCsrPhi(sistemaType& sistema, ramoType& ramo);