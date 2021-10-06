#pragma once

#ifndef FLUMEN_GPU
	//#include "../externo/cuComplex.h"
    #include <cuComplex.h>
#else
    #include <cuComplex.h>
#endif

#include "../global.h"
#include "sistema.h"

// indexador simples (iterador): No do ramo
struct ramoType {
	complex_type* z; //[NL];
	float_type* bsh; //[NL];// = 0;
	complex_type* tap; //[NL];// = 0;

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor> eigen_phi;
	//int csrColIndPhi; // GPU only
	//int csrRowPtrPhi;

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	int* de; //[NL];
	int* para; //[NL];
};

void initBranch(sistemaType &sistema, ramoType &ramo);
void finBranch(ramoType &ramo);

void InitCsrPhi(sistemaType& sistema, ramoType& ramo);