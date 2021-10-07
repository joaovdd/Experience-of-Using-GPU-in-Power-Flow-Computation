#pragma once

#ifndef FLUMEN_GPU

	#include <cuComplex.h>
	#include "../externo/Eigen/Dense"
	#include "../externo/Eigen/Sparse"
#else
    #include <cuComplex.h>
#endif

#include <complex>

#include "../global.h"

enum subMatJ {
	H = 1,
	L = 2,
	M = 3,
	N = 4
};

struct sistemaType {
	int* barrasPV;
	int nPV;
	int* barrasPQ;
	int nPQ;
	int barraVO;
	int nB;
	int nL;
	float_type baseMVA;

	complex_type* Y;

		Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>* spY;

	complex_type* spYval; 
	std::complex<float_type>* spYvalE; 
	int* csrRowPtrY;
	int* csrColIndY;
	int* cooRowIndY;
	int nnzY;

	#ifndef FLUMEN_GPU
			       std::vector<float_type> spJval;
		 std::vector<subMatJ> spJsubmatType;
		 std::vector<int> cooColIndSubMatJ;
		 std::vector<int> cooRowIndSubMatJ;
		 std::vector<int> cooColIndJ;
		 std::vector<int> cooRowIndJ;
		std::vector<int> csrRowPtrJ;
		int nnzJ;
	#endif

	float_type* limQinf; 
	float_type* limQsup; 

	float_type* VfixadaPV; 
};

void initSistema(sistemaType &sistema);
void finSistema(sistemaType &sistema);