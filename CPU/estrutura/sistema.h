#pragma once

#ifndef FLUMEN_GPU
	//#include "../externo/cuComplex.h"
	#include <cuComplex.h>
	#include "../externo/Eigen/Dense"
	#include "../externo/Eigen/Sparse"
#else
    #include <cuComplex.h>
#endif

#include <complex>

#include "../global.h"

enum subMatJ {
	//nda = 0,
	H = 1,
	L = 2,
	M = 3,
	N = 4
};

struct sistemaType {
	int* barrasPV;//[NPV]; // indexador simples (iterador)
	int nPV;// = 0;
	int* barrasPQ;//[NPQ]; // indexador simples (iterador)
	int nPQ;// = 0;
	int barraVO;// = 0;
	int nB;
	int nL;
	float_type baseMVA;

	complex_type* Y;//[NB*NB]; // elementos da matriz admitância (via IDX2F)
	
	//#ifndef FLUMEN_GPU
		Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>* spY;
	//#endif

	complex_type* spYval; // programa da GPU
	std::complex<float_type>* spYvalE; //usado pela biblioteca eigen, complex_type no ponteiro acima na versão GPU
	int* csrRowPtrY;
	int* csrColIndY;
	int* cooRowIndY;
	int nnzY;

	#ifndef FLUMEN_GPU
			       std::vector<float_type> spJval;
		/* int* */ std::vector<subMatJ> spJsubmatType;
		/* int* */ std::vector<int> cooColIndSubMatJ;
		/* int* */ std::vector<int> cooRowIndSubMatJ;
		/* int* */ std::vector<int> cooColIndJ;
		/* int* */ std::vector<int> cooRowIndJ;
		std::vector<int> csrRowPtrJ;
		int nnzJ;
	#endif

	//limite de injeção de reativos
	float_type* limQinf; // [NPV] // limite inferior de injeção de reativos de cada barra PV
	float_type* limQsup; // [NPV] // limite superior de injeção de reativos de cada barra PV

	float_type* VfixadaPV; // [NPV] // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat
};

void initSistema(sistemaType &sistema);
void finSistema(sistemaType &sistema);