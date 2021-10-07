#define EIGEN_USE_MKL_ALL
#include "Eigen/PardisoSupport"

#include "Eigen/Sparse"
#include <cuComplex.h>

#include <iostream>

#include "opcoesDeCompilacao.h"

struct sistema {
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

	float_type* limQinf; 
	float_type* limQsup; 

	float_type* VfixadaPV; 
};

struct barra {
	int* id; 

	float_type* V;
	float_type* theta;

	float_type* Pliq;
	float_type* Qliq;
	float_type* Pload;
	float_type* Qload;
	float_type* Pg;
	float_type* Qg;

	float_type* Vbase;

	float_type* gsh;
	float_type* bsh;

	float_type* phi;
};

struct ramo {
	complex_type* z;
	float_type* bsh;
	complex_type* tap;

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>* eigen_phi;
	int* d_csrColIndPhi; 
	int* d_csrRowPtrPhi; 
	float_type* phiVal;               
	int nnzPhi;          

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	int* de;
	int* para;
};

struct iterativo {
	float_type* Pcalc;
	float_type* Qcalc;

	int iteracao, noMax;

	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	bool* limQ; 
				
	int* barrasPVlim; 
	int nPVlim; 
	int* barrasPQlim; 
	int nPQlim; 

	float_type* QliqLim;       
	bool flgMudancaLimInjReat;
};

enum subMatJ {
	H = 1,
	L = 2,
	M = 3,
	N = 4
};

struct h_sparse {
	std::vector<float_type> spJval;
	 std::vector<subMatJ> spJsubmatType;
	 std::vector<int> cooColIndSubMatJ;
	 std::vector<int> cooRowIndSubMatJ;
	 std::vector<int> cooColIndJ;
	 std::vector<int> cooRowIndJ;
	std::vector<int> csrRowPtrJ;
	int nnzJ;

	std::vector<int> Hpos;
	std::vector<int> Lpos;
	std::vector<int> Mpos;
	std::vector<int> Npos;
};

enum class metodo {
	esparso,

	denso,
	hibridoA,
	hibridoB,

	nda
};

enum class metodoDeDecomposicao {
	LU,
	QR,
	nda
};

enum class metodoDeCalculoDeYbus {
	dnCPU,
	spCPU,
	nda
};

enum class output_benchmarkType {
	all,
	screen,
	file
};

namespace global {
	extern float_type v_inicial, theta_inicial;
	extern int no_max_iter;
	extern float_type tol;
	extern std::string arq_entrada;
	extern bool verbose_mode;
	extern bool lim_inj_reat;
	extern bool laconic_mode;
	extern bool openmp;
	extern metodo metodo;
	extern metodoDeDecomposicao metodoDeDecomposicao;
	extern metodoDeCalculoDeYbus metodoDeCalculoDeYbus;

	extern bool temporizador;
	extern output_benchmarkType output_benchmark;
	extern bool output_ans;
	extern bool output_processo_iterativo;
	extern bool streams;

	extern bool isMatpower;
}

bool spSolve(sistema* sistema, barra* barra, ramo* ramo, iterativo* iterativo, h_sparse* h_sparse) {
		
		using namespace Eigen;
		using namespace std;

		int aux = sistema->nB - 1 + iterativo->nPQlim;
		Map<SparseMatrix<float_type, Eigen::RowMajor>> spMap( aux,  aux,  h_sparse->cooRowIndJ.size(),  h_sparse->csrRowPtrJ.data(),  h_sparse->cooColIndJ.data(),  h_sparse->spJval.data(), 0);
		SparseMatrix<float_type, Eigen::RowMajor> spJ = spMap.eval();
		
		if (global::verbose_mode) {
			printf("Jacobiano =\n");
			std::cout << spJ << std::endl;
		}

		Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);
		
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			B(i - 1) = iterativo->gLim[i - 1];
		}

		Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
		
		PardisoLU<Eigen::SparseMatrix<float_type>> solver(spJ);

		X = solver.solve(B);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			iterativo->gLim[i - 1] = X(i - 1);
		}

	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MKL_INT int

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

bool spSolveMKL(sistema* sistema, barra* barra, ramo* ramo, iterativo* iterativo, h_sparse* h_sparse) {
	MKL_INT n = sistema->nB - 1 + iterativo->nPQlim;
    MKL_INT* ia = h_sparse->csrRowPtrJ.data();
    MKL_INT* ja = h_sparse->cooColIndJ.data();
    float_type* a = h_sparse->spJval.data();

    MKL_INT mtype = 11; 
    
    struct matrix_descr descrA;
    
    float_type *b = iterativo->gLim, *x = (float_type*)malloc((sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type)); 
	memset(x, 0, (sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type));
	float_type res, res0;

    MKL_INT nrhs = 1; 
    
    void *pt[64];
    
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    
    MKL_INT i, j;
    float_type ddum;  
    MKL_INT idum; 
                  
    for (i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;   
    iparm[1] = 2;   
    iparm[3] = 0;   
    iparm[4] = 0;   
    iparm[5] = 0;   
    iparm[6] = 0;   
    iparm[7] = 2;   
    iparm[8] = 0;   
    iparm[9] = 13;  
    iparm[10] = 1;  
    iparm[11] = 0;  
    iparm[12] = 1;  
    iparm[13] = 0;  
    iparm[14] = 0;  
    iparm[15] = 0;  
    iparm[16] = 0;  
    iparm[17] = -1; 
    iparm[18] = -1; 
    iparm[19] = 0;  
    maxfct = 1;     
    mnum = 1;       
    msglvl = 0;     
    error = 0;      
    iparm[27] = !DOUBLE_MODE;  
    iparm[34] = 1;  
                    
    for (i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }
    
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
	if(global::verbose_mode && !global::laconic_mode) {
		printf("\nReordering completed ... ");
		printf("\nNumber of nonzeros in factors = %d", iparm[17]);
		printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	}
    
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    if(global::verbose_mode && !global::laconic_mode) printf("\nFactorization completed ... ");
    
    phase = 33;

    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    
        i = 0;
        iparm[11] = i; 
        
        if(global::verbose_mode && !global::laconic_mode) printf("\n\nSolving system with iparm[11] = %d ...\n", (int)iparm[11]);
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
        if (error != 0)
        {
            printf("\nERROR during solution: %d", error);
            exit(3);
        }

        if(global::verbose_mode && !global::laconic_mode) printf("\n");

    phase = -1; 
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

	iterativo->gLim = x;

	if(b) free(b);

    return 0;
}