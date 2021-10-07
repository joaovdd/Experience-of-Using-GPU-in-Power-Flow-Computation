#include "spSolver.h"

void spEigenSolverNaive(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	using namespace Eigen;

	Matrix<float_type, Dynamic, Dynamic> dA(sistema->nB - 1 + iterativo->nPQlim, sistema->nB - 1 + iterativo->nPQlim);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			for (int j = 1; j <= sistema->nB - 1 + iterativo->nPQlim; j++) {
				dA(i - 1, j - 1) = iterativo->Jlim[IDX2F(i, j, sistema->nB - 1 + iterativo->nPQlim)];
			}
		}

	SparseMatrix<float_type> sA = dA.sparseView(); 

	Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			B(i - 1) = iterativo->gLim[IDX1F(i)];
		}

	Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);

	SparseLU<Eigen::SparseMatrix<float_type>, Eigen::COLAMDOrdering<int>> solver;

	solver.analyzePattern(sA);

	solver.factorize(sA);

	X = solver.solve(B);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			iterativo->gLim[IDX1F(i)] = X(IDX1F(i));
		}
}

bool spSolve(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) { 

	using namespace Eigen;

	int aux = sistema->nB - 1 + iterativo->nPQlim;
	Map<SparseMatrix<float_type, Eigen::RowMajor>> spMap( aux,  aux,  sistema->cooRowIndJ.size(),  sistema->csrRowPtrJ.data(),  sistema->cooColIndJ.data(),  sistema->spJval.data(), 0);
	SparseMatrix<float_type, Eigen::RowMajor> spJ = spMap.eval();

	if (global::verbose_mode) {
		printf("Jacobiano =\n");
		std::cout << spJ << std::endl;
	}

	Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			B(i - 1) = iterativo->gLim[IDX1F(i)];
		}

	Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);

	PardisoLU<Eigen::SparseMatrix<float_type>> solver(spJ);

	X = solver.solve(B);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			iterativo->gLim[IDX1F(i)] = X(IDX1F(i));
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

bool spSolveMKL(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	MKL_INT n = sistema->nB - 1 + iterativo->nPQlim;
    MKL_INT* ia = sistema->csrRowPtrJ.data();
    MKL_INT* ja = sistema->cooColIndJ.data();
    float_type* a = sistema->spJval.data();

    MKL_INT mtype = 11; 
    
    float_type *b = iterativo->gLim, *x = (float_type*)malloc((sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type)); 
	memset(x, 0, (sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type));

    MKL_INT nrhs = 1; 
    
    void *pt[64];
    
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    
    MKL_INT i; 
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