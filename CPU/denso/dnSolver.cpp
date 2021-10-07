#include "dnSolver.h"

void dnEigenSolver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	using namespace Eigen;
	    Matrix<float_type, Dynamic, Dynamic> A(sistema->nB - 1 + iterativo->nPQlim, sistema->nB - 1 + iterativo->nPQlim);
		
		#pragma omp parallel for if (global::openmp)
			for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++){
				for (int j = 1; j <= sistema->nB - 1 + iterativo->nPQlim; j++) {
					A((int)i-1,(int)j-1) = iterativo->Jlim[IDX2F(i, j, sistema->nB - 1 + iterativo->nPQlim)];
				}
			}

	    Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);
		
		#pragma omp parallel for if (global::openmp)
			for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
				B((int)i - 1) = iterativo->gLim[IDX1F(i)];
			}

		Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
		
		X = A.partialPivLu().solve(B);

		#pragma omp parallel for if (global::openmp)
			for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
				iterativo->gLim[IDX1F(i)] = X(IDX1F((int)i));
			}
}