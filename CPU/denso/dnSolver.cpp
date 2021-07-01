#include "dnSolver.h"

//inline void LAPACKE_solver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
//	int n = sistema->nPQ + sistema->nPV + sistema->nPQ, nrhs = 1, lda = n, ldb = 1, info;
//	int* ipiv;
//	ipiv = (int*)malloc(n);
//	//              (n_rows, n_cols, matriz, lda, ipiv, info)
//	info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, iterativo->J, n, ipiv); // http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
//
//	if (0 != info) {
//		fprintf(stderr, "Error: LU factorization failed\n");
//	}
//
//	//                   (forma-do-sistema, ordem, colunas na atriz B, A, lda, ipiv, B, ldb, info) //n_rows, n_cols, matriz, lda, )
//	info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, 1, iterativo->J, n, ipiv, iterativo->g, n); // http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga58e332cb1b8ab770270843221a48296d.html#ga58e332cb1b8ab770270843221a48296d
//
////      /* Check for the exact singularity */
//	if (info > 0) {
//		printf("The diagonal element of the triangular factor of A,\n");
//		printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
//		printf("the solution could not be computed.\n");
//		exit(1);
//	}
//	free(ipiv);
//}

void dnEigenSolver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	using namespace Eigen;
	    Matrix<float_type, Dynamic, Dynamic> A(sistema->nB - 1 + iterativo->nPQlim, sistema->nB - 1 + iterativo->nPQlim);
		//Matrix<float_type, Dynamic, Dynamic> A(sistema->nB - 1 + sistema->nB - 1, sistema->nB - 1 + sistema->nB - 1);

		#pragma omp parallel for if (global::openmp)
			for (unsigned int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++){
				for (unsigned int j = 1; j <= sistema->nB - 1 + iterativo->nPQlim; j++) {
					A((unsigned int)i-1,(unsigned int)j-1) = iterativo->Jlim[IDX2F(i, j, sistema->nB - 1 + iterativo->nPQlim)];
				}
			}

		//std::cout << "Here is the matrix A:\n" << A << std::endl;

	    Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);
		//Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + sistema->nB - 1, 1);

		#pragma omp parallel for if (global::openmp)
			for (unsigned int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
				B((unsigned int)i - 1) = iterativo->gLim[IDX1F(i)];
			}

		//std::cout << "Here is the matrix b:\n" << B << std::endl;

		Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
		//Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + sistema->nB - 1, 1);

		//X = A.fullPivLu().solve(B);
		X = A.partialPivLu().solve(B);

		//std::cout << "Here is the matrix x:\n"/* << X */<< std::endl;

		#pragma omp parallel for if (global::openmp)
			for (unsigned int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
				iterativo->gLim[IDX1F(i)] = X(IDX1F((unsigned int)i));
			}
}