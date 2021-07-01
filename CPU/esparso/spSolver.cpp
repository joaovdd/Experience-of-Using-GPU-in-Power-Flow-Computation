#include "spSolver.h"

void spEigenSolverNaive(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	using namespace Eigen;
	// Cria matriz densa do eigen
	Matrix<float_type, Dynamic, Dynamic> dA(sistema->nB - 1 + iterativo->nPQlim, sistema->nB - 1 + iterativo->nPQlim);
	//Matrix<float_type, Dynamic, Dynamic> A(sistema->nB - 1 + sistema->nB - 1, sistema->nB - 1 + sistema->nB - 1);
	
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			for (int j = 1; j <= sistema->nB - 1 + iterativo->nPQlim; j++) {
				dA(i - 1, j - 1) = iterativo->Jlim[IDX2F(i, j, sistema->nB - 1 + iterativo->nPQlim)];
			}
		}

	SparseMatrix<float_type> sA = dA.sparseView(); // converte para esparsa

	//std::cout << "Here is the matrix A:\n" << A << std::endl;

	Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);
	//Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + sistema->nB - 1, 1);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			B(i - 1) = iterativo->gLim[IDX1F(i)];
		}

	//std::cout << "Here is the matrix b:\n" << B << std::endl;

	Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
	//Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + sistema->nB - 1, 1);

	//X = A.fullPivLu().solve(B);

	SparseLU<Eigen::SparseMatrix<float_type>, Eigen::COLAMDOrdering<int>> solver;
	// Compute the ordering permutation vector from the structural pattern of A
	solver.analyzePattern(sA);
	// Compute the numerical factorization 
	solver.factorize(sA);
	//Use the factors to solve the linear system 
	X = solver.solve(B);

	//std::cout << "Here is the matrix x:\n"/* << X */<< std::endl;

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			iterativo->gLim[IDX1F(i)] = X(IDX1F(i));
		}
}

// soluciona o sistema linear com abordagem esparsa. Se bem sucedida retorna 0.
bool /*buildEigenJandSolve*/spSolve(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) { // apenas 0 based
	// https://eigen.tuxfamily.org/dox/classEigen_1_1Map_3_01SparseMatrixType_01_4.html
	// https://stackoverflow.com/questions/42481785/how-to-set-sparsematrix-valueptr-sparsematrix-outerindexptr-and-sparsematri#comment72127756_42481785

	using namespace Eigen;

	//// coo to csr
	//int /*acc = 0,*/ linha = 0;
	//sistema->csrRowPtrJ.clear();
	//sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	//sistema->csrRowPtrJ[0] = 0;
	//sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;
	////linha++; // se é 1 based
	//for (int j = 1; j < sistema->nnzJ; j++) {
	//	if (sistema->cooRowIndJ[j] != linha /*+1...*/) {
	//		sistema->csrRowPtrJ[linha + 1] = j;
	//		linha++;
	//	}
	//}

	///CSR format: nonZeroArray, rowIndex, colIndex
	unsigned int aux = sistema->nB - 1 + iterativo->nPQlim;
	Map<SparseMatrix<float_type, Eigen::RowMajor>> spMap(/*rowCount*/ aux, /*colCount*/ aux, /*nonZeroCount*/ sistema->cooRowIndJ.size(), /*rowIndex*/ sistema->csrRowPtrJ.data(), /*colIndex*/ sistema->cooColIndJ.data(), /*nonZeroArray*/ sistema->spJval.data(), 0);
	SparseMatrix<float_type, Eigen::RowMajor> spJ = spMap.eval();
	//spJ.reserve(sistema->cooRowIndJ.size());

	if (global::verbose_mode) {
		printf("Jacobiano =\n");
		std::cout << spJ << std::endl;
	}

	//Solver
	Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + iterativo->nPQlim, 1);
	//Matrix<float_type, Dynamic, Dynamic> B(sistema->nB - 1 + sistema->nB - 1, 1);

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			B(i - 1) = iterativo->gLim[IDX1F(i)];
		}

	Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
	//SparseLU<Eigen::SparseMatrix<float_type>, Eigen::COLAMDOrdering<int>> solver;
	PardisoLU<Eigen::SparseMatrix<float_type>> solver(spJ);
	
	// Compute the ordering permutation vector from the structural pattern of A
	// solver.analyzePattern(spJ);
	// // Compute the numerical factorization 
	// solver.factorize(spJ);
	// //Use the factors to solve the linear system 
	// if (solver.info()) {
	// 	// erro
	// 	return 1;
	// }
	// else {
	// 	X = solver.solve(B);
	// }
	X = solver.solve(B);

	//std::cout << "Here is the matrix x:\n"/* << X */<< std::endl;

	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i <= sistema->nB - 1 + iterativo->nPQlim; i++) {
			iterativo->gLim[IDX1F(i)] = X(IDX1F(i));
		}

	return 0;
}