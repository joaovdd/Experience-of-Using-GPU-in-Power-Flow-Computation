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
	int aux = sistema->nB - 1 + iterativo->nPQlim;
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MKL_INT int

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

// soluciona o sistema linear com abordagem esparsa. Se bem sucedida retorna 0.
bool spSolveMKL(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	// rowPointers:  h_sparse->csrRowPtrJ.data()
	// colIndex:     h_sparse->cooColIndJ.data()
	// nonZeroArray: h_sparse->spJval.data()

	// b:            iterativo->gLim
	// size:         sistema->nB - 1 + iterativo->nPQlim

	// if (global::verbose_mode) {
	// 	printf("Jacobiano =\n");
	// 	std::cout << spJ << std::endl;
	// }

    /* Matrix data. */
    // MKL_INT n = 5;
    // MKL_INT ia[6] = {0, 3, 5, 8, 11, 13};
    // MKL_INT ja[13] = {0, 1, 3, 0, 1, 2, 3, 4, 0, 2, 3, 1, 4};
    // double a[13] = {1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 7.0, 8.0, -5.0};

	MKL_INT n = sistema->nB - 1 + iterativo->nPQlim;
    MKL_INT* ia = sistema->csrRowPtrJ.data();
    MKL_INT* ja = sistema->cooColIndJ.data();
    float_type* a = sistema->spJval.data();


    MKL_INT mtype = 11; /* Real unsymmetric matrix */
    // Descriptor of main sparse matrix properties
    // struct matrix_descr descrA;
    // Structure with sparse matrix stored in CSR format
    // sparse_matrix_t csrA; // dd [&& csrA]
    // sparse_operation_t transA; // dd [&& res]

    /* RHS and solution vectors. */
    // double b[5], x[5], bs[5], res, res0;
    // double *b = iterativo->gLim, *x = (double*)malloc(n * sizeof(double)), *bs = (double*)malloc(n * sizeof(double)); // dd
    float_type *b = iterativo->gLim, *x = (float_type*)malloc((sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type)); // , *bs = (float_type*)malloc(n * sizeof(float_type)); // dd
	memset(x, 0, (sistema->nPV + sistema->nPV + sistema->nPQ + sistema->nPQ) * sizeof(float_type));
	// float_type res, res0;

    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i; // , j;
    float_type ddum;  /* Double dummy */
    MKL_INT idum; /* Integer dummy. */
                  /* -------------------------------------------------------------------- */
                  /* .. Setup Pardiso control parameters. */
                  /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;   /* No solver default */
    iparm[1] = 2;   /* Fill-in reordering from METIS */
    iparm[3] = 0;   /* No iterative-direct algorithm */
    iparm[4] = 0;   /* No user fill-in reducing permutation */
    iparm[5] = 0;   /* Write solution into x */
    iparm[6] = 0;   /* Not in use */
    iparm[7] = 2;   /* Max numbers of iterative refinement steps */
    iparm[8] = 0;   /* Not in use */
    iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Conjugate transposed/transpose solve */
    iparm[12] = 1;  /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    maxfct = 1;     /* Maximum number of numerical factorizations. */
    mnum = 1;       /* Which factorization to use. */
    msglvl = 0;     /* Print statistical information  */
    error = 0;      /* Initialize error flag */
    iparm[27] = !DOUBLE_MODE;  /* dd: floatType */
    iparm[34] = 1;  /* dd: zero-based */
                    /* -------------------------------------------------------------------- */
                    /* .. Initialize the internal solver memory pointer. This is only */
                    /* necessary for the FIRST call of the PARDISO solver. */
                    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
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
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    if(global::verbose_mode && !global::laconic_mode) printf("\nFactorization completed ... ");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;

    // descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    // descrA.mode = SPARSE_FILL_MODE_UPPER;
    // descrA.diag = SPARSE_DIAG_NON_UNIT;
    // mkl_sparse_d_create_csr ( &csrA, SPARSE_INDEX_BASE_ONE, n, n, ia, ia+1, ja, a );
    // mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, n, n, ia, ia + 1, ja, a); // dd [&& csrA]

    /* Set right hand side to one. */
    // for (i = 0; i < n; i++)
    // {
    //     b[i] = 1;
    // }
    //  Loop over 3 solving steps: Ax=b, AHx=b and ATx=b
    // for (i = 0; i < 3; i++)
    // {
        i = 0;
        iparm[11] = i; /* Conjugate transposed/transpose solve */
        // if (i == 0)
        //     transA = SPARSE_OPERATION_NON_TRANSPOSE; // dd [&& res]
        // else if (i == 1)
        //     transA = SPARSE_OPERATION_CONJUGATE_TRANSPOSE;
        // else
        //     transA = SPARSE_OPERATION_TRANSPOSE;

        if(global::verbose_mode && !global::laconic_mode) printf("\n\nSolving system with iparm[11] = %d ...\n", (int)iparm[11]);
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
        if (error != 0)
        {
            printf("\nERROR during solution: %d", error);
            exit(3);
        }

        // printf("\nThe solution of the system is: ");
        // for (j = 0; j < n; j++)
        // {
        //     printf("\n x [%d] = % f", j, x[j]);
        // }
        if(global::verbose_mode && !global::laconic_mode) printf("\n");

        // Compute residual [&& res]
        // mkl_sparse_d_mv(transA, 1.0, csrA, descrA, x, 0.0, bs);
        // res = 0.0;
        // res0 = 0.0;
        // for (j = 1; j <= n; j++)
        // {
        //     res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
        //     res0 += b[j - 1] * b[j - 1];
        // }
        // res = sqrt(res) / sqrt(res0);
        // if(global::verbose_mode && !global::laconic_mode) printf("\nRelative residual = %e", res);
        // // Check residual
        // if (res > 1e-10)
        // {
        //     printf("Error: residual is too high!\n");
        //     exit(10 + i);
        // }
    // }
    // mkl_sparse_destroy(csrA); // dd [&& csrA]

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

	// put ans on glim and glim on x to be freed
	iterativo->gLim = x;
	// x = b;

	if(b) free(b);
	// if(bs) free(bs);
    return 0;
}