#define EIGEN_USE_MKL_ALL
#include "Eigen/PardisoSupport"

#include "Eigen/Sparse"
#include <cuComplex.h>
//#include "cuComplex.h"
#include <iostream>

#include "opcoesDeCompilacao.h"

//// defina como verdadeiro para que seja usado double,
//// caso contr�rio ser� utilizado float 
//#define DOUBLE_MODE false 
//
//// *****************************************************************************************
//// N�o alterar a partir daqui***************************************************************
//// *****************************************************************************************
//
//// defini��o do tipo de ponto flutuante
//#if DOUBLE_MODE
//typedef double float_type;
//#else
//typedef float float_type;
//#endif
//
//// defini��o do tipo de n�mero complexo
//#if DOUBLE_MODE
//typedef cuDoubleComplex complex_type;
//#else
//typedef cuFloatComplex complex_type;
//#endif
//
//// defini��o do constritor complexo
//#if DOUBLE_MODE
//#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
//#else
//#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)
//#endif

struct sistema {
	unsigned int* barrasPV;//[NPV]; // indexador simples (iterador)
	unsigned int nPV;// = 0;
	unsigned int* barrasPQ;//[NPQ]; // indexador simples (iterador)
	unsigned int nPQ;// = 0;
	unsigned int barraVO;// = 0;
	unsigned int nB;
	unsigned int nL;
	float_type baseMVA;

	complex_type* Y;//[NB*NB]; // elementos da matriz admit�ncia (via IDX2F)

	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>* spY;
	complex_type* spYval;
	std::complex<float_type>* spYvalE; // uso do eigen
	int* csrRowPtrY;
	int* csrColIndY;
	int* cooRowIndY;
	int nnzY;

	//limite de inje��o de reativos
	float_type* limQinf; // [NPV] // limite inferior de inje��o de reativos de cada barra PV
	float_type* limQsup; // [NPV] // limite superior de inje��o de reativos de cada barra PV

	float_type* VfixadaPV; // [NPV] // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat
};

struct barra {
	unsigned int* id; // bus identifier number

	float_type* V;//[NB];
	float_type* theta;//[NB];

	float_type* Pliq;//[NB];
	float_type* Qliq;//[NB];
	float_type* Pload;//[NB];
	float_type* Qload;//[NB];
	float_type* Pg;//[NB];
	float_type* Qg;//[NB];

	float_type* Vbase;//[NB];

	float_type* gsh;//[NB];
	float_type* bsh;//[NB];

	float_type* phi;
};

// indexador simples (iterador): No do ramo
struct ramo {
	complex_type* z;//[NL];
	float_type* bsh;//[NL];// = 0;
	complex_type* tap;//[NL];// = 0;

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>* eigen_phi;
	unsigned int* d_csrColIndPhi; // gpu only
	unsigned int* d_csrRowPtrPhi; // gpu only
	float_type* phiVal;               // gpu only
	unsigned int nnzPhi;          // gpu only

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	unsigned int* de;//[NL];
	unsigned int* para;//[NL];
};

struct iterativo {
	float_type* Pcalc;//[NB];
	float_type* Qcalc;//[NB];

	//float_type* J;//[(NB-1+NPQ)*(NB-1+NPQ)]; -> Jlim
	unsigned int iteracao, noMax;

	//float_type* deltaP;
	//float_type* deltaQ;
	//float_type* g; -> gLim

	// limite de inje��o de reativos
	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	//limite de inje��o de reativos
	bool* limQ; // FUTURO?:dois bits: o 1o � 0 se o limite atingido � o superior, 1 se o inferior. O segundo �
				// 0 se o limite de inje��o de reativos n�o foi atingido, 1 caso contr�rio
				// nPV entradas
	unsigned short* barrasPVlim; //[NPV]; // vetor contendo as barras de tipo PV ap�s a an�lise de limite de reativos
	unsigned short nPVlim; // n�mero de barras de tipo PV ap�s a an�lise de limite de reativos
	unsigned short* barrasPQlim; //[NPV]; // vetor contendo as barras de tipo PV ap�s a an�lise de limite de reativos
	unsigned short nPQlim; // n�mero de barras de tipo PV ap�s a an�lise de limite de reativos

	float_type* QliqLim;       // atualiza��o do limite do valor de Qliq para as novas barras de tipo PQ
	bool flgMudancaLimInjReat;
};

enum subMatJ {
	//nda = 0,
	H = 1,
	L = 2,
	M = 3,
	N = 4
};


struct h_sparse {
	std::vector<float_type> spJval;
	/* int* */ std::vector<subMatJ> spJsubmatType;
	/* int* */ std::vector<int> cooColIndSubMatJ;
	/* int* */ std::vector<int> cooRowIndSubMatJ;
	/* int* */ std::vector<int> cooColIndJ;
	/* int* */ std::vector<int> cooRowIndJ;
	std::vector<int> csrRowPtrJ;
	int nnzJ;

	std::vector<int> Hpos;
	std::vector<int> Lpos;
	std::vector<int> Mpos;
	std::vector<int> Npos;
};

enum class metodo {
	esparso,
	//esparsoSimples,
	denso,
	hibridoA,
	hibridoB,
	//denso_LAPACKE,
	// paralelo,
	// paraleloSimples,
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
	extern unsigned short no_max_iter;
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
}

// soluciona o sistema linear com abordagem esparsa. Se bem sucedida retorna 0.
bool spSolve(sistema* sistema, barra* barra, ramo* ramo, iterativo* iterativo, h_sparse* h_sparse) {
		// https://eigen.tuxfamily.org/dox/classEigen_1_1Map_3_01SparseMatrixType_01_4.html
		// https://stackoverflow.com/questions/42481785/how-to-set-sparsematrix-valueptr-sparsematrix-outerindexptr-and-sparsematri#comment72127756_42481785
	
		using namespace Eigen;
		using namespace std;
	
		///CSR format: nonZeroArray, rowIndex, colIndex
		unsigned int aux = sistema->nB - 1 + iterativo->nPQlim;
		Map<SparseMatrix<float_type, Eigen::RowMajor>> spMap(/*rowCount*/ aux, /*colCount*/ aux, /*nonZeroCount*/ h_sparse->cooRowIndJ.size(), /*rowIndex*/ h_sparse->csrRowPtrJ.data(), /*colIndex*/ h_sparse->cooColIndJ.data(), /*nonZeroArray*/ h_sparse->spJval.data(), 0);
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
			B(i - 1) = iterativo->gLim[i - 1];
		}
	
		Matrix<float_type, Dynamic, Dynamic> X(sistema->nB - 1 + iterativo->nPQlim, 1);
		// SparseLU<Eigen::SparseMatrix<float_type>, Eigen::COLAMDOrdering<int>> solver;
		PardisoLU<Eigen::SparseMatrix<float_type>> solver(spJ);


		// // Compute the ordering permutation vector from the structural pattern of A
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
			iterativo->gLim[i - 1] = X(i - 1);
		}

	return 0;
}