#pragma once
#include "cusolverDn.h"
#include "opcoesDeCompilacao_2.h"


class cuSolverStuff_type
{
public:
	cuSolverStuff_type();
	~cuSolverStuff_type();
	// declaração das variáveis**************************************************************

	// metodo denso
	cusolverDnHandle_t DnHandle = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t cuSOLVERstatus;
	cudaError_t cuStatus;

	int bufferSize = 0;
	int* info = NULL;
	float_type* buffer = NULL;
	int* ipiv = NULL; // pivoting sequence
	int h_info = 0;

	// metodo esparso
	cusolverSpHandle_t SpHandle = NULL;
	cusparseHandle_t cusparseHandle = NULL; /* used in residual evaluation */
	//cudaStream_t stream = NULL;
	cusparseMatDescr_t descrJ = NULL;

	//cusolverStatus_t cuSOLVERstatus;
	cusparseStatus_t cuSPARSEstatus;

	int singularity = 0; /* -1 if A is invertible under tol. */

	// 	int n = sistPon.nPQ + sistPon.nPV + sistPon.nPQ;

private:

};

cuSolverStuff_type::cuSolverStuff_type()
{

}

cuSolverStuff_type::~cuSolverStuff_type()
{
}