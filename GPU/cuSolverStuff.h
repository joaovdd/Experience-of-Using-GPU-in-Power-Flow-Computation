#pragma once
#include "cusolverDn.h"
#include "opcoesDeCompilacao_2.h"

class cuSolverStuff_type
{
public:
	cuSolverStuff_type();
	~cuSolverStuff_type();

	cusolverDnHandle_t DnHandle = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t cuSOLVERstatus;
	cudaError_t cuStatus;

	int bufferSize = 0;
	int* info = NULL;
	float_type* buffer = NULL;
	int* ipiv = NULL; 
	int h_info = 0;

	cusolverSpHandle_t SpHandle = NULL;
	cusparseHandle_t cusparseHandle = NULL; 

	cusparseMatDescr_t descrJ = NULL;

	cusparseStatus_t cuSPARSEstatus;

	int singularity = 0; 

private:

};

cuSolverStuff_type::cuSolverStuff_type()
{
}

cuSolverStuff_type::~cuSolverStuff_type()
{
}