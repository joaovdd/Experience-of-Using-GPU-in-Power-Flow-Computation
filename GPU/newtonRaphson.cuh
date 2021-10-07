#pragma once

#include <time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <cuda_runtime.h>

#include "cublas_v2.h"
#include "cusolverDn.h"

#include "cusparse.h"
#include "cusolverSp.h"

#include "esparso.h"
#include "limiteInjReat.cuh"

#include "solverCPU.h"

float_type maxi(float_type* vet, int dim) {
	float_type aux = 0;
	for (int i = 0; i < dim; i++) {
		if (abs(vet[i]) > aux) {
			aux = abs(vet[i]);
			
		}
	}
	return aux;
}

__global__ void calcResPLim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 
	bool offx = (idx >= sistPon.barraVO);
	idx += offx;

	if (idx <= sistPon.nB) {
		
		iterPon.gLim[IDX1F(idx - offx)] = barraPon.Pliq[IDX1F(idx)]
			- iterPon.Pcalc[IDX1F(idx)];
	}
}

__global__ void calcResQLim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x; 
	const int offset = sistPon.nB - 1;

	if (idx < iterPon.nPQlim) {
		
		iterPon.gLim[idx + offset] = iterPon.QliqLim[IDX1F(iterPon.barrasPQlim[idx])]
			- iterPon.Qcalc[IDX1F(iterPon.barrasPQlim[idx])];
	}
}

void calcResLimf(const sistema sistPon, const barra barraPon, const ramo ramoPon,
	const iterativo iterPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nB - 1; 
	dim3 dimBlock(3 * deviceprop.warpSize, 1);
	dim3 dimGridP((int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	calcResPLim <<<dimGridP, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

	tamanho = iterPon.nPQlim; 
	dim3 dimGridQ((int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	calcResQLim <<<dimGridQ, dimBlock>>> (sistPon, barraPon, ramoPon, iterPon);
}

__global__ void attOlim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 

	if (idx <= sistPon.nB - 1) {
		
		bool offx = (idx >= sistPon.barraVO);
		idx += offx;

		barraPon.theta[IDX1F(idx)] += iterPon.gLim[IDX1F(idx - offx)]; 
	}
}

__global__ void attVlim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 

	if (idx <= iterPon.nPQlim) {
		
		const int offset = sistPon.nB - 1;

		barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] += iterPon.gLim[IDX1F(
			idx + offset)]; 
	}
}

void attVOlimf(const sistema sistPon, const barra barraPon, const ramo ramoPon,
	const iterativo iterPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nB - 1; 

	dim3 dimBlock(3 * deviceprop.warpSize, 1);
	dim3 dimGridO(
		(int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	attOlim <<<dimGridO, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

	tamanho = iterPon.nPQlim; 
	dim3 dimGridV(
		(int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	attVlim <<<dimGridV, dimBlock>>> (sistPon, barraPon, ramoPon, iterPon);
}

void mCalculePQ(sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, barra& barraPon, iterativo& iterPon, cudaStream_t* streams) {
	switch (global::metodo)
	{
	case metodo::esparso:
	case metodo::hibridoA:
	case metodo::hibridoB:
		d_calculePQ_0based(d_sistema, d_ramo, sistPon, barraPon, iterPon, streams);
		break;
	case metodo::denso:
		d_dnCalculePQ(d_sistema, d_barra, d_ramo, sistPon, barraPon, iterPon);
		break;
	default:
		break;
	}
}

void mCalcJac(sistema& h_sistema, barra& h_barra, ramo& h_ramo, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceProp, h_sparse& h_sparse, d_sparse& sparsePon, cudaStream_t* streams) {
	switch (global::metodo)
	{
	case metodo::esparso:
	case metodo::hibridoA:
	case metodo::hibridoB:
		
	{ BENCHMARK_JACOBIANOSTENCIL_FILL
		SpCalcJac(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceProp, streams);
		
		BENCHMARK_SYNC
	}
		
		if (global::verbose_mode) {
			printf("J^(%d) =\n", iterPon.iteracao);
			float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
			checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
			showVec(aux, sparsePon.nnzJ, 5);
			free(aux);
		}
		break;

	case metodo::denso:
	{ BENCHMARK_JACOBIANO
		checkCudaErrors(cudaMemset(iterPon.Jlim, 0, (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * sizeof(float_type))); 
		calcJacLim(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp);
		BENCHMARK_SYNC
	}
		break;
	default:
		break;
	}
}

void mNR(sistema& h_sistema, barra& h_barra, ramo& h_ramo, iterativo& h_iterativo, sistema* d_sistema,
	barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon,
	barra& barraPon, ramo& ramoPon, iterativo& iterPon,
	const cudaDeviceProp& deviceProp, cudaStream_t* streams) {
	cusolverDnHandle_t DnHandle = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t cuSOLVERstatus;
	cudaError_t cuStatus;

	int n = sistPon.nPQ + sistPon.nPV + sistPon.nPQ;

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

	{ BENCHMARK_INIT_LIB
		
		switch (global::metodo)
		{
		case metodo::esparso:
		
			cuSOLVERstatus = cusolverSpCreate(&SpHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSOLVERstatus);

			cuSPARSEstatus = cusparseCreate(&cusparseHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			checkCudaErrors(cudaStreamCreate(&stream));

			cuSOLVERstatus = cusolverSpSetStream(SpHandle, stream);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuSPARSEstatus = cusparseSetStream(cusparseHandle, stream);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			cuSPARSEstatus = cusparseCreateMatDescr(&descrJ);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			cuSPARSEstatus = cusparseSetMatType(descrJ, CUSPARSE_MATRIX_TYPE_GENERAL); 
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			cuSPARSEstatus = cusparseSetMatIndexBase(descrJ, CUSPARSE_INDEX_BASE_ZERO); 
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			break;
		case metodo::hibridoB:
			if (global::CPUsolverFlg) {
				break;
			}
			else {
				
				cuSPARSEstatus = cusparseCreate(&cusparseHandle);
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				cuSPARSEstatus = cusparseCreateMatDescr(&descrJ);
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				cuSPARSEstatus = cusparseSetMatType(descrJ, CUSPARSE_MATRIX_TYPE_GENERAL); 
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				cuSPARSEstatus = cusparseSetMatIndexBase(descrJ, CUSPARSE_INDEX_BASE_ZERO); 
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
			}
		case metodo::denso:

			cuSOLVERstatus = cusolverDnCreate(&DnHandle);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuStatus = cudaStreamCreate(&stream);
			assert(cudaSuccess == cuStatus);

			cuSOLVERstatus = cusolverDnSetStream(DnHandle, stream);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuSOLVERstatus = _cusolverDnGetrf_bufferSize(DnHandle, n, n, (float_type*)iterPon.Jlim, n, &bufferSize); 
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			break;
		default:
			break;
		}
		BENCHMARK_SYNC
	}

	h_sparse h_sparse; 

	initLim(&h_sistema, &h_sparse, &h_barra, iterPon, &h_iterativo);
	if (global::verbose_mode && global::lim_inj_reat) {
		printAlLim(h_sistema, h_barra, h_ramo, h_iterativo);
	}
	d_sparse sparsePon; 
	{ BENCHMARK_INITGPU
		if (global::metodo == metodo::esparso || global::metodo == metodo::hibridoA || global::metodo == metodo::hibridoB) {
			sparsePon = d_initSparse(h_sistema, h_sparse); 
			SpJacCpyH2D(&h_sistema, &h_sparse, sparsePon); 
		}
		BENCHMARK_SYNC
	}

	{ BENCHMARK_PROCESSOITERATIVO
		for (iterPon.iteracao = 0; iterPon.iteracao < global::no_max_iter; iterPon.iteracao++) {
			if (global::verbose_mode) {
				printf("\n****************************************************************** ITERACAO %d ******************************************************************\n", iterPon.iteracao);
			}
			{ BENCHMARK_CALCPQ
				mCalculePQ(d_sistema, d_barra, d_ramo, sistPon, barraPon, iterPon, streams);
				
				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					
					cudaMemcpy(h_iterativo.Pcalc, iterPon.Pcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					cudaMemcpy(h_iterativo.Qcalc, iterPon.Qcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					break;
				case metodo::denso:
					break;
				default:
					break;
				}
				BENCHMARK_SYNC
			}

			if ((iterPon.iteracao >= 1) && global::lim_inj_reat) {
				
				checkCudaErrors(cudaMemcpy(h_iterativo.gLim, iterPon.gLim, (iterPon.ngLim) * sizeof(float_type), cudaMemcpyDeviceToHost));
			}

			if ((iterPon.iteracao >= 1) && global::lim_inj_reat && maxi(&(h_iterativo.gLim[0]), iterPon.ngLim) < global::tol * 10.) { 

				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					break;
				case metodo::denso:
					
					cudaMemcpy(h_iterativo.Pcalc, iterPon.Pcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					cudaMemcpy(h_iterativo.Qcalc, iterPon.Qcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					
					break;
				default:
					break;
				}

				chkLimQ(&h_sistema, &h_barra, &h_iterativo);   
				geraVetPV(&h_sistema, &h_barra, &h_iterativo); 
				iterPon.nPVlim = h_iterativo.nPVlim;
				geraVetPQ(&h_sistema, &h_barra, &h_iterativo); 
				iterPon.nPQlim = h_iterativo.nPQlim;
				iterPon.ngLim = iterPon.nPQlim + iterPon.nPQlim + iterPon.nPVlim; 
				checkCudaErrors(cudaMemcpy(iterPon.barrasPQlim, h_iterativo.barrasPQlim, iterPon.nPQlim * sizeof(int), cudaMemcpyHostToDevice));
				checkCudaErrors(cudaMemcpy(iterPon.barrasPVlim, h_iterativo.barrasPVlim, iterPon.nPVlim * sizeof(int), cudaMemcpyHostToDevice)); 
				checkCudaErrors(cudaMemcpy(iterPon.QliqLim, h_iterativo.QliqLim, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));

				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					
				{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
					
					h_sparse.spJval.clear();
					h_sparse.spJsubmatType.clear();
					h_sparse.cooColIndSubMatJ.clear();
					h_sparse.cooRowIndSubMatJ.clear();
					h_sparse.cooColIndJ.clear();
					h_sparse.cooRowIndJ.clear();
					
					h_sparse.Hpos.clear();
					h_sparse.Lpos.clear();
					h_sparse.Mpos.clear();
					h_sparse.Npos.clear();

					Jstencil0based(&h_sistema, &h_sparse, &h_barra, &h_iterativo);
					
					d_attSparse(h_sistema, h_sparse, sparsePon); 
					
					BENCHMARK_SYNC
				}
					break;
				case metodo::denso:
					
				{ BENCHMARK_JACOBIANO
					
					cuSOLVERstatus = _cusolverDnGetrf_bufferSize(DnHandle,  iterPon.ngLim,  iterPon.ngLim, (float_type*)iterPon.Jlim,  iterPon.ngLim, &bufferSize);
					assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
					BENCHMARK_SYNC
				}

					break;
				default:
					break;
				}
			}

			checkCudaErrors(cudaDeviceSynchronize()); 

			calcResLimf(sistPon, barraPon, ramoPon, iterPon, deviceProp); 

			checkCudaErrors(cudaMemcpy(h_iterativo.gLim, iterPon.gLim, iterPon.ngLim * sizeof(float_type), cudaMemcpyDeviceToHost));

			if (global::verbose_mode) {
				
				printf("Pcalc^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.Pcalc, h_sistema.nB, 4);
				printf("\nQcalc^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.Qcalc, h_sistema.nB, 4);
				printf("\ndeltaP^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.gLim, h_sistema.nB - 1, 4);
				printf("\ndeltaQ^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.gLim + h_sistema.nB - 1, h_iterativo.nPQlim, 4);

				cout << "\nlimQ:" << endl;
				showVec(h_iterativo.limQ, h_sistema.nPV, 4);
				cout << "\nbarrasPVlim:" << endl;
				showVec(h_iterativo.barrasPVlim, h_iterativo.nPVlim, 4);
				cout << "\nbarrasPQlim:" << endl;
				showVec(h_iterativo.barrasPQlim, h_iterativo.nPQlim, 4);
				cout << "\nlimQinf:" << endl;
				showVec(h_sistema.limQinf, h_sistema.nPV, 4);
				cout << "\nlimQsup:" << endl;
				showVec(h_sistema.limQsup, h_sistema.nPV, 4);
				cout << "\nQliqLim:" << endl;
				showVec(h_iterativo.QliqLim, h_sistema.nB, 4);
			}

			float_type maxG = maxi(h_iterativo.gLim, iterPon.ngLim); 
			if (global::verbose_mode || global::output_processo_iterativo && !global::laconic_mode) {
				printf("\nmax(g) = %e", maxG);
			}
			if (maxG <= global::tol) {
				if (global::verbose_mode) {
					printf("\n\nbreak!\nmax(g) = %e\n\n", maxi(h_iterativo.gLim, iterPon.ngLim));
				}
				return; 
			}
			if (global::verbose_mode) {
				printf("     ...     continua:\n\n");
			}

			mCalcJac(h_sistema, h_barra, h_ramo, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp, h_sparse, sparsePon, streams);

			if (global::verbose_mode) {
				
				cudaMemcpy(h_barra.V, barraPon.V, sizeof(float_type) * (sistPon.nB),
					cudaMemcpyDeviceToHost);
				printf("V^(%d) = \n", iterPon.iteracao);
				showVec(h_barra.V, sistPon.nB, 5);

				cudaMemcpy(h_barra.theta, barraPon.theta, sizeof(float_type) * (sistPon.nB),
					cudaMemcpyDeviceToHost);
				printf("theta^(%d) = \n", iterPon.iteracao);
				showVec(h_barra.theta, sistPon.nB, 5);

				cudaMemcpy(h_iterativo.gLim, iterPon.gLim,
					sizeof(float_type) * (sistPon.nB - 1 + sistPon.nPQ),
					cudaMemcpyDeviceToHost);
				printf("deltaPQ = \n");
				showVec(h_iterativo.gLim, sistPon.nB - 1 + iterPon.nPQlim, 5);
			}

			{ BENCHMARK_SISTEMALINEAR
				switch (global::metodo)
				{
				case metodo::esparso:
					{
						int szJ = (h_sistema.nB - 1 + iterPon.nPQlim);

						switch (global::metodoDeDecomposicao) {
						case metodoDeDecomposicao::LU:
							
							assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
							break;
						case metodoDeDecomposicao::QR:
							
							cuSOLVERstatus = _cusolverSpCsrlsvqr(SpHandle, szJ, sparsePon.nnzJ,
																 descrJ, sparsePon.spJval, sparsePon.csrRowPtrJ, sparsePon.cooColIndJ,
																 iterPon.gLim, global::tol, 0,  iterPon.gLim, &singularity);
							assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
							break;
						}

						if (singularity != -1) {
							printf("\n\n[ERRO] J singular para a tolerancia.\n\n");
							assert(1);
						}
					}
					break;
				
				case metodo::hibridoB:

					if (global::CPUsolverFlg)
					{
						h_sparse.spJval.resize(sparsePon.nnzJ);
						checkCudaErrors(cudaMemcpy(h_sparse.spJval.data(), sparsePon.spJval, sizeof(float_type) * (sparsePon.nnzJ),
							cudaMemcpyDeviceToHost));

						spSolveMKL(&h_sistema, &h_barra, &h_ramo, &h_iterativo, &h_sparse);
						
						checkCudaErrors(cudaMemcpy(iterPon.gLim, h_iterativo.gLim, sizeof(float_type) * (iterPon.ngLim),
							cudaMemcpyHostToDevice));
						break;
					}
					else
					{
						_cusparseCsr2dense(cusparseHandle,
							iterPon.ngLim,
							iterPon.ngLim,
							descrJ,
							sparsePon.spJval,
							sparsePon.csrRowPtrJ,
							sparsePon.cooColIndJ,
							iterPon.Jlim,
							iterPon.ngLim);
					}

				case metodo::denso:

					checkCudaErrors(cudaMalloc(&info, sizeof(int)));
					checkCudaErrors(cudaMalloc(&buffer, sizeof(float_type) * bufferSize));
					checkCudaErrors(cudaMalloc(&ipiv, sizeof(int) *  iterPon.ngLim));

					cudaMemset(info, 0, sizeof(int));

					_cusolverDnGetrf(DnHandle,  iterPon.ngLim,  iterPon.ngLim, iterPon.Jlim,  iterPon.ngLim, buffer, ipiv, info);
					checkCudaErrors(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

					if (0 != h_info) {
						fprintf(stderr, "Error: LU factorization failed\n");
					}

					_cusolverDnGetrs(DnHandle, CUBLAS_OP_N,  iterPon.ngLim, 1, iterPon.Jlim,  iterPon.ngLim, ipiv, iterPon.gLim,
						 iterPon.ngLim, info);
					checkCudaErrors(cudaDeviceSynchronize());

					if (info) { checkCudaErrors(cudaFree(info)); }
					if (buffer) { checkCudaErrors(cudaFree(buffer)); }
					if (ipiv) { checkCudaErrors(cudaFree(ipiv)); }

					break;
				case metodo::nda:
					break;
				default:
					break;
				}
				BENCHMARK_SYNC
			}

			if (global::verbose_mode) {
				checkCudaErrors(cudaMemcpy(h_iterativo.gLim, iterPon.gLim, iterPon.ngLim * sizeof(float_type), cudaMemcpyDeviceToHost));

				printf("\n (deltaX)g^(%d) = \n", iterPon.iteracao);
				checkCudaErrors(cudaDeviceSynchronize());
				showVec(h_iterativo.gLim, iterPon.ngLim, 5);
			}

			attVOlimf(sistPon, barraPon, ramoPon, iterPon, deviceProp); 

			{ BENCHMARK_CUDAMEMCPY
				checkCudaErrors(cudaMemcpy(h_barra.V, barraPon.V, sizeof(float_type)* (sistPon.nB),
					cudaMemcpyDeviceToHost));
				checkCudaErrors(cudaMemcpy(h_barra.theta, barraPon.theta, sizeof(float_type)* (sistPon.nB),
					cudaMemcpyDeviceToHost));
				BENCHMARK_SYNC
			}

			if (global::verbose_mode) {
				
				printf("V^(%d) = \n", iterPon.iteracao + 1);
				showVec(h_barra.V, sistPon.nB, 6);

				printf("theta^(%d) = \n", iterPon.iteracao + 1);
				showVec(h_barra.theta, sistPon.nB, 6);
			}
		}
		BENCHMARK_SYNC
	}
	   
	switch (global::metodo)
	{
	case metodo::esparso:
		d_finSparse(sparsePon); 

		if (SpHandle) {
			cuSOLVERstatus = cusolverSpDestroy(SpHandle);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
		}
		if (cusparseHandle) {
			cuSPARSEstatus = cusparseDestroy(cusparseHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
		}
		if (stream) { checkCudaErrors(cudaStreamDestroy(stream)); }
		if (descrJ) {
			cuSPARSEstatus = cusparseDestroyMatDescr(descrJ);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
		}

		break;
	case metodo::hibridoB:
		d_finSparse(sparsePon);

		if (SpHandle) {
			cuSOLVERstatus = cusolverSpDestroy(SpHandle);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
		}
		if (cusparseHandle) {
			cuSPARSEstatus = cusparseDestroy(cusparseHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
		}
		if (stream) { checkCudaErrors(cudaStreamDestroy(stream)); }
		if (descrJ) {
			cuSPARSEstatus = cusparseDestroyMatDescr(descrJ);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
		}
	case metodo::denso:

		if (DnHandle) {
			cuSOLVERstatus = cusolverDnDestroy(DnHandle);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
		}

		if (stream) {
			cuStatus = cudaStreamDestroy(stream);
			assert(cudaSuccess == cuStatus);
		}

		break;
	case metodo::nda:
		break;
	default:
		break;
	}
}

