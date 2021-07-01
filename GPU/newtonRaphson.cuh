#pragma once
//#include "sistInfo.h"
//#include "sistema.h"
//#include "dim.h"
//#include "PQcalc.cuh"
//#include "Jacobiano.cuh"

//#include "lapacke.h"

#include <time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include <cuda_runtime.h>

#include "cublas_v2.h"
#include "cusolverDn.h"

//#include "helper_cuda.h"
//#include "helper_cusolver.h"

#include "cusparse.h"
#include "cusolverSp.h"

#include "esparso.h"
#include "limiteInjReat.cuh"

#include "solverCPU.h"

float_type maxi(float_type* vet, unsigned int dim) {
	float_type aux = 0;
	for (unsigned int i = 0; i < dim; i++) {
		if (abs(vet[i]) > aux) {
			aux = abs(vet[i]);
			//printf("oi %d\n", aux);
		}
	}
	return aux;
}

// idx = 1:nPQ+nPV
__global__ void calcResPLim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	bool offx = (idx >= sistPon.barraVO);
	idx += offx;

	if (idx <= sistPon.nB) {
		//printf("idx = %d Pliq = %f Pcalc = %f\n", idx, barraPon.Pliq[IDX1F(idx)], iterPon.Pcalc[IDX1F(idx)]);

		iterPon.gLim[IDX1F(idx - offx)] = barraPon.Pliq[IDX1F(idx)]
			- iterPon.Pcalc[IDX1F(idx)];
	}
}

// idx = 1:iterativo->nPQLim
__global__ void calcResQLim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x; //começa de 0
	const unsigned int offset = sistPon.nB - 1;
	//printf("idx = %d sistPon.nPQ = %d ; %d\n", idx, sistPon.nPQ, idx < sistPon.nPQ);

	if (idx < iterPon.nPQlim) {
		//printf("idx = %d\n", idx);

		//printf("idx = %d offset = %d IDX1F(sistPon.barrasPQ[idx]) = %d IDX1F(sistPon.barrasPQ[idx]) = %d\n", idx, offset, IDX1F(sistPon.barrasPQ[idx]), IDX1F(sistPon.barrasPQ[idx]));
		iterPon.gLim[idx + offset] = iterPon.QliqLim[IDX1F(iterPon.barrasPQlim[idx])]
			- iterPon.Qcalc[IDX1F(iterPon.barrasPQlim[idx])];
	}
}

// Calcula o resíduo
void calcResLimf(const sistema sistPon, const barra barraPon, const ramo ramoPon,
	const iterativo iterPon, cudaDeviceProp deviceprop) {
	// deltaPs
	unsigned int tamanho = sistPon.nB - 1; // nPQ+nPV
	dim3 dimBlock(3 * deviceprop.warpSize, 1);
	dim3 dimGridP((unsigned int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	calcResPLim <<<dimGridP, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

	//	cudaDeviceSynchronize();
		//std::cin.get();

		// deltaQs
	tamanho = iterPon.nPQlim; 
	dim3 dimGridQ((unsigned int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	calcResQLim <<<dimGridQ, dimBlock>>> (sistPon, barraPon, ramoPon, iterPon);
}

// idx = 1:nPQ+nPV
__global__ void attOlim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	//printf("attO: idx=%d bloco=%d\n", idx, blockIdx.x);
	if (idx <= sistPon.nB - 1) {
		//printf("CALC_attO: idx=%d bloco=%d\n", idx, blockIdx.x);
		bool offx = (idx >= sistPon.barraVO);
		idx += offx;

		barraPon.theta[IDX1F(idx)] += iterPon.gLim[IDX1F(idx - offx)]; // atualiza valor de theta; g é igual a -deltaX!
	}
}

// idx = 1:nPQ+nPV
__global__ void attVlim(const sistema sistPon, const barra barraPon,
	const ramo ramoPon, const iterativo iterPon) {
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	//printf("attV: idx=%d bloco=%d\n", idx, blockIdx.x);
	if (idx <= iterPon.nPQlim) {
		//printf("CALC_attV: idx=%d bloco=%d\n", idx, blockIdx.x);
		const unsigned int offset = sistPon.nB - 1;

		barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] += iterPon.gLim[IDX1F(
			idx + offset)]; // atualiza valor de V; g é igual a -deltaX!
	}
}

// Calcula o resíduo
void attVOlimf(const sistema sistPon, const barra barraPon, const ramo ramoPon,
	const iterativo iterPon, cudaDeviceProp deviceprop) {
	// theta
	unsigned int tamanho = sistPon.nB - 1; // nPQ+nPV
	//printf("4/3*32) = %f",ceil((float) tamanho/ ((float) (3*deviceprop.warpSize))));
	dim3 dimBlock(3 * deviceprop.warpSize, 1);
	dim3 dimGridO(
		(unsigned int)ceil(
		((float)tamanho) / (float)(3 * deviceprop.warpSize)), 1);

	attOlim <<<dimGridO, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

	// V
	tamanho = iterPon.nPQlim; // nPQ
	dim3 dimGridV(
		(unsigned int)ceil(
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
		//calculo na CPU
		//fillJstencil0based(&h_sistema, &h_sparse, &h_barra, &h_ramo, &h_iterativo);
		//checkCudaErrors(cudaMemcpy(sparsePon.spJval, h_sparse.spJval.data(), sparsePon.nnzJ * sizeof(float_type), cudaMemcpyHostToDevice));
	{ BENCHMARK_JACOBIANOSTENCIL_FILL
		SpCalcJac(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceProp, streams);
		//SpCalcJUnof(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceProp, streams);
		BENCHMARK_SYNC
	}
		// print J
		if (global::verbose_mode) {
			printf("J^(%d) =\n", iterPon.iteracao);
			float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
			checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
			showVec(aux, sparsePon.nnzJ, 5);
			free(aux);
		}
		break;
	//case esparsoSimples:
	//	checkCudaErrors(cudaMemset(iterPon.Jlim, 0, (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * sizeof(float_type))); // inicializa como 0
	//	calcJacLim(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp);
	//	break;
	case metodo::denso:
	{ BENCHMARK_JACOBIANO
		checkCudaErrors(cudaMemset(iterPon.Jlim, 0, (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * sizeof(float_type))); // inicializa como 0
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

	// declaração das variáveis**************************************************************

	// metodo denso
	cusolverDnHandle_t DnHandle = NULL;
	cudaStream_t stream = NULL;

	cusolverStatus_t cuSOLVERstatus;
	cudaError_t cuStatus;

	int n = sistPon.nPQ + sistPon.nPV + sistPon.nPQ;

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

	//***************************************************************************************
	{ BENCHMARK_INIT_LIB
		// Inicialização das variáveis
		switch (global::metodo)
		{
		case metodo::esparso:
		//case esparsoSimples:

			cuSOLVERstatus = cusolverSpCreate(&SpHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSOLVERstatus);

			cuSPARSEstatus = cusparseCreate(&cusparseHandle);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			checkCudaErrors(cudaStreamCreate(&stream));

			/* bind stream to cusparse and cusolver*/
			cuSOLVERstatus = cusolverSpSetStream(SpHandle, stream);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuSPARSEstatus = cusparseSetStream(cusparseHandle, stream);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			/* configure matrix descriptor*/
			cuSPARSEstatus = cusparseCreateMatDescr(&descrJ);
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			cuSPARSEstatus = cusparseSetMatType(descrJ, CUSPARSE_MATRIX_TYPE_GENERAL); // a matriz J não é simetrica, hermitiana, OU triangular
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			/* One based index */
			cuSPARSEstatus = cusparseSetMatIndexBase(descrJ, CUSPARSE_INDEX_BASE_ZERO); // antes CUSPARSE_INDEX_BASE_ONE
			assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

			break;
		case metodo::hibridoB:
			if (global::CPUsolverFlg) {
				break;
			}
			else {
				//cuSOLVERstatus = cusolverSpCreate(&SpHandle);
				//assert(CUSPARSE_STATUS_SUCCESS == cuSOLVERstatus);

				cuSPARSEstatus = cusparseCreate(&cusparseHandle);
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				//checkCudaErrors(cudaStreamCreate(&stream));

				/* bind stream to cusparse and cusolver*/
				//cuSOLVERstatus = cusolverSpSetStream(SpHandle, stream);
				//assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

				//cuSPARSEstatus = cusparseSetStream(cusparseHandle, stream);
				//assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				/* configure matrix descriptor*/
				cuSPARSEstatus = cusparseCreateMatDescr(&descrJ);
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				cuSPARSEstatus = cusparseSetMatType(descrJ, CUSPARSE_MATRIX_TYPE_GENERAL); // a matriz J não é simetrica, hermitiana, OU triangular
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

				/* One based index */
				cuSPARSEstatus = cusparseSetMatIndexBase(descrJ, CUSPARSE_INDEX_BASE_ZERO); // antes CUSPARSE_INDEX_BASE_ONE
				assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
			}
		case metodo::denso:

			cuSOLVERstatus = cusolverDnCreate(&DnHandle);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuStatus = cudaStreamCreate(&stream);
			assert(cudaSuccess == cuStatus);

			cuSOLVERstatus = cusolverDnSetStream(DnHandle, stream);
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			cuSOLVERstatus = _cusolverDnGetrf_bufferSize(DnHandle, n, n, (float_type*)iterPon.Jlim, n, &bufferSize); // necessária adaptação aqui para a rotina com limInjReat: passar criação e testruição para as rotinas de initLim e chk 
			assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);

			break;
		default:
			break;
		}
		BENCHMARK_SYNC
	}
	

	// Init ***********************************************************
	
	h_sparse h_sparse; // Inicialização das estruturas do problema esparso

	// Inicialização das estruturas do limite de injecao de reativos
	initLim(&h_sistema, &h_sparse, &h_barra, iterPon, &h_iterativo);
	if (global::verbose_mode && global::lim_inj_reat) {
		printAlLim(h_sistema, h_barra, h_ramo, h_iterativo);
	}
	d_sparse sparsePon; // Inicialização das estruturas do problema esparso
	{ BENCHMARK_INITGPU
		if (global::metodo == metodo::esparso || global::metodo == metodo::hibridoA || global::metodo == metodo::hibridoB) {
			sparsePon = d_initSparse(h_sistema, h_sparse); // armazena J esparso
			SpJacCpyH2D(&h_sistema, &h_sparse, sparsePon); // copia stencil criado
		}
		BENCHMARK_SYNC
	}


	//Jstencil0based(&h_sistema, &h_sparse, &h_barra, &h_iterativo);

	//*****************************************************************
	{ BENCHMARK_PROCESSOITERATIVO
		for (iterPon.iteracao = 0; iterPon.iteracao < global::no_max_iter; iterPon.iteracao++) {
			if (global::verbose_mode) {
				printf("\n****************************************************************** ITERACAO %d ******************************************************************\n", iterPon.iteracao);
			}
			{ BENCHMARK_CALCPQ
				mCalculePQ(d_sistema, d_barra, d_ramo, sistPon, barraPon, iterPon, streams);
				//checkCudaErrors(cudaDeviceSynchronize()); // debug

				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					// precisa estar fora da região verbose para uso da função host para preencher o stencil
					// e para controle de reativos
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

			//// precisa estar fora da região verbose para uso da função host para preencher o stencil
			//// e para controle de reativos
			//cudaMemcpy(h_iterativo.Pcalc, iterPon.Pcalc, sizeof(float_type) * (sistPon.nB),
			//	cudaMemcpyDeviceToHost);
			//cudaMemcpy(h_iterativo.Qcalc, iterPon.Qcalc, sizeof(float_type) * (sistPon.nB),
			//	cudaMemcpyDeviceToHost);

			if ((iterPon.iteracao >= 1) && global::lim_inj_reat) {
				// necessário para a proposição lógica da condicional de controle (limInjReat)
				checkCudaErrors(cudaMemcpy(h_iterativo.gLim, iterPon.gLim, (iterPon.ngLim) * sizeof(float_type), cudaMemcpyDeviceToHost));
			}

			// controle (limInjReat)
			if ((iterPon.iteracao >= 1) && global::lim_inj_reat && maxi(&(h_iterativo.gLim[0]), iterPon.ngLim) < global::tol * 10.) { // pula a primeira e segunda iterações e deltaX deve ser inferior a um centésimo da tolerancia

				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					break;
				case metodo::denso:
					// precisa estar fora da região verbose para uso da função host para preencher o stencil
					// e para controle de reativos
					cudaMemcpy(h_iterativo.Pcalc, iterPon.Pcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					cudaMemcpy(h_iterativo.Qcalc, iterPon.Qcalc, sizeof(float_type) * (sistPon.nB),
						cudaMemcpyDeviceToHost);
					//checkCudaErrors(cudaDeviceSynchronize()); // debug
					break;
				default:
					break;
				}

				chkLimQ(&h_sistema, &h_barra, &h_iterativo);   // atualiza sistema->limQ
				geraVetPV(&h_sistema, &h_barra, &h_iterativo); // atualiza sistema->barrasPVlim e iterativo->nPQlim
				iterPon.nPVlim = h_iterativo.nPVlim;
				geraVetPQ(&h_sistema, &h_barra, &h_iterativo); // atualiza sistema->barrasPQlim e iterativo->nPQlim  [[ver openmp while]]
				iterPon.nPQlim = h_iterativo.nPQlim;
				iterPon.ngLim = iterPon.nPQlim + iterPon.nPQlim + iterPon.nPVlim; // atualiza o tamanho de g com os limites
				checkCudaErrors(cudaMemcpy(iterPon.barrasPQlim, h_iterativo.barrasPQlim, iterPon.nPQlim * sizeof(unsigned short), cudaMemcpyHostToDevice));
				checkCudaErrors(cudaMemcpy(iterPon.barrasPVlim, h_iterativo.barrasPVlim, iterPon.nPVlim * sizeof(unsigned short), cudaMemcpyHostToDevice)); // ver transferencia device-device
				checkCudaErrors(cudaMemcpy(iterPon.QliqLim, h_iterativo.QliqLim, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));

				switch (global::metodo)
				{
				case metodo::esparso:
				case metodo::hibridoA:
				case metodo::hibridoB:
					// verificar se ocorreu mudanca no chkLimQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
					//zera Jac
					h_sparse.spJval.clear();
					h_sparse.spJsubmatType.clear();
					h_sparse.cooColIndSubMatJ.clear();
					h_sparse.cooRowIndSubMatJ.clear();
					h_sparse.cooColIndJ.clear();
					h_sparse.cooRowIndJ.clear();
					//h_sparse.csrRowPtrJ.clear();

					h_sparse.Hpos.clear();
					h_sparse.Lpos.clear();
					h_sparse.Mpos.clear();
					h_sparse.Npos.clear();

					Jstencil0based(&h_sistema, &h_sparse, &h_barra, &h_iterativo);
					//sparsePon.nnzJ = h_sparse.nnzJ;
					d_attSparse(h_sistema, h_sparse, sparsePon); // ineficiência!!!!
					//SpJacCpyH2D(&h_sistema, &h_sparse, sparsePon); // copia stencil criado // tá errado!!!!
					BENCHMARK_SYNC
				}
					break;
				case metodo::denso:
					// verificar se ocorreu mudanca no chkLimQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
				{ BENCHMARK_JACOBIANO
					// atualiza o buffer
					cuSOLVERstatus = _cusolverDnGetrf_bufferSize(DnHandle, /*n*/ iterPon.ngLim, /*n*/ iterPon.ngLim, (float_type*)iterPon.Jlim, /*n*/ iterPon.ngLim, &bufferSize);
					assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
					BENCHMARK_SYNC
				}

					break;
				default:
					break;
				}
			}

			checkCudaErrors(cudaDeviceSynchronize()); // deve-se terminar o calculo de p e q antes de prosseguir

			calcResLimf(sistPon, barraPon, ramoPon, iterPon, deviceProp); // calcula deltaP's e deltaQ's

			//checkCudaErrors(cudaMemcpy(h_iterativo.g, iterPon.g, (h_sistema.nB - 1 + h_sistema.nPQ) * sizeof(float_type), cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_iterativo.gLim, iterPon.gLim, iterPon.ngLim * sizeof(float_type), cudaMemcpyDeviceToHost));

			if (global::verbose_mode) {
				//checkCudaErrors(cudaMemcpy(h_iterativo.QliqLim, iterPon.QliqLim, sistPon.nB * sizeof(float_type), cudaMemcpyDeviceToHost));

				printf("Pcalc^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.Pcalc, h_sistema.nB, 4);
				printf("\nQcalc^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.Qcalc, h_sistema.nB, 4);
				printf("\ndeltaP^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.gLim, h_sistema.nB - 1, 4);
				printf("\ndeltaQ^(%d) =\n", iterPon.iteracao);
				showVec(h_iterativo.gLim + h_sistema.nB - 1, h_iterativo.nPQlim, 4);

				//printf("\nQg_barrasPV^(%d) =\n", h_iterativo.iteracao);
				//cout << '[' << h_iterativo.Qcalc[IDX1F(h_sistema.barrasPV[0])] + h_barra.Qload[IDX1F(h_sistema.barrasPV[0])] << endl;
				//for (unsigned short i = 1; i < (h_sistema.nPV - 1); i++) {
				//	cout << ' ' << h_iterativo.Qcalc[IDX1F(h_sistema.barrasPV[i])] + h_barra.Qload[IDX1F(h_sistema.barrasPV[i])] << endl;
				//}
				//cout << ' ' << h_iterativo.Qcalc[IDX1F(h_sistema.barrasPV[h_sistema.nPV - 1])] + h_barra.Qload[IDX1F(h_sistema.barrasPV[h_sistema.nPV - 1])] << ']' << endl;

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

			float_type maxG = maxi(h_iterativo.gLim, iterPon.ngLim); // retorna valor absoluto do elto de maior módulo em g
			if (global::verbose_mode || global::output_processo_iterativo && !global::laconic_mode) {
				printf("\nmax(g) = %e", maxG);
			}
			if (maxG <= global::tol) {
				if (global::verbose_mode) {
					printf("\n\nbreak!\nmax(g) = %e\n\n", maxi(h_iterativo.gLim, iterPon.ngLim));
				}
				return; //break; // terminou!!!
			}
			if (global::verbose_mode) {
				printf("     ...     continua:\n\n");
			}

			//checkCudaErrors(cudaDeviceSynchronize());
			//assert(cudaGetLastError() == cudaSuccess);

			//fillJstencil0based(&h_sistema, &h_sparse, &h_barra, &h_ramo, &h_iterativo);
			//checkCudaErrors(cudaMemcpy(sparsePon.spJval, h_sparse.spJval.data(), sparsePon.nnzJ * sizeof(float_type), cudaMemcpyHostToDevice));

			//mCalcJac(h_sistema, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp, h_sparse, sparsePon);
			mCalcJac(h_sistema, h_barra, h_ramo, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp, h_sparse, sparsePon, streams);

			if (global::verbose_mode) {
				// Impressão dos valores atualizados
	
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

			//  solução do sistema***************************************************************************
			{ BENCHMARK_SISTEMALINEAR
				switch (global::metodo)
				{
				case metodo::esparso:
					{
						int szJ = (h_sistema.nB - 1 + iterPon.nPQlim);

						// Calculo da forma de armazenamento esparsa:

						// parou de funcionar
						//cusparseCreateCsr(descrJ, h_iterativo.ngLim, h_iterativo.ngLim, h_sparse.nnzJ, sparsePon.csrRowPtrJ,
						//	sparsePon.cooColIndJ, sparsePon.spJval, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
						//	CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

						switch (global::metodoDeDecomposicao) {
						case metodoDeDecomposicao::LU:
							//cusolverSpDcsrlsvqr(SpHandle, szJ, nnzTotal,
							//	descrJ, d_csrValJ, d_csrRowPtrJ, d_csrColIndJ,
							//	iterPon.g, tol, 0, /*ans*/ iterPon.g, &singularity);

							//ESTA ERRADO!!!! Variaveis devem ser do host!
							//cuSOLVERstatus = cusolverSpDcsrlsvluHost(SpHandle, szJ, sparsePon.nnzJ,
							//							   			 descrJ, sparsePon.spJval, sparsePon.csrRowPtrJ, sparsePon.cooColIndJ,
							//							   			 iterPon.gLim, global::tol, 0, /*ans*/ iterPon.gLim, &singularity);     //  Remark 2: only CPU (Host) path is provided. 
							assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
							break;
						case metodoDeDecomposicao::QR:
							//cusolverSpDcsrlsvqr(SpHandle, szJ, nnzTotal,
							//	descrJ, d_csrValJ, d_csrRowPtrJ, d_csrColIndJ,
							//	iterPon.g, tol, 0, /*ans*/ iterPon.g, &singularity);
					
							cuSOLVERstatus = _cusolverSpCsrlsvqr(SpHandle, szJ, sparsePon.nnzJ,
																 descrJ, sparsePon.spJval, sparsePon.csrRowPtrJ, sparsePon.cooColIndJ,
																 iterPon.gLim, global::tol, 0, /*ans*/ iterPon.gLim, &singularity);
							assert(CUSOLVER_STATUS_SUCCESS == cuSOLVERstatus);
							break;
						}

						if (singularity != -1) {
							printf("\n\n[ERRO] J singular para a tolerancia.\n\n");
							assert(1);
						}
					}
					break;
				//case esparsoSimples:
				//	{
				//		int szJ = (h_sistema.nB - 1 + h_sistema.nPQ);

				//		// Calculo da forma de armazenamento esparsa:

				//		int* d_csrRowPtrJ = NULL;
				//		int* d_csrColIndJ = NULL;
				//		float_type* d_csrValJ = NULL;

				//		int* d_nnzPerRowColumn;
				//		int nnzTotal;

				//		checkCudaErrors(cudaMalloc((void**)&d_csrRowPtrJ, sizeof(int)* (szJ + 1)));
				//		checkCudaErrors(cudaMalloc((void**)& (d_nnzPerRowColumn), szJ * sizeof(int)));

				//		// calcula numero de eltos nao nulos
				//		cusparseDnnz(cusparseHandle, CUSPARSE_DIRECTION_ROW,
				//			szJ, szJ, descrJ, iterPon.J, szJ,
				//			d_nnzPerRowColumn, &nnzTotal);

				//		checkCudaErrors(cudaMalloc((void**)&d_csrColIndJ, sizeof(int)* (nnzTotal)));
				//		checkCudaErrors(cudaMalloc((void**)&d_csrValJ, sizeof(float_type)* (nnzTotal)));

				//		cusparseDdense2csr(cusparseHandle, szJ, szJ, descrJ, iterPon.Jlim, szJ,
				//			d_nnzPerRowColumn, d_csrValJ, d_csrRowPtrJ, d_csrColIndJ);

				//		switch (global::metodoDeDecomposicao)
				//		{
				//		case metodoLU:
				//			cusolverSpDcsrlsvqr(SpHandle, szJ, nnzTotal,
				//				descrJ, d_csrValJ, d_csrRowPtrJ, d_csrColIndJ,
				//				iterPon.g, global::tol, 0, /*ans*/ iterPon.g, &singularity);
				//			break;
				//		case metodoQR:
				//			cusolverSpDcsrlsvqr(SpHandle, szJ, nnzTotal,
				//				descrJ, d_csrValJ, d_csrRowPtrJ, d_csrColIndJ,
				//				iterPon.g, global::tol, 0, /*ans*/ iterPon.g, &singularity);
				//			break;
				//		}

				//		if (singularity != -1) {
				//			printf("\n\n[ERRO] J singular para a tolerancia.\n\n");
				//			assert(1);
				//		}

				//		// destruição das variaveis
				//		if (d_csrRowPtrJ) { cudaFree(d_csrRowPtrJ); }
				//		if (d_nnzPerRowColumn) { cudaFree(d_nnzPerRowColumn); }
				//		if (d_csrColIndJ) { cudaFree(d_csrColIndJ); }
				//		if (d_csrValJ) { cudaFree(d_csrValJ); }
				//	}
				//	break;
				case metodo::hibridoB:

					if (global::CPUsolverFlg)
					{
						h_sparse.spJval.resize(sparsePon.nnzJ);
						checkCudaErrors(cudaMemcpy(h_sparse.spJval.data(), sparsePon.spJval, sizeof(float_type) * (sparsePon.nnzJ),
							cudaMemcpyDeviceToHost));

						spSolve(&h_sistema, &h_barra, &h_ramo, &h_iterativo, &h_sparse);
						
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
					checkCudaErrors(cudaMalloc(&ipiv, sizeof(int) * /*n*/ iterPon.ngLim));

					cudaMemset(info, 0, sizeof(int));

					// a funcao getrf sobrescrevera a matriz jacobiano com L
					_cusolverDnGetrf(DnHandle, /*n*/ iterPon.ngLim, /*n*/ iterPon.ngLim, iterPon.Jlim, /*n*/ iterPon.ngLim, buffer, ipiv, info);
					checkCudaErrors(cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost));

					if (0 != h_info) {
						fprintf(stderr, "Error: LU factorization failed\n");
					}

					_cusolverDnGetrs(DnHandle, CUBLAS_OP_N, /*n*/ iterPon.ngLim, 1, iterPon.Jlim, /*n*/ iterPon.ngLim, ipiv, iterPon.gLim,
						/*n*/ iterPon.ngLim, info);
					checkCudaErrors(cudaDeviceSynchronize());

					// destruição das variaveis
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

			//  metodo***********************************************************************************

			attVOlimf(sistPon, barraPon, ramoPon, iterPon, deviceProp); // atualização das variáveis de estado

			{ BENCHMARK_CUDAMEMCPY
				checkCudaErrors(cudaMemcpy(h_barra.V, barraPon.V, sizeof(float_type)* (sistPon.nB),
					cudaMemcpyDeviceToHost));
				checkCudaErrors(cudaMemcpy(h_barra.theta, barraPon.theta, sizeof(float_type)* (sistPon.nB),
					cudaMemcpyDeviceToHost));
				BENCHMARK_SYNC
			}

			if (global::verbose_mode) {
				// Impressão dos valores atualizados

				printf("V^(%d) = \n", iterPon.iteracao + 1);
				showVec(h_barra.V, sistPon.nB, 6);

				printf("theta^(%d) = \n", iterPon.iteracao + 1);
				showVec(h_barra.theta, sistPon.nB, 6);
			}
		}
		BENCHMARK_SYNC
	}
	   
	// Finalização das variáveis
	switch (global::metodo)
	{
	case metodo::esparso:
		d_finSparse(sparsePon); // apaga variáveis esparsas da GPU
	//case esparsoSimples:

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

		// esparso (apagar se não for mais usar biblioteca cusparse)
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

//void nR(sistema &h_sistema, barra &h_barra, ramo &h_ramo, iterativo &h_iterativo, sistema* d_sistema,
//		barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema &sistPon,
//		barra &barraPon, ramo &ramoPon, iterativo &iterPon,
//		const cudaDeviceProp &deviceProp) {
//
//}