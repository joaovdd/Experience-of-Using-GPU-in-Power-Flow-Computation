#pragma once
#include <assert.h>

__global__ void calcHlim_Eficiente(const sistema sistPon, const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	int id = threadIdx.x + blockDim.x * blockIdx.x + 1; 

	if (id <= (sistPon.nnzY)) { 
		int idx = sistPon.cooRowIndY[IDX1F(id)];
		int idy = sistPon.csrColIndY[IDX1F(id)];

		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			bool offm = (idx >= sistPon.barraVO); 
			
			bool offk = (idy >= sistPon.barraVO);
			
			if (idy != idx) {
				
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				d_iterativo->Jlim[IDX2F(idy - offk, idx - offm, szJ)] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.Y[IDX2F(idy, idx, sistPon.nB)].x * sin(aux) - sistPon.Y[IDX2F(idy, idx, sistPon.nB)].y * cos(aux));
				
			}
			else {
				
				d_iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].y;
				
			}
		}
	}
}

__global__ void calcHlim(const sistema sistPon, const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) { 
	int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 
	int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; 

	if (idx <= (d_sistema->nB - 1) && idy <= (d_sistema->nB - 1)) { 
		
		bool offm = (idx >= d_sistema->barraVO); 
		idx += offm; 

		bool offk = (idy >= d_sistema->barraVO);
		idy += offk; 
		
		if (idy != idx) {
			
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			d_iterativo->Jlim[IDX2F(idy - offk, idx - offm, szJ)] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].x * sin(aux) - d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].y * cos(aux));
			
		}
		else {
			
			d_iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * d_sistema->Y[IDX2F(idx, idx, d_sistema->nB)].y;
			
		}
	}
}

__global__ void calcLlim(const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 
	int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; 
	int offset = d_sistema->nB - 1; 

	if (idx <= (iterPon.nPQlim) && idy <= (iterPon.nPQlim)) { 

		if (iterPon.barrasPQlim[IDX1F(idy)] != iterPon.barrasPQlim[IDX1F(idx)]) {
			
			float_type aux = d_barra->theta[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] - d_barra->theta[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];

			d_iterativo->Jlim[IDX2F(idy + offset, idx + offset, szJ)] = d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] * (d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].x * sin(aux) - d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].y * cos(aux));
		}
		else {
			
			d_iterativo->Jlim[IDX2F(idx + offset, idx + offset, szJ)] = (d_iterativo->Qcalc[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] - d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] * d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] * d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idx)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].y) / d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];
		}

	}

	if (idx == 1 && idy == 1) {
		
	}

}

__global__ void calcMlim(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon) {
	int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 
	int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; 

	if (idx <= (sistPon.nB - 1) && idy <= (iterPon.nPQlim)) { 
		
		bool offm = (idx >= sistPon.barraVO); 
		idx += offm;

		int offset = sistPon.nB - 1; 

		if (iterPon.barrasPQlim[IDX1F(idy)] != idx) {
			
			float_type aux = barraPon.theta[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] - barraPon.theta[IDX1F(idx)];
			
			iterPon.Jlim[IDX2F(idy + offset, idx - offm, szJ)] = -barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] * barraPon.V[IDX1F(idx)] * (sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], idx, sistPon.nB)].x * cos(aux) + sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], idx, sistPon.nB)].y * sin(aux));

		}
		else {
			
			iterPon.Jlim[IDX2F(idy + offset, idx - offm, szJ)] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].x;
			
		}
	}
}

__global__ void calcNlim(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon){
	int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; 
	int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; 

	if(idx<=(iterPon.nPQlim) && idy<=(sistPon.nB-1)){ 
		
		int offset = sistPon.nB-1; 

		bool offk = (idy >= sistPon.barraVO);
		idy += offk;
		if (idy != iterPon.barrasPQlim[IDX1F(idx)]){
			
			float_type aux = barraPon.theta[IDX1F(idy)] -
					barraPon.theta[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];

			iterPon.Jlim[IDX2F(idy-offk, idx+offset, szJ)] = barraPon.V[IDX1F(idy)] *
					(sistPon.Y[IDX2F(idy, iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].x*cos(aux) +
					 sistPon.Y[IDX2F(idy, iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].y*sin(aux));
			
		}
		else{
			
			iterPon.Jlim[IDX2F(idy-offk, idx+offset, szJ)] =
					(iterPon.Pcalc[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] +
					 barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])]*barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] *
					 sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idx)], iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].x)/
					 barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];
		}
	}
}

void calcHfLim_eficiente(sistema& sistPon, iterativo& iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nnzY;

	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / (float)16 * deviceprop.warpSize), 1);

	calcHlim_Eficiente <<<dimGrid, dimBlock >>> (sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcHfLim(sistema& sistPon, iterativo& iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nB - 1;

	dim3 dimBlock(deviceprop.warpSize, 16);
	dim3 dimGrid((int)ceil(((float)tamanho) / (float)deviceprop.warpSize), (int)ceil(((float)tamanho) / (float)16));

	calcHlim <<<dimGrid, dimBlock >>> (sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcLfLim(iterativo& iterPon, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	int tamanho = iterPon.nPQlim;

	dim3 dimBlock(deviceprop.warpSize, deviceprop.warpSize);
	dim3 dimGrid((int)ceil(((float)tamanho) / (float)deviceprop.warpSize), (int)ceil(((float)tamanho) / (float)deviceprop.warpSize));

	calcLlim <<<dimGrid, dimBlock >>> (iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcMfLim(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	int tamanhoX = sistPon.nB - 1;
	int tamanhoY = iterPon.nPQlim;

	dim3 dimBlock(32, 16);
	dim3 dimGrid((int)ceil(((float)tamanhoX) / (float)32), (int)ceil(((float)tamanhoY) / (float)16));

	calcMlim <<<dimGrid, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

}

void calcNfLim(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	int tamanhoX = iterPon.nPQlim;
	int tamanhoY = sistPon.nB - 1;

	dim3 dimBlock(16, 32); 
	dim3 dimGrid((int)ceil(((float)tamanhoX) / (float)16), (int)ceil(((float)tamanhoY) / (float)32));

	calcNlim <<<dimGrid, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);

}

void calcJacLim(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcHfLim(sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcLfLim(iterPon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcMfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcNfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	if (global::verbose_mode) {
		printf("JacobianoN =\n");
		checkCudaErrors(cudaMemcpy(h_iterativo.Jlim, iterPon.Jlim, sizeof(float_type) * (h_sistema.nB - 1 + iterPon.nPQlim) * (h_sistema.nB - 1 + iterPon.nPQlim), cudaMemcpyDeviceToHost));
		showMat(h_iterativo.Jlim, h_sistema.nB - 1 + iterPon.nPQlim);
	}

}

void calcJacLim_H_eficiente(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcHfLim_eficiente(sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcLfLim(iterPon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcMfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcNfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	if (global::verbose_mode) {
		printf("JacobianoN =\n");
		checkCudaErrors(cudaMemcpy(h_iterativo.Jlim, iterPon.Jlim, sizeof(float_type) * (h_sistema.nB - 1 + iterPon.nPQlim) * (h_sistema.nB - 1 + iterPon.nPQlim), cudaMemcpyDeviceToHost));
		showMat(h_iterativo.Jlim, h_sistema.nB - 1 + iterPon.nPQlim);
	}

}

__global__ void SpCalcH(const sistema sistPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < (sparsePon.nnzH)) { 
		
		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Hpos[id]] + 1; 
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Hpos[id]] + 1; 
		
		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			
			if (idy != idx) {
				
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				sparsePon.spJval[sparsePon.Hpos[id]] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
			}
			else {
				
				sparsePon.spJval[sparsePon.Hpos[id]] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y;
			}
		}
	}
}

__global__ void SpCalcL(const sistema sistPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < sparsePon.nnzL) { 
		
		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Lpos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Lpos[id]] + 1;
		
		if (idy != idx) {
			
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			sparsePon.spJval[sparsePon.Lpos[id]] = d_barra->V[IDX1F(idy)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Lpos[id]] = (d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y) / d_barra->V[IDX1F(idx)];
		}
	}
}

__global__ void SpCalcM(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon, const d_sparse sparsePon) {
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < sparsePon.nnzM) { 

		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Mpos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Mpos[id]] + 1;

		if (idy != idx) {
			
			float_type aux = barraPon.theta[IDX1F(idy)] - barraPon.theta[IDX1F(idx)];
			
			sparsePon.spJval[sparsePon.Mpos[id]] = -barraPon.V[IDX1F(idy)] * barraPon.V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) + sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Mpos[id]] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x;
		}
	}
}

__global__ void SpCalcN(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon, const d_sparse sparsePon) {
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < sparsePon.nnzN) { 
		
		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Npos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Npos[id]] + 1;
		
		if (idy != idx) {
			
			float_type aux = barraPon.theta[IDX1F(idy)] -
				barraPon.theta[IDX1F(idx)];

			sparsePon.spJval[sparsePon.Npos[id]] = barraPon.V[IDX1F(idy)] *
				(sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) +
					sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Npos[id]] =
				(iterPon.Pcalc[IDX1F(idx)] +
					barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] *
					sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x) /
				barraPon.V[IDX1F(idx)];
		}
	}
}

void SpCalcHf(sistema& sistPon, iterativo& iterPon, d_sparse& sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sparsePon.nnzH;

	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / ((float)16 * deviceprop.warpSize)), 1);

	SpCalcH <<<dimGrid, dimBlock, 0 , streams[0]>>> (sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void SpCalcLf(sistema& sistPon, iterativo& iterPon, d_sparse& sparsePon, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sparsePon.nnzL;

	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / ((float)16 * deviceprop.warpSize)), 1);

	SpCalcL <<<dimGrid, dimBlock, 0, streams[1]>>> (sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void SpCalcMf(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sparsePon.nnzM;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcM <<<dimGrid, dimBlock, 0, streams[2]>>> (sistPon, barraPon, ramoPon, iterPon, sparsePon);
}

void SpCalcNf(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sparsePon.nnzN;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcN <<<dimGrid, dimBlock, 0, streams[3]>>> (sistPon, barraPon, ramoPon, iterPon, sparsePon);
}

void SpCalcJac(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	SpCalcHf(sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop, streams);

	if (global::verbose_mode) {
		printf("H^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}

	SpCalcLf(sistPon, iterPon, sparsePon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop, streams);

	if (global::verbose_mode) {
		printf("L^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}

	SpCalcMf(sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceprop, streams);

	if (global::verbose_mode) {
		printf("M^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}

	SpCalcNf(sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceprop, streams);

}

__global__ void SpCalcJUno(const sistema sistPon, const barra barraPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < (sparsePon.nnzH)) { 
		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Hpos[id]] + 1; 
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Hpos[id]] + 1; 

		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			
			if (idy != idx) {
				
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				sparsePon.spJval[sparsePon.Hpos[id]] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
			}
			else {
				
				sparsePon.spJval[sparsePon.Hpos[id]] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y;
			}
		}
	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL)) { 
		id -= sparsePon.nnzH;

		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Lpos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Lpos[id]] + 1;

		if (idy != idx) {
			
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			sparsePon.spJval[sparsePon.Lpos[id]] = d_barra->V[IDX1F(idy)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Lpos[id]] = (d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y) / d_barra->V[IDX1F(idx)];
		}

	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM)) { 
		id -= sparsePon.nnzH + sparsePon.nnzL;

		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Mpos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Mpos[id]] + 1;

		if (idy != idx) {
			
			float_type aux = barraPon.theta[IDX1F(idy)] - barraPon.theta[IDX1F(idx)];
			
			sparsePon.spJval[sparsePon.Mpos[id]] = -barraPon.V[IDX1F(idy)] * barraPon.V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) + sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Mpos[id]] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x;
		}
	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM + sparsePon.nnzN)) { 
		id -= sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM;

		int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Npos[id]] + 1;
		int idx = sparsePon.cooColIndSubMatJ[sparsePon.Npos[id]] + 1;

		if (idy != idx) {
			
			float_type aux = barraPon.theta[IDX1F(idy)] -
				barraPon.theta[IDX1F(idx)];

			sparsePon.spJval[sparsePon.Npos[id]] = barraPon.V[IDX1F(idy)] *
				(sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) +
					sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			
			sparsePon.spJval[sparsePon.Npos[id]] =
				(iterPon.Pcalc[IDX1F(idx)] +
					barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] *
					sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x) /
				barraPon.V[IDX1F(idx)];
		}
	}
}

void SpCalcJUnof(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sparsePon.nnzJ;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcJUno <<<dimGrid, dimBlock>>> (sistPon, barraPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

