#pragma once
#include <assert.h>

// Jacobiano
// ***************************************************************************************************************

// percorre uma coluna de H
// idx = 1:nnzY
__global__ void calcHlim_Eficiente(const sistema sistPon, const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	//unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1

	if (id <= (sistPon.nnzY)) { // só faz calculos dentro das dimenções de H
		unsigned int idx = sistPon.cooRowIndY[IDX1F(id)];
		unsigned int idy = sistPon.csrColIndY[IDX1F(id)];

		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			bool offm = (idx >= sistPon.barraVO); // barra idx, posicao idx+offswing
			//idx += offm; // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2

			// idx percorre colunas de H diretamente
			// idx=1:nPQ
			// idy percorre linhas de H diretamente
			// idy=1:nPQ+nPV

			bool offk = (idy >= sistPon.barraVO);
			//idy += offk; // pula barra swing
			//printf("H(%d,%d) -> J(%d,%d)\n",idx,idy,idx-offm,idy-offk);
			if (idy != idx) {
				//fórmula pra km
				//float_type aux = phif(idx, idy, d_sistema, d_ramo);
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				d_iterativo->Jlim[IDX2F(idy - offk, idx - offm, szJ)] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.Y[IDX2F(idy, idx, sistPon.nB)].x * sin(aux) - sistPon.Y[IDX2F(idy, idx, sistPon.nB)].y * cos(aux));
				//printf("dif: H(%d,%d) -> J(%d,%d) - %f == %f\n",idx,idy,idx-offm,idy-offk,d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * ( sistPon.Y[IDX2F(idy, idx, sistPon.nB)].x*sin(aux) - sistPon.Y[IDX2F(idy, idx, sistPon.nB)].y*cos(aux)), d_iterativo->J[IDX2F(idy-offk, idx-offm, szJ)]);
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				d_iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].y;
				//printf("eq:  H(%d,%d) -> J(%d,%d) - %f == %f\n",idx,idy,idx-offm,idx-offm,-d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)]*d_barra->V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].y);
			}
		}
	}
}

// percorre uma coluna de H
// idx = 1:nB-1
__global__ void calcHlim(const sistema sistPon, const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) { //(sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1

	if (idx <= (d_sistema->nB - 1) && idy <= (d_sistema->nB - 1)) { // só faz calculos dentro das dimenções de H
		//printf("(%d,%d) -> \n",idx,idy);

		bool offm = (idx >= d_sistema->barraVO); // barra idx, posicao idx+offswing
		idx += offm; // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2

		// idx percorre colunas de H diretamente
		// idx=1:nPQ
		// idy percorre linhas de H diretamente
		// idy=1:nPQ+nPV

		bool offk = (idy >= d_sistema->barraVO);
		idy += offk; // pula barra swing
		//printf("H(%d,%d) -> J(%d,%d)\n",idx,idy,idx-offm,idy-offk);
		if (idy != idx) {
			//fórmula pra km
			//float_type aux = phif(idx, idy, d_sistema, d_ramo);
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
			d_iterativo->Jlim[IDX2F(idy - offk, idx - offm, szJ)] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].x * sin(aux) - d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].y * cos(aux));
			//printf("dif: H(%d,%d) -> J(%d,%d) - %f == %f\n",idx,idy,idx-offm,idy-offk,d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * ( d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].x*sin(aux) - d_sistema->Y[IDX2F(idy, idx, d_sistema->nB)].y*cos(aux)), d_iterativo->J[IDX2F(idy-offk, idx-offm, szJ)]);
		}
		else {
			// Hkk = -Qk - V^2_k*Bkk
			d_iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * d_sistema->Y[IDX2F(idx, idx, d_sistema->nB)].y;
			//printf("eq:  H(%d,%d) -> J(%d,%d) - %f == %f\n",idx,idy,idx-offm,idx-offm,-d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)]*d_barra->V[IDX1F(idx)] * d_sistema->Y[IDX2F(idx, idx, d_sistema->nB)].y);
		}
	}
}

// idx = 1:nPQ
__global__ void calcLlim(const iterativo iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1
	unsigned int offset = d_sistema->nB - 1; // deslocamento para armazenar elementos na matriz J

	if (idx <= (iterPon.nPQlim) && idy <= (iterPon.nPQlim)) { // só faz calculos dentro das dimenções de L


		// idy   é iterador do vetor barrasPQ, percorre linhas  de L
		// idx é iterador do vetor barrasPQ, percorre colunas de L

		//bool offm = (barrasPQ[IDX1F(idx)] > swing); // se idx>swing => armazenara elto uma posição a menos em J

		if (iterPon.barrasPQlim[IDX1F(idy)] != iterPon.barrasPQlim[IDX1F(idx)]) {
			//fórmula pra km
			float_type aux = d_barra->theta[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] - d_barra->theta[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];

			// Vk * (Gkm*sen(theta) - Bkmcos(theta))
			d_iterativo->Jlim[IDX2F(idy + offset, idx + offset, szJ)] = d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] * (d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].x * sin(aux) - d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].y * cos(aux));
		}
		else {
			// Lkk= (Qk - V^2_k*Bkk)/Vk
			d_iterativo->Jlim[IDX2F(idx + offset, idx + offset, szJ)] = (d_iterativo->Qcalc[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] - d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] * d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] * d_sistema->Y[IDX2F(iterPon.barrasPQlim[IDX1F(idx)], iterPon.barrasPQlim[IDX1F(idx)], d_sistema->nB)].y) / d_barra->V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];
		}


	}

	if (idx == 1 && idy == 1) {
		//printf("H(3,2)=%f\n",d_iterativo->J[IDX2F(2, 1, d_sistema->nB-1+d_sistema->nPQ)]);
	}


}

__global__ void calcMlim(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon) {
	unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1
	//printf("total: (idx, idy) = (%d,%d)\n", idx, idy);
	if (idx <= (sistPon.nB - 1) && idy <= (iterPon.nPQlim)) { // só faz calculos dentro das dimenções de M
		//printf("local: (idy, idx) = (%d,%d)\n", idy, idx);
		bool offm = (idx >= sistPon.barraVO); // barra idx, posicao idx+offswing
		idx += offm;

		unsigned int offset = sistPon.nB - 1; // =nPV+nPQ deslocamento para armazenar elementos na matriz J

		// idy é iterador do vetor barrasPQ, percorre linhas de M
		// idy=1:nPQ
		// idx percorre colunas de M diretamente
		// idx=1:nB - swing

		if (iterPon.barrasPQlim[IDX1F(idy)] != idx) {
			//fórmula pra km
			float_type aux = barraPon.theta[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] - barraPon.theta[IDX1F(idx)];
			//printf("dif: (idy, idx) = (%d,%d)\n", idy, idx);
			// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
			iterPon.Jlim[IDX2F(idy + offset, idx - offm, szJ)] = -barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idy)])] * barraPon.V[IDX1F(idx)] * (sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], idx, sistPon.nB)].x * cos(aux) + sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idy)], idx, sistPon.nB)].y * sin(aux));

			//printf("dif: M(%d,%d) -> J(%d,%d)\n", sistPon.barrasPQ[IDX1F(idy)], idx, idy+offset, idx-offm); //ESSE AQUI

			//printf("dif: M(%d,%d) -> J(%d,%d) - %f == %f\n"
			//	   "      (idy, idx) = (%d,%d);\n",
//				   "     V(%d - sistPon.barrasPQ[IDX1F(idy)]) = %f;\n"
//				   "     V(%d - idx) = %f;\n"
//				   "     Y(%d - sistPon.barrasPQ[IDX1F(idy)],%d - idx) = %f+j%f\n",
			//		sistPon.barrasPQ[IDX1F(idy)],idx,idy+offset, idx-offm,(-barraPon.V[IDX1F(sistPon.barrasPQ[IDX1F(idy)])]*barraPon.V[IDX1F(idx)] * ( sistPon.Y[IDX2F(sistPon.barrasPQ[IDX1F(idy)], idx, sistPon.nB)].x*cos(aux) + sistPon.Y[IDX2F(sistPon.barrasPQ[IDX1F(idy)],idx, sistPon.nB)].y*sin(aux))),iterPon.J[IDX2F(idy+offset, idx-offm, szJ)],
			//		idy, idx
					//sistPon.barrasPQ[IDX1F(idy)], barraPon.V[IDX1F(sistPon.barrasPQ[IDX1F(idy)])],
//					idx, barraPon.V[IDX1F(idx)],
//					sistPon.barrasPQ[IDX1F(idy)], idx, (sistPon.Y[IDX2F(sistPon.barrasPQ[IDX1F(idy)], idx, sistPon.nB)].x), (sistPon.Y[IDX2F(sistPon.barrasPQ[IDX1F(idy)], idx, sistPon.nB)].y));
			//		);
			//printf("dif: J-> %p\n", iterPon.J);
			//printf("dif: V(%d - sistPon.barrasPQ[IDX1F(idy)]) = %f; V(%d - idx) = %f == %f; Y(%d - sistPon.barrasPQ[IDX1F(idy)],%d - idx) = %f\n",sistPon.barrasPQ[IDX1F(idy)], barraPon.V[IDX1F(sistPon.barrasPQ[IDX1F(idy)])], idx, sistPon.barrasPQ[IDX1F(idy)], idx, sistPon.Y[IDX2F(sistPon.barrasPQ[IDX1F(idy)], idx, sistPon.nB)].x);
		}
		else {
			// elemento 2,2 da não fica gravado..........
			// Mkk= (Pk - V^2_k*Gkk)
			iterPon.Jlim[IDX2F(idy + offset, idx - offm, szJ)] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].x;
			//printf("eq:  M(%d,%d) -> J(%d,%d) - %f == %f\n",idx,idx,idy+offset, idx-offm, iterPon.J[IDX2F(idx+offset, idx-offm, szJ)], iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)]*barraPon.V[IDX1F(idx)] * sistPon.Y[IDX2F(idx, idx, sistPon.nB)].x);

			//printf("eq:  M(%d,%d) -> J(%d,%d)\n", idx, idx, idy+offset, idx-offm); //ESSE AQUI
		}
	}
}

// idx = 1:nPQ
__global__ void calcNlim(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon){
	unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1
	//printf("total: (idx, idy) = (%d,%d)\n", idx, idy);
	if(idx<=(iterPon.nPQlim) && idy<=(sistPon.nB-1)){ // só faz calculos dentro das dimenções de N
		//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1 // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2
		//printf("loacal: (idx, idy) = (%d,%d)\n", idx, idy);
		unsigned int offset = sistPon.nB-1; // =nPV+nPQ deslocamento para armazenar elementos na matriz J (colunas)

		// idx é iterador do vetor barrasPQ, percorre colunas de N
		// idx=1:nPQ
		// idy percorre linhas de N diretamente
		// idy=1:nPQ+nPV

		bool offk = (idy >= sistPon.barraVO);
		idy += offk;
		if (idy != iterPon.barrasPQlim[IDX1F(idx)]){
			//fórmula pra km=iidx
			float_type aux = barraPon.theta[IDX1F(idy)] -
					barraPon.theta[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];

			// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
			iterPon.Jlim[IDX2F(idy-offk, idx+offset, szJ)] = barraPon.V[IDX1F(idy)] *
					(sistPon.Y[IDX2F(idy, iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].x*cos(aux) +
					 sistPon.Y[IDX2F(idy, iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].y*sin(aux));
			//printf("dif: (idy, idx) = (%d,%d); N(%d,%d) -> J(%d,%d) - %f == %f\n",idy,idx,idy,sistPon.barrasPQ[IDX1F(idx)],idy-offk, idx+offset, (iterPon.J[IDX2F(idy-offk, idx+offset, szJ)]), (barraPon.V[IDX1F(idy)] * (sistPon.Y[IDX2F(idy, sistPon.barrasPQ[IDX1F(idx)], sistPon.nB)].x*cos(aux) + sistPon.Y[IDX2F(idy, sistPon.barrasPQ[IDX1F(idx)], sistPon.nB)].y*sin(aux))));
		}
		else{
			// Nkk= (Qk - V^2_k*Bkk)/Vkk
			iterPon.Jlim[IDX2F(idy-offk, idx+offset, szJ)] =
					(iterPon.Pcalc[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] +
					 barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])]*barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])] *
					 sistPon.Y[IDX2F(iterPon.barrasPQlim[IDX1F(idx)], iterPon.barrasPQlim[IDX1F(idx)], sistPon.nB)].x)/
					 barraPon.V[IDX1F(iterPon.barrasPQlim[IDX1F(idx)])];
		}
	}
}

void calcHfLim_eficiente(sistema& sistPon, iterativo& iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	unsigned int tamanho = sistPon.nnzY;

	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / (float)16 * deviceprop.warpSize), 1);

	calcHlim_Eficiente <<<dimGrid, dimBlock >>> (sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcHfLim(sistema& sistPon, iterativo& iterPon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	unsigned int tamanho = sistPon.nB - 1;

	dim3 dimBlock(deviceprop.warpSize, 16);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / (float)deviceprop.warpSize), (unsigned int)ceil(((float)tamanho) / (float)16));

	calcHlim <<<dimGrid, dimBlock >>> (sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcLfLim(iterativo& iterPon, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop) {
	unsigned int tamanho = iterPon.nPQlim;


	dim3 dimBlock(deviceprop.warpSize, deviceprop.warpSize);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / (float)deviceprop.warpSize), (unsigned int)ceil(((float)tamanho) / (float)deviceprop.warpSize));

	calcLlim <<<dimGrid, dimBlock >>> (iterPon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void calcMfLim(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	unsigned int tamanhoX = sistPon.nB - 1;
	unsigned int tamanhoY = iterPon.nPQlim;

	// printf("size      = (%d,%d)\n", tamanhoX, tamanhoY);
	// printf("blocksize = (%d,%d)\n", deviceprop.warpSize, deviceprop.warpSize);
	// printf("gridsize  = (%d,%d)\n", (unsigned int)ceil(( (float)tamanhoX )/(float)deviceprop.warpSize), (unsigned int)ceil(( (float)tamanhoY )/(float)deviceprop.warpSize) );

	// dim3 dimBlock(deviceprop.warpSize, deviceprop.warpSize);
	// dim3 dimGrid((unsigned int)ceil(((float)tamanhoX)/(float)deviceprop.warpSize), (unsigned int)ceil(((float)tamanhoY)/(float)deviceprop.warpSize));

	dim3 dimBlock(32, 16);
	dim3 dimGrid((unsigned int)ceil(((float)tamanhoX) / (float)32), (unsigned int)ceil(((float)tamanhoY) / (float)16));

	//calcM<<<dimGrid, dimBlock>>>(sistPon, barraPon, ramoPon, iterPon);
	calcMlim <<<dimGrid, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);
	//cudaDeviceSynchronize();
}

void calcNfLim(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	unsigned int tamanhoX = iterPon.nPQlim;
	unsigned int tamanhoY = sistPon.nB - 1;

	// dim3 dimBlock(deviceprop.warpSize, deviceprop.warpSize); // cada bloco é 32x32 em GPU típica
	// dim3 dimGrid((unsigned int)ceil(((float)tamanhoX)/(float)deviceprop.warpSize), (unsigned int)ceil(((float)tamanhoY)/(float)deviceprop.warpSize));

	dim3 dimBlock(16, 32); // cada bloco é 32x32 em GPU típica
	dim3 dimGrid((unsigned int)ceil(((float)tamanhoX) / (float)16), (unsigned int)ceil(((float)tamanhoY) / (float)32));

	calcNlim <<<dimGrid, dimBlock >>> (sistPon, barraPon, ramoPon, iterPon);
	//	cudaDeviceSynchronize();
}

void calcJacLim(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	//printf("Jacobiano_0 =\n");
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//calcHf(h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);

	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcHfLim(sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	// cudaDeviceSynchronize();
	// printf("JacobianoH =\n");
	// float_type h_J[h_sistema.nB*h_sistema.nB]; // copia de d_J
	// checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	// showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************
	calcLfLim(iterPon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	//printf("JacobianoL =\n");
	//checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************
	calcMfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//TESTE***************************************
	//printf("JacobianoM =\n");
	//checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************

	calcNfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	if (global::verbose_mode) {
		printf("JacobianoN =\n");
		checkCudaErrors(cudaMemcpy(h_iterativo.Jlim, iterPon.Jlim, sizeof(float_type) * (h_sistema.nB - 1 + iterPon.nPQlim) * (h_sistema.nB - 1 + iterPon.nPQlim), cudaMemcpyDeviceToHost));
		showMat(h_iterativo.Jlim, h_sistema.nB - 1 + iterPon.nPQlim);
	}
	//********************************************
}

void calcJacLim_H_eficiente(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, cudaDeviceProp deviceprop) {
	//printf("Jacobiano_0 =\n");
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//calcHf(h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);

	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);

	calcHfLim_eficiente(sistPon, iterPon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	// cudaDeviceSynchronize();
	// printf("JacobianoH =\n");
	// float_type h_J[h_sistema.nB*h_sistema.nB]; // copia de d_J
	// checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	// showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************
	calcLfLim(iterPon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	//printf("JacobianoL =\n");
	//checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************
	calcMfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//TESTE***************************************
	//printf("JacobianoM =\n");
	//checkCudaErrors(cudaMemcpy(h_iterativo.J, iterPon.J, sizeof(float_type)*(h_sistema.nB-1+h_sistema.nPQ)*(h_sistema.nB-1+h_sistema.nPQ), cudaMemcpyDeviceToHost));
	//showMat(h_iterativo.J, h_sistema.nB-1 + h_sistema.nPQ);
	//********************************************

	calcNfLim(sistPon, barraPon, ramoPon, iterPon, deviceprop);
	checkCudaErrors(cudaDeviceSynchronize());
	assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	if (global::verbose_mode) {
		printf("JacobianoN =\n");
		checkCudaErrors(cudaMemcpy(h_iterativo.Jlim, iterPon.Jlim, sizeof(float_type) * (h_sistema.nB - 1 + iterPon.nPQlim) * (h_sistema.nB - 1 + iterPon.nPQlim), cudaMemcpyDeviceToHost));
		showMat(h_iterativo.Jlim, h_sistema.nB - 1 + iterPon.nPQlim);
	}
	//********************************************
}


// *********************************************************************************************************************************************************************************************************************************************
// * versão esparsa ****************************************************************************************************************************************************************************************************************************
// *********************************************************************************************************************************************************************************************************************************************

//// idx para cada nnzH
__global__ void SpCalcH(const sistema sistPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
	//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);
	//__syncthreads();
	if (id < (sparsePon.nnzH)) { // só faz calculos dentro das dimenções de H
		//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);

		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based
		//printf("(%d, %d)\n\n", idx, idy);

		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			//bool offm = (idx >= sistPon.barraVO); // barra idx, posicao idx+offswing
			//idx += offm; // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2

			// idx percorre colunas de H diretamente
			// idx=1:nPQ
			// idy percorre linhas de H diretamente
			// idy=1:nPQ+nPV

			//bool offk = (idy >= sistPon.barraVO);
			//idy += offk; // pula barra swing
			//printf("H(%d,%d) -> J(%d,%d)\n",idx,idy,idx-offm,idy-offk);
			if (idy != idx) {
				//fórmula pra km
				//float_type aux = phif(idx, idy, d_sistema, d_ramo);
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				//printf("H_(%d,%d) = %f*%f*(%f*sin(%f-%f) - %f*cos(%f-%f))\nH_(%d,%d) = %f\n", idy, idx, d_barra->V[IDX1F(idy)], d_barra->V[IDX1F(idx)], sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x, d_barra->theta[IDX1F(idy)], d_barra->theta[IDX1F(idx)], sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y, d_barra->theta[IDX1F(idy)], d_barra->theta[IDX1F(idx)], idy, idx, d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux)));
				sparsePon.spJval[sparsePon.Hpos[id]] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				sparsePon.spJval[sparsePon.Hpos[id]] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y;
			}
		}
	}
}

// uso de memoria compartilhada (deu errado a cópia...)
//__global__ void SpCalcH(const sistema sistPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
//	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
//	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
//	// unsigned int idx = threadIdx.x;
//	//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);
//	//__syncthreads();
//	if (id < (sparsePon.nnzH)) { // só faz calculos dentro das dimenções de H
//		//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);
//
//		extern __shared__ int s_csrRowPtrY[];
//		int* s_csrColIndY = s_csrRowPtrY + sistPon.nB;
//
//		
//		for (size_t i = threadIdx.x; i <= sistPon.nB; i+= blockDim.x)	{
//			s_csrRowPtrY[i] = sistPon.csrRowPtrY[i];
//		}
//
//		for (size_t i = threadIdx.x; i < sistPon.nnzY; i += blockDim.x) {
//			s_csrColIndY[i] = sistPon.csrColIndY[i];
//		}
//		
//
//		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based
//		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based
//		//printf("(%d, %d)\n\n", idx, idy);
//
//		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
//			//bool offm = (idx >= sistPon.barraVO); // barra idx, posicao idx+offswing
//			//idx += offm; // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2
//
//			// idx percorre colunas de H diretamente
//			// idx=1:nPQ
//			// idy percorre linhas de H diretamente
//			// idy=1:nPQ+nPV
//
//			//bool offk = (idy >= sistPon.barraVO);
//			//idy += offk; // pula barra swing
//			//printf("H(%d,%d) -> J(%d,%d)\n",idx,idy,idx-offm,idy-offk);
//			if (idy != idx) {
//				//fórmula pra km
//				float_type aux = phif(idx, idy, d_sistema, d_ramo);
//				aux += d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];
//
//				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
//				//printf("H_(%d,%d) = %f*%f*(%f*sin(%f-%f) - %f*cos(%f-%f))\nH_(%d,%d) = %f\n", idy, idx, d_barra->V[IDX1F(idy)], d_barra->V[IDX1F(idx)], sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x, d_barra->theta[IDX1F(idy)], d_barra->theta[IDX1F(idx)], sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y, d_barra->theta[IDX1F(idy)], d_barra->theta[IDX1F(idx)], idy, idx, d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux)));
//				sparsePon.spJval[sparsePon.Hpos[id]] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), s_csrRowPtrY, s_csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), s_csrRowPtrY, s_csrColIndY, sistPon.nnzY)].y * cos(aux));
//			}
//			else {
//				// Hkk = -Qk - V^2_k*Bkk
//				sparsePon.spJval[sparsePon.Hpos[id]] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), s_csrRowPtrY, s_csrColIndY, sistPon.nnzY)].y;
//			}
//		}
//	}
//}

__global__ void SpCalcL(const sistema sistPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;

	//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);
	//__syncthreads();
	if (id < sparsePon.nnzL) { // só faz calculos dentro das dimenções de L
		//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);
		//__syncthreads();
		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Lpos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Lpos[id]] + 1;
		//printf("(%d, %d)\n\n", idx, idy);
		//__syncthreads();

		// idy   é iterador do vetor barrasPQ, percorre linhas  de L
		// idx é iterador do vetor barrasPQ, percorre colunas de L

		if (idy != idx) {
			//fórmula pra km
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			// Vk * (Gkm*sen(theta) - Bkmcos(theta))
			sparsePon.spJval[sparsePon.Lpos[id]] = d_barra->V[IDX1F(idy)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
		}
		else {
			// Lkk= (Qk - V^2_k*Bkk)/Vk
			sparsePon.spJval[sparsePon.Lpos[id]] = (d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y) / d_barra->V[IDX1F(idx)];
		}
	}
}

__global__ void SpCalcM(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon, const d_sparse sparsePon) {
	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	//unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1

	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;

	//printf("<<<%d, %d, %d>>> (%d)\n", threadIdx.x, blockDim.x, blockIdx.x, id);

	if (id < sparsePon.nnzM) { // só faz calculos dentro das dimenções de M

		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Mpos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Mpos[id]] + 1;

		//bool offm = (idx >= sistPon.barraVO); // barra idx, posicao idx+offswing
		//idx += offm;

		// unsigned int offset = sistPon.nB - 1; // =nPV+nPQ deslocamento para armazenar elementos na matriz J

		// idy é iterador do vetor barrasPQ, percorre linhas de M
		// idy=1:nPQ
		// idx percorre colunas de M diretamente
		// idx=1:nB - swing

		if (idy != idx) {
			//fórmula pra km
			float_type aux = barraPon.theta[IDX1F(idy)] - barraPon.theta[IDX1F(idx)];
			// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
			sparsePon.spJval[sparsePon.Mpos[id]] = -barraPon.V[IDX1F(idy)] * barraPon.V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) + sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			// Mkk= (Pk - V^2_k*Gkk)
			sparsePon.spJval[sparsePon.Mpos[id]] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x;
		}
	}
}

__global__ void SpCalcN(const sistema sistPon, const barra barraPon, const ramo ramoPon, const iterativo iterPon, const d_sparse sparsePon) {
	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
	//unsigned int idy = threadIdx.y + blockDim.y * blockIdx.y + 1; //começa de 1

	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < sparsePon.nnzN) { // só faz calculos dentro das dimenções de N
		
		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Npos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Npos[id]] + 1;
		
		//printf("(%d, %d)\n\n", idx, idy);

		// unsigned int offset = sistPon.nB - 1;

		// idx é iterador do vetor barrasPQ, percorre colunas de N
		// idx=1:nPQ
		// idy percorre linhas de N diretamente
		// idy=1:nPQ+nPV

		//bool offk = (idy >= sistPon.barraVO);
		//idy += offk;
		if (idy != idx) {
			//fórmula pra km=iidx
			float_type aux = barraPon.theta[IDX1F(idy)] -
				barraPon.theta[IDX1F(idx)];

			// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
			sparsePon.spJval[sparsePon.Npos[id]] = barraPon.V[IDX1F(idy)] *
				(sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) +
					sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			// Nkk= (Qk - V^2_k*Bkk)/Vkk
			sparsePon.spJval[sparsePon.Npos[id]] =
				(iterPon.Pcalc[IDX1F(idx)] +
					barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] *
					sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x) /
				barraPon.V[IDX1F(idx)];
		}
	}
}

void SpCalcHf(sistema& sistPon, iterativo& iterPon, d_sparse& sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	unsigned int tamanho = sparsePon.nnzH;

	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / ((float)16 * deviceprop.warpSize)), 1);

	SpCalcH <<<dimGrid, dimBlock, 0 /*(sistPon.nnzY + sistPon.nB + 1) * sizeof(unsigned int)*/, streams[0]>>> (sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void SpCalcLf(sistema& sistPon, iterativo& iterPon, d_sparse& sparsePon, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	unsigned int tamanho = sparsePon.nnzL;


	dim3 dimBlock(16 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / ((float)16 * deviceprop.warpSize)), 1);

	SpCalcL <<<dimGrid, dimBlock, 0, streams[1]>>> (sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

void SpCalcMf(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	unsigned int tamanho = sparsePon.nnzM;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcM <<<dimGrid, dimBlock, 0, streams[2]>>> (sistPon, barraPon, ramoPon, iterPon, sparsePon);
}

void SpCalcNf(sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	unsigned int tamanho = sparsePon.nnzN;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcN <<<dimGrid, dimBlock, 0, streams[3]>>> (sistPon, barraPon, ramoPon, iterPon, sparsePon);
}

void SpCalcJac(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	
	//checkCudaErrors(cudaDeviceSynchronize());
	//assert(cudaGetLastError() == cudaSuccess);

	SpCalcHf(sistPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop, streams);
	//checkCudaErrors(cudaDeviceSynchronize());
	//assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	if (global::verbose_mode) {
		printf("H^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}
	//********************************************
	SpCalcLf(sistPon, iterPon, sparsePon, h_sistema, d_sistema, d_barra, d_ramo, d_iterativo, deviceprop, streams);
	//checkCudaErrors(cudaDeviceSynchronize());
	//assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
	//TESTE***************************************
	if (global::verbose_mode) {
		printf("L^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}
	//********************************************
	SpCalcMf(sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceprop, streams);
	//checkCudaErrors(cudaDeviceSynchronize());
	//assert(cudaGetLastError() == cudaSuccess);
	//TESTE***************************************
	if (global::verbose_mode) {
		printf("M^(%d) =\n", iterPon.iteracao);
		float_type* aux = (float_type*)malloc(sparsePon.nnzJ * sizeof(float_type));
		checkCudaErrors(cudaMemcpy(aux, sparsePon.spJval, sparsePon.nnzJ * sizeof(float_type), cudaMemcpyDeviceToHost));
		showVec(aux, sparsePon.nnzJ, 5);
		free(aux);
	}
	//********************************************

	SpCalcNf(sistPon, barraPon, ramoPon, iterPon, sparsePon, deviceprop, streams);
	//checkCudaErrors(cudaDeviceSynchronize());
	//assert(cudaGetLastError() == cudaSuccess);
	//cudaDeviceSynchronize();
}

// **********************************************

__global__ void SpCalcJUno(const sistema sistPon, const barra barraPon, const iterativo iterPon, const d_sparse sparsePon, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo) {
	//unsigned int szJ = iterPon.nPVlim + iterPon.nPQlim + iterPon.nPQlim;
	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;

	if (id < (sparsePon.nnzH)) { // H
		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Hpos[id]] + 1; // 1 based

		if ((idx != sistPon.barraVO) && (idy != sistPon.barraVO)) {
			// idx percorre colunas de H diretamente
			// idx=1:nPQ
			// idy percorre linhas de H diretamente
			// idy=1:nPQ+nPV

			if (idy != idx) {
				//fórmula pra km
				//float_type aux = phif(idx, idy, d_sistema, d_ramo);
				float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				sparsePon.spJval[sparsePon.Hpos[id]] = d_barra->V[IDX1F(idy)] * d_barra->V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				sparsePon.spJval[sparsePon.Hpos[id]] = -d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y;
			}
		}
	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL)) { // L
		id -= sparsePon.nnzH;

		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Lpos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Lpos[id]] + 1;

		// idy   é iterador do vetor barrasPQ, percorre linhas  de L
		// idx é iterador do vetor barrasPQ, percorre colunas de L

		if (idy != idx) {
			//fórmula pra km
			float_type aux = d_barra->theta[IDX1F(idy)] - d_barra->theta[IDX1F(idx)];

			// Vk * (Gkm*sen(theta) - Bkmcos(theta))
			sparsePon.spJval[sparsePon.Lpos[id]] = d_barra->V[IDX1F(idy)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * sin(aux) - sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * cos(aux));
		}
		else {
			// Lkk= (Qk - V^2_k*Bkk)/Vk
			sparsePon.spJval[sparsePon.Lpos[id]] = (d_iterativo->Qcalc[IDX1F(idx)] - d_barra->V[IDX1F(idx)] * d_barra->V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y) / d_barra->V[IDX1F(idx)];
		}

	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM)) { // M
		id -= sparsePon.nnzH + sparsePon.nnzL;

		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Mpos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Mpos[id]] + 1;

		// idy é iterador do vetor barrasPQ, percorre linhas de M
		// idy=1:nPQ
		// idx percorre colunas de M diretamente
		// idx=1:nB - swing

		if (idy != idx) {
			//fórmula pra km
			float_type aux = barraPon.theta[IDX1F(idy)] - barraPon.theta[IDX1F(idx)];
			// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
			sparsePon.spJval[sparsePon.Mpos[id]] = -barraPon.V[IDX1F(idy)] * barraPon.V[IDX1F(idx)] * (sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) + sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			// Mkk= (Pk - V^2_k*Gkk)
			sparsePon.spJval[sparsePon.Mpos[id]] = iterPon.Pcalc[IDX1F(idx)] - barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] * sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x;
		}
	}
	else if (id < (sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM + sparsePon.nnzN)) { // N
		id -= sparsePon.nnzH + sparsePon.nnzL + sparsePon.nnzM;

		unsigned int idy = sparsePon.cooRowIndSubMatJ[sparsePon.Npos[id]] + 1;
		unsigned int idx = sparsePon.cooColIndSubMatJ[sparsePon.Npos[id]] + 1;

		// idx é iterador do vetor barrasPQ, percorre colunas de N
		// idx=1:nPQ
		// idy percorre linhas de N diretamente
		// idy=1:nPQ+nPV

		//idy += offk;
		if (idy != idx) {
			//fórmula pra km=iidx
			float_type aux = barraPon.theta[IDX1F(idy)] -
				barraPon.theta[IDX1F(idx)];

			// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
			sparsePon.spJval[sparsePon.Npos[id]] = barraPon.V[IDX1F(idy)] *
				(sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux) +
					sistPon.spYval[coeffPos(IDX1F(idy), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		}
		else {
			// Nkk= (Qk - V^2_k*Bkk)/Vkk
			sparsePon.spJval[sparsePon.Npos[id]] =
				(iterPon.Pcalc[IDX1F(idx)] +
					barraPon.V[IDX1F(idx)] * barraPon.V[IDX1F(idx)] *
					sistPon.spYval[coeffPos(IDX1F(idx), IDX1F(idx), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x) /
				barraPon.V[IDX1F(idx)];
		}
	}
}

void SpCalcJUnof(sistema& h_sistema, iterativo& h_iterativo, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, sistema& sistPon, barra& barraPon, ramo& ramoPon, iterativo& iterPon, d_sparse& sparsePon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	unsigned int tamanho = sparsePon.nnzJ;

	dim3 dimBlock(8 * deviceprop.warpSize, 1);
	dim3 dimGrid((unsigned int)ceil(((float)tamanho) / ((float)8 * deviceprop.warpSize)), 1);

	SpCalcJUno <<<dimGrid, dimBlock>>> (sistPon, barraPon, iterPon, sparsePon, d_sistema, d_barra, d_ramo, d_iterativo);
}

// **********************************************

//__global__ void stencilGPU(const sistema sistPon, const barra barraPon, const iterativo iterPon, const d_sparse sparsePon) {
//	unsigned int id = threadIdx.x + blockDim.x * blockIdx.x;
//	if (id < (sistPon.nB - 1 + sistPon.nPQ)) {
//
//	}
//}

//void stencilGPUf(sistema& sistPon, barra& barraPon, iterativo& iterPon, d_sparse& sparsePon) {
//	checkCudaErrors(cudaMalloc(&(sparsePon.cooColIndJ), 4 * sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.cooColIndJ, 0, 4 * sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.cooRowIndJ), 4 * sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.cooRowIndJ, 0, 4 * sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.csrRowPtrJ), sistPon.nB * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.csrRowPtrJ, 0, sistPon.nB * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.cooColIndSubMatJ), 4 * sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.cooColIndSubMatJ, 0, 4 * sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.cooRowIndSubMatJ), 4 * sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.cooRowIndSubMatJ, 0, 4 * sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.Hpos), sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.Hpos, 0, sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.Lpos), sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.Lpos, 0, sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.Mpos), sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.Mpos, 0, sistPon.nnzY * sizeof(unsigned int)));
//
//	checkCudaErrors(cudaMalloc(&(sparsePon.Npos), sistPon.nnzY * sizeof(unsigned int)));
//	checkCudaErrors(cudaMemset(sparsePon.Npos, 0, sistPon.nnzY * sizeof(unsigned int)));
//}