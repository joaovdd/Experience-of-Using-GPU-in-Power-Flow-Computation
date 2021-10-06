#pragma once
#include "sistema.h"
#include "esparso.h"

// checa limite de inje��o de reativos e cria vetor limQ
void chkLimQ(sistema* sistema, barra* barra, iterativo* iterativo) {
	// prints não modificados...
	printf("\n");
	#pragma omp parallel for if (openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
			// Qg = Qliq - Qload
			if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - global::tol)) {
				if (!iterativo->limQ[i]) {
					printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
				}
				else if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] != sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
					printf("Barra %d [%d]: PQ (limite inferior) -> PQ (limite superior)\n", sistema->barrasPV[i], i);
				}
				// limite reativo inferior violado: Barra se torna PQ
				iterativo->limQ[i] = 1;
				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
			}
			else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + global::tol)) {
				if (!iterativo->limQ[i]) {
					printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
				}
				else if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] != sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
					printf("Barra %d [%d]: PQ (limite superior) -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
				}
				// limite reativo superior violado: barra se torna PQ
				iterativo->limQ[i] = 1;
				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
			}
			else {
				if (iterativo->limQ[i]) {
					printf("Barra %d [%d]: PQ -> PV\n", sistema->barrasPV[i], i);
				}
				iterativo->limQ[i] = 0;
				//iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = 0;
			}
		}
}

// preenche vetor sistema->barrasPQlim com as barras PQ e as PV que atingiram os limites de 
// reativos. Sempre em ordem crescente.
void geraVetPQ(sistema* sistema, barra* barra, iterativo* iterativo) {
	int pntPQ = 0; // aponta para a pr�xima entrada do vetor sistema->nPQ a ser copiada para sistema->barrasPQlim
	int pntPQlim = 0; // aponta para a pr�xima entrada do vetor sistema->barrasPQlim a ser preenchida
	//#pragma omp parallel for if (openmp)
	for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
		if (iterativo->limQ[i]) { // encontrou o �ndice da barra PQ "novata"
			while ((sistema->barrasPQ[pntPQ] < sistema->barrasPV[i]) && (pntPQ < sistema->nPQ)) { // j percorre barras PQ diretamente.
																								  // enquanto forem inferiores � novata
																								  // e n�o tiver terminado de percorrer o vetor barrasPQ.
				iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; // Preenche vetor barrasPQlim com barras inferiores � "novata".
				pntPQ++;
				pntPQlim++;
			}
			iterativo->barrasPQlim[pntPQlim] = sistema->barrasPV[i]; // sistema->barrasPQ[aux] > sistema->barrasPV[i]
			pntPQlim++;
		}
	}
	while (pntPQ < sistema->nPQ) { // se n�o tiver terminado de percorrer o vetor barrasPQ anteriormente
		iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; // Preenche vetor barrasPQlim com barras inferiores � "novata".
		pntPQ++;
		pntPQlim++;
	}
	iterativo->nPQlim = pntPQlim;
}

// preenche vetor sistema->barrasPVlim com barras PV que n�o atingiram os limites de reativos
void geraVetPV(sistema* sistema, barra* barra, iterativo* iterativo) {
	iterativo->nPVlim = 0;
	//#pragma omp parallel for if (openmp)
	for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
		if (!(iterativo->limQ[i])) {
			iterativo->barrasPVlim[iterativo->nPVlim] = sistema->barrasPV[i];
			iterativo->nPVlim++;
		}
	}
}

d_sparse d_initSparse(sistema& h_sistema, h_sparse& h_sparse) {

	d_sparse aux; // estrutura-ponteiro para d_barra

	checkCudaErrors(cudaMalloc(&(aux.csrRowPtrJ), 2 * (h_sistema.nB - 1) * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.csrRowPtrJ, 0, 2 * (h_sistema.nB - 1) * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.cooColIndJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooColIndJ, 0, h_sparse.nnzJ * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.spJval), h_sparse.nnzJ * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.spJval, 0, h_sparse.nnzJ * sizeof(float_type)));

	aux.nnzJ = h_sparse.nnzJ;

	// stencil
	checkCudaErrors(cudaMalloc(&(aux.cooColIndSubMatJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooColIndSubMatJ, 0, h_sparse.nnzJ * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.cooRowIndSubMatJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooRowIndSubMatJ, 0, h_sparse.nnzJ * sizeof(int)));

	// HLMNpos
	aux.nnzH = h_sparse.Hpos.size();
	aux.nnzL = h_sparse.Lpos.size();
	aux.nnzM = h_sparse.Mpos.size();
	aux.nnzN = h_sparse.Npos.size();

	checkCudaErrors(cudaMalloc(&(aux.Hpos), aux.nnzH * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Hpos, 0, aux.nnzH * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Lpos), aux.nnzL * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Lpos, 0, aux.nnzL * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Mpos), aux.nnzM * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Mpos, 0, aux.nnzM * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Npos), aux.nnzN * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Npos, 0, aux.nnzN * sizeof(int)));

	return aux;
}

void d_finSparse(d_sparse aux) {
	if (aux.csrRowPtrJ) { cudaFree(aux.csrRowPtrJ); }
	if (aux.cooColIndJ) { cudaFree(aux.cooColIndJ); }
	if (aux.spJval) { cudaFree(aux.spJval); }
	if (aux.cooColIndSubMatJ) { cudaFree(aux.cooColIndSubMatJ); }
	if (aux.cooRowIndSubMatJ) { cudaFree(aux.cooRowIndSubMatJ); }
	if (aux.Hpos) { cudaFree(aux.Hpos); }
	if (aux.Lpos) { cudaFree(aux.Lpos); }
	if (aux.Mpos) { cudaFree(aux.Mpos); }
	if (aux.Npos) { cudaFree(aux.Npos); }
}

void SpJacCpyH2D(sistema* sistema, h_sparse* h_sparse, d_sparse& sparsePon) {
	checkCudaErrors(cudaMemcpy(sparsePon.cooColIndJ, h_sparse->cooColIndJ.data(), sparsePon.nnzJ * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(sparsePon.csrRowPtrJ, h_sparse->csrRowPtrJ.data(), h_sparse->csrRowPtrJ.size() * sizeof(int), cudaMemcpyHostToDevice));

	// stencil
	checkCudaErrors(cudaMemcpy(sparsePon.cooColIndSubMatJ, h_sparse->cooColIndSubMatJ.data(), sparsePon.nnzJ * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(sparsePon.cooRowIndSubMatJ, h_sparse->cooRowIndSubMatJ.data(), sparsePon.nnzJ * sizeof(int), cudaMemcpyHostToDevice));

	// HLMNpos
	checkCudaErrors(cudaMemcpy(sparsePon.Hpos, h_sparse->Hpos.data(), h_sparse->Hpos.size() * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(sparsePon.Lpos, h_sparse->Lpos.data(), h_sparse->Lpos.size() * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(sparsePon.Mpos, h_sparse->Mpos.data(), h_sparse->Mpos.size() * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(sparsePon.Npos, h_sparse->Npos.data(), h_sparse->Npos.size() * sizeof(int), cudaMemcpyHostToDevice));
}

void d_attSparse(sistema& h_sistema, h_sparse& h_sparse, d_sparse& aux) {

	if (aux.cooColIndJ) { cudaFree(aux.cooColIndJ); }
	if (aux.spJval) { cudaFree(aux.spJval); }

	checkCudaErrors(cudaMalloc(&(aux.cooColIndJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooColIndJ, 0, h_sparse.nnzJ * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.spJval), h_sparse.nnzJ * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.spJval, 0, h_sparse.nnzJ * sizeof(float_type)));

	if (aux.cooRowIndSubMatJ) { cudaFree(aux.cooRowIndSubMatJ); }
	if (aux.cooColIndSubMatJ) { cudaFree(aux.cooColIndSubMatJ); }

	checkCudaErrors(cudaMalloc(&(aux.cooRowIndSubMatJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooRowIndSubMatJ, 0, h_sparse.nnzJ * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.cooColIndSubMatJ), h_sparse.nnzJ * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.cooColIndSubMatJ, 0, h_sparse.nnzJ * sizeof(int)));

	aux.nnzJ = h_sparse.nnzJ;

	// HLMNpos
	if (aux.Hpos) { cudaFree(aux.Hpos); }
	if (aux.Lpos) { cudaFree(aux.Lpos); }
	if (aux.Mpos) { cudaFree(aux.Mpos); }
	if (aux.Npos) { cudaFree(aux.Npos); }

	checkCudaErrors(cudaMalloc(&(aux.Hpos), h_sparse.Hpos.size() * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Hpos, 0, h_sparse.Hpos.size() * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Lpos), h_sparse.Lpos.size() * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Lpos, 0, h_sparse.Lpos.size() * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Mpos), h_sparse.Mpos.size() * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Mpos, 0, h_sparse.Mpos.size() * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.Npos), h_sparse.Npos.size() * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.Npos, 0, h_sparse.Npos.size() * sizeof(int)));

	aux.nnzH = h_sparse.Hpos.size();
	aux.nnzL = h_sparse.Lpos.size();
	aux.nnzM = h_sparse.Mpos.size();
	aux.nnzN = h_sparse.Npos.size();

	// atualiza
	SpJacCpyH2D(&h_sistema, &h_sparse, aux);
}

void initLim(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo& iterPon, iterativo* iterativo) {
	// Barras PQ
	#pragma omp parallel for if (openmp)
		for (int i = 0; i < sistema->nPQ; i++) { // percorre barras PQ diretamente
			iterativo->barrasPQlim[i] = sistema->barrasPQ[i];
		}
	iterativo->nPQlim = sistema->nPQ;
	iterPon.nPQlim = sistema->nPQ;
	checkCudaErrors(cudaMemcpy(iterPon.barrasPQlim, iterativo->barrasPQlim, (sistema->nPV + sistema->nPQ) * sizeof(int), cudaMemcpyHostToDevice)); // ver transferencia device-device

	// Barras PV
	#pragma omp parallel for if (openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV diretamente
			iterativo->barrasPVlim[i] = sistema->barrasPV[i];
		}
	iterativo->nPVlim = sistema->nPV;
	iterPon.nPVlim = sistema->nPV;
	checkCudaErrors(cudaMemcpy(iterPon.barrasPVlim, iterativo->barrasPVlim, sistema->nPV * sizeof(int), cudaMemcpyHostToDevice)); // ver transferencia device-device

	iterativo->ngLim = sistema->nPQ + sistema->nPQ + sistema->nPV;
	iterPon.ngLim = sistema->nPQ + sistema->nPQ + sistema->nPV;

	#pragma omp parallel for if (openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV diretamente
			iterativo->limQ[i] = 0;
		}
	//checkCudaErrors(cudaMalloc(&(iterPon.limQ), h_sistema.nPV * sizeof(float_type)));
	//checkCudaErrors(cudaMemset(iterPon.limQ, 0, h_sistema.nPV * sizeof(float_type)));

	#pragma omp parallel for if (openmp)
		for (int i = 0; i < sistema->nB; i++) {
			iterativo->QliqLim[i] = barra->Qliq[i];
		}

	checkCudaErrors(cudaMemcpy(iterPon.QliqLim, iterativo->QliqLim, sistema->nB * sizeof(float_type), cudaMemcpyHostToDevice));

	if (global::metodo == metodo::esparso || global::metodo == metodo::hibridoA || global::metodo == metodo::hibridoB) { BENCHMARK_JACOBIANOSTENCIL_BUILD
		//Jstencil0based(sistema, h_sparse, barra, iterativo); // calcula o estêncil da matriz jacobiano
		//reordJac(sistema , h_sparse, barra, iterativo);
		//Jstencil(sistema, h_sparse, barra, iterativo); // calcula o estêncil da matriz jacobiano com listas de trabalho já ordenadas para minimização da divergência
		// Jstencil_eficiente_(sistema, h_sparse, barra, iterativo);
		Jstencil_eficiente(sistema, h_sparse, barra, iterativo);
	}
}

void printAlLim(sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo) {
	cout << "nPQlim = " << iterativo.nPQlim << endl;
	cout << "nPVlim = " << iterativo.nPVlim << endl;

	cout << "limQ:" << endl;
	showVec(iterativo.limQ, sistema.nPV, 4);
	cout << "barrasPVlim:" << endl;
	showVec(iterativo.barrasPVlim, iterativo.nPVlim, 4);
	cout << "barrasPQlim:" << endl;
	showVec(iterativo.barrasPQlim, iterativo.nPQlim, 4);
	cout << "limQinf:" << endl;
	showVec(sistema.limQinf, sistema.nPV, 4);
	cout << "limQsup:" << endl;
	showVec(sistema.limQsup, sistema.nPV, 4);
	cout << "QliqLim:" << endl;
	showVec(iterativo.QliqLim, sistema.nB, 4);
}