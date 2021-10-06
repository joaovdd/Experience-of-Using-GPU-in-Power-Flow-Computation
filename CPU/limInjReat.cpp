#include "limInjReat.h"

// preenche vetor sistema->barrasPVlim com barras PV que não atingiram os limites de reativos
void geraVetPV(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	iterativo->nPVlim = 0;
	//#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
			if (!(iterativo->limQ[i])) {
				iterativo->barrasPVlim[iterativo->nPVlim] = sistema->barrasPV[i];
				iterativo->nPVlim++;
			}
		}	
}

// preenche vetor sistema->barrasPQlim com as barras PQ e as PV que atingiram os limites de 
// reativos. Sempre em ordem crescente.
void geraVetPQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int pntPQ = 0; // aponta para a próxima entrada do vetor sistema->nPQ a ser copiada para sistema->barrasPQlim
	int pntPQlim = 0; // aponta para a próxima entrada do vetor sistema->barrasPQlim a ser preenchida
	//#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
			if (iterativo->limQ[i]) { // encontrou o índice da barra PQ "novata"
				while ((sistema->barrasPQ[pntPQ] < sistema->barrasPV[i]) && (pntPQ < sistema->nPQ)) { // j percorre barras PQ diretamente.
																									  // enquanto forem inferiores à novata
																									  // e não tiver terminado de percorrer o vetor barrasPQ.
					iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; // Preenche vetor barrasPQlim com barras inferiores à "novata".
					pntPQ++;
					pntPQlim++;
				}
				iterativo->barrasPQlim[pntPQlim] = sistema->barrasPV[i]; // sistema->barrasPQ[aux] > sistema->barrasPV[i]
				pntPQlim++;
			}
		}
	while (pntPQ < sistema->nPQ) { // se não tiver terminado de percorrer o vetor barrasPQ anteriormente
		iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; // Preenche vetor barrasPQlim com barras inferiores à "novata".
		pntPQ++;
		pntPQlim++;
	}
	iterativo->nPQlim = pntPQlim;
}

//// checa limite de injeção de reativos e cria vetor limQ
//void chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
//	//printf("\n");
//#pragma omp parallel for if (global::openmp)
//	for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
//		// Qg = Qliq - Qload
//		if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - global::tol)) {
//			if (!iterativo->limQ[i]) {
//				printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
//			}
//			else if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] != sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
//				printf("Barra %d [%d]: PQ (limite inferior) -> PQ (limite superior)\n", sistema->barrasPV[i], i);
//			}
//			// limite reativo inferior violado: Barra se torna PQ
//			iterativo->limQ[i] = 1;
//			iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
//		}
//		else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + global::tol)) {
//			if (!iterativo->limQ[i]) {
//				printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
//			}
//			else if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] != sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
//				printf("Barra %d [%d]: PQ (limite superior) -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
//			}
//			// limite reativo superior violado: barra se torna PQ
//			iterativo->limQ[i] = 1;
//			iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
//		}
//		else {
//			//if (iterativo->limQ[i]) {
//			//	//printf("Barra %d [%d]: PQ -> PV\n", sistema->barrasPV[i], i);
//			//}
//			//else {
//			//	iterativo->limQ[i] = 0;
//			//}
//			if (iterativo->limQ[i]) {
//				printf("Barra %d [%d]: PQ -> PV\n", sistema->barrasPV[i], i);
//			}
//			iterativo->limQ[i] = 0;
//			//iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = 0;
//		}
//	}
//}

// checa limite de injeção de reativos e cria vetor limQ
void chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	//printf("\n");
	auto folga = global::tol_limInjReat; // 0.000001; ///*10**/global::tol;
	//int cnt = 0;
	//int max = 9999;

	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
			if (!iterativo->limQ[i]) { // se iterativo->limQ[i] == 0
				// Qg = Qliq + Qload
				if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - folga)) {
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
					}

					// limite reativo inferior violado: Barra se torna PQ
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
					iterativo->limQ[i] = 1;
					//cnt++;
				}
				else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + folga)) {
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
					}

					// limite reativo superior violado: barra se torna PQ
					iterativo->limQ[i] = 1;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
					//cnt++;
				}
			}
			else { // iterativo->limQ[i] == 1
				// verificar se a tensão não precisa mais ser controlada:
				if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] == sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
					// se o limite de injeção de reativos atingido foi o superior
					if (barra->V[IDX1F(sistema->barrasPV[i])] > sistema->VfixadaPV[i] + folga) {
						// e se a tensão voltou a patamar superior ao fixado
						if (global::output_processo_iterativo && !global::laconic_mode) {
							printf("Barra %d [%d]: PQ (limite superior) -> PV\n", sistema->barrasPV[i], i);
						}
						iterativo->limQ[i] = 0;
						iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
						//cnt++;
					}
				}
				else {
					// se o limite de injeção de reativos atingido foi o inferior (i.i., não é o sup)
					if (barra->V[IDX1F(sistema->barrasPV[i])] < sistema->VfixadaPV[i] - folga) {
						// e se a tensão voltou a patamar inferior ao fixado
						if (global::output_processo_iterativo && !global::laconic_mode) {
							printf("Barra %d [%d]: PQ (limite inferior) -> PV\n", sistema->barrasPV[i], i);
						}
						iterativo->limQ[i] = 0;
						iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
						//cnt++;
					}
				}
			}
			//if (cnt == max) {
			//	break;
			//}
		}
}

// checa limite de injeção de reativos de barras que já violaram e, se necessário, atualiza vetor limQ 
bool undo_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	
	bool ans = 0;
	auto folga = global::tol_limInjReat;//0.00001; ///*10*/global::tol;
	//int cnt = 0;
	//int max = 9999;

#pragma omp parallel for if (global::openmp)
	for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
		if (iterativo->limQ[i]) { // se iterativo->limQ[i] == 1
			// verificar se a tensão não precisa mais ser controlada:
			if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] == sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
				// se o limite de injeção de reativos atingido foi o superior
				if (barra->V[IDX1F(sistema->barrasPV[i])] > sistema->VfixadaPV[i] + folga) {
					// e se a tensão voltou a patamar inferior ao fixado
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PQ (limite superior) -> PV\n", sistema->barrasPV[i], i);
					}
					iterativo->limQ[i] = 0;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
					ans = 1;
					//cnt++;
				}
			}
			else {
				// se o limite de injeção de reativos atingido foi o inferior (i.i., não é o sup)
				if (barra->V[IDX1F(sistema->barrasPV[i])] < sistema->VfixadaPV[i] - folga) {
					// e se a tensão voltou a patamar inferior ao fixado
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PQ (limite inferior) -> PV\n", sistema->barrasPV[i], i);
					}
					iterativo->limQ[i] = 0;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
					ans = 1;
					//cnt++;
				}
			}
		}
		//if (cnt == max) {
		//	break;
		//}
	}
	return ans;
}

void do_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	//printf("\n");
#pragma omp parallel for if (global::openmp)
	for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV
		if (!iterativo->limQ[i]) { // se iterativo->limQ[i] == 0
			// Qg = Qliq - Qload
			if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - global::tol)) {
				if (global::output_processo_iterativo && !global::laconic_mode) {
					printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
				}

				// limite reativo inferior violado: Barra se torna PQ
				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
				iterativo->limQ[i] = 1;
			}
			else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + global::tol)) {
				if (global::output_processo_iterativo && !global::laconic_mode) {
					printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
				}

				// limite reativo superior violado: barra se torna PQ
				iterativo->limQ[i] = 1;
				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
			}
		}
	}
}

// Calcula o resíduo com limites de injeção de reativos
void calcResLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	// deltaPs
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i < sistema->barraVO; i++) { //antes da swing; i percorre barras
			// deltaP = |Pcalc - Pesp|
			iterativo->gLim[IDX1F(i)] = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
			//printf("g[%d] = ", i-1);
			//printf("%e = |%e - %e| ", iterativo->g[IDX1F(i)], iterativo->Pcalc[IDX1F(i)], barra->Pliq[IDX1F(i)]);
			//printf("g[%d] = %e\n", i, iterativo->g[IDX1F(i)]);
		}

	#pragma omp parallel for if (global::openmp)
		for (int i = sistema->barraVO + 1; i <= sistema->nB; i++) { // depois da swing; i percorre barras
			// deltaP = |Pcalc - Pesp|
			float_type aux = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
			iterativo->gLim[IDX1F(i) - 1] = aux;
			//printf("g[%d] = ", i-1-1);
			//printf("%e = |%14.6e - %14.6e| ", iterativo->g[IDX1F(i)-1], iterativo->Pcalc[IDX1F(i)], barra->Pliq[IDX1F(i)]);
			//printf("g[%d] = %e\n", i-1, iterativo->g[IDX1F(i) - 1]);
		}

	// deltaQs
	const int offset = sistema->nB - 1;
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < iterativo->nPQlim; i++) { // i é iterador das barras PQ
			// deltaQ = |Qcalc - Qesp|
			iterativo->gLim[i + offset] = iterativo->QliqLim[IDX1F(iterativo->barrasPQlim[i])] - iterativo->Qcalc[IDX1F(iterativo->barrasPQlim[i])];
			//printf("g[%d] = ", i+sistema->nPQ+sistema->nPV);
			//printf("%e = |%14.6e - %14.6e| ", iterativo->g[i+sistema->nPQ+sistema->nPV], iterativo->Qcalc[IDX1F(sistema->barrasPQ[i])], barra->Qliq[IDX1F(sistema->barrasPQ[i])]);
			//printf("g[%d] = %14.6e\n", i+sistema->nPQ+sistema->nPV, iterativo->g[i+sistema->nPQ+sistema->nPV]);
		}
}

// **********************************************************************************************************************************************************************************

void initLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	// Barras PQ
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPQ; i++) { // percorre barras PV diretamente
			iterativo->barrasPQlim[i] = sistema->barrasPQ[i];
		}
	iterativo->nPQlim = sistema->nPQ;
	// Barras PV
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV diretamente
			iterativo->barrasPVlim[i] = sistema->barrasPV[i];
		}
	iterativo->nPVlim = sistema->nPV;

	iterativo->ngLim = sistema->nPQ + sistema->nPQ + sistema->nPV;

	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { // percorre barras PV diretamente
			iterativo->limQ[i] = 0;
		}

	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nB; i++) { // percorre barras PV diretamente
			iterativo->QliqLim[i] = barra->Qliq[i];
		}

#if JACOBIANO_ESPARSO_STENCIL
	if (global::metodo == esparso) {
		{ BENCHMARK_JACOBIANOSTENCIL_BUILD
			Jstencil0based(sistema, barra, ramo, iterativo);
		}
	}
#endif
}

// **********************************************************************************************************************************************************************************

void attVOlim(sistemaType* sistema, barraType* barra, iterativoType* iterativo) {
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i < sistema->barraVO; i++) {  // i percorre as barras diretamente e pula a barra swing
			barra->theta[IDX1F(i)] += iterativo->gLim[IDX1F(i)]; // atualiza valor de theta; g é igual a -deltaX!
			//printf("attVO: theta[%d] = %f\n", i, barra->theta[IDX1F(i)]);
		}
	#pragma omp parallel for if (global::openmp)
		for (int i = sistema->barraVO + 1; i <= sistema->nB; i++) { // i percorre as barras diretamente e pula a barra swing
			barra->theta[IDX1F(i)] += iterativo->gLim[IDX1F(i) - 1];           // atualiza valor de theta; g é igual a -deltaX!
			//printf("attVO: theta[%d] = %f\n", i, barra->theta[IDX1F(i)]);
		}

	const int offset = sistema->nB - 1;
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < iterativo->nPQlim; i++) {                 // i é iterador das barras PQ
			barra->V[IDX1F(iterativo->barrasPQlim[i])] += iterativo->gLim[i + offset]; // atualiza valor de V; g é igual a -deltaX!
			//printf("attVO: V[%d] = %f\n", sistema->barrasPQ[i], barra->V[IDX1F(sistema->barrasPQ[i])]);
		}

	if (global::verbose_mode) {
		printf("\nV^(%d):\n", iterativo->iteracao);
		showVec(barra->V, sistema->nB, 6);
		printf("\ntheta^(%d):\n", iterativo->iteracao);
		showVec(barra->theta, sistema->nB, 6);
	}

}

void printAlLim(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	std::cout << "nPQlim = " << iterativo.nPQlim << std::endl;
	std::cout << "nPVlim = " << iterativo.nPVlim << std::endl;
	
	std::cout << "limQ:" << std::endl;
	showVec(iterativo.limQ, sistema.nPV, 4);
	std::cout << "barrasPVlim:" << std::endl;
	showVec(iterativo.barrasPVlim, iterativo.nPVlim, 4);
	std::cout << "barrasPQlim:" << std::endl;
	showVec(iterativo.barrasPQlim, iterativo.nPQlim, 4);
	std::cout << "limQinf:" << std::endl;
	showVec(sistema.limQinf, sistema.nPV, 4);
	std::cout << "limQsup:" << std::endl;
	showVec(sistema.limQsup, sistema.nPV, 4);
	std::cout << "QliqLim:" << std::endl;
	showVec(iterativo.QliqLim, sistema.nB, 4);
}
