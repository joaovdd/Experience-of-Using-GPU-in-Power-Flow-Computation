#include "limInjReat.h"

void geraVetPV(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	iterativo->nPVlim = 0;

		for (int i = 0; i < sistema->nPV; i++) { 
			if (!(iterativo->limQ[i])) {
				iterativo->barrasPVlim[iterativo->nPVlim] = sistema->barrasPV[i];
				iterativo->nPVlim++;
			}
		}	
}

void geraVetPQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int pntPQ = 0; 
	int pntPQlim = 0; 

		for (int i = 0; i < sistema->nPV; i++) { 
			if (iterativo->limQ[i]) { 
				while ((sistema->barrasPQ[pntPQ] < sistema->barrasPV[i]) && (pntPQ < sistema->nPQ)) { 
																									  
					iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; 
					pntPQ++;
					pntPQlim++;
				}
				iterativo->barrasPQlim[pntPQlim] = sistema->barrasPV[i]; 
				pntPQlim++;
			}
		}
	while (pntPQ < sistema->nPQ) { 
		iterativo->barrasPQlim[pntPQlim] = sistema->barrasPQ[pntPQ]; 
		pntPQ++;
		pntPQlim++;
	}
	iterativo->nPQlim = pntPQlim;
}

void chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	auto folga = global::tol_limInjReat;
	
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) { 
			if (!iterativo->limQ[i]) { 
				if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - folga)) {
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
					}

					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
					iterativo->limQ[i] = 1;
				}
				else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + folga)) {
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
					}

					iterativo->limQ[i] = 1;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
					
				}
			}
			else { 
				
				if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] == sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
					
					if (barra->V[IDX1F(sistema->barrasPV[i])] > sistema->VfixadaPV[i] + folga) {
						
						if (global::output_processo_iterativo && !global::laconic_mode) {
							printf("Barra %d [%d]: PQ (limite superior) -> PV\n", sistema->barrasPV[i], i);
						}
						iterativo->limQ[i] = 0;
						iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
						
					}
				}
				else {
					
					if (barra->V[IDX1F(sistema->barrasPV[i])] < sistema->VfixadaPV[i] - folga) {
						
						if (global::output_processo_iterativo && !global::laconic_mode) {
							printf("Barra %d [%d]: PQ (limite inferior) -> PV\n", sistema->barrasPV[i], i);
						}
						iterativo->limQ[i] = 0;
						iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
						
					}
				}
			}
			
		}
}

bool undo_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	bool ans = 0;
	auto folga = global::tol_limInjReat;

#pragma omp parallel for if (global::openmp)
	for (int i = 0; i < sistema->nPV; i++) { 
		if (iterativo->limQ[i]) { 
			
			if (iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] == sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])]) {
				
				if (barra->V[IDX1F(sistema->barrasPV[i])] > sistema->VfixadaPV[i] + folga) {
					
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PQ (limite superior) -> PV\n", sistema->barrasPV[i], i);
					}
					iterativo->limQ[i] = 0;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
					ans = 1;
					
				}
			}
			else {
				
				if (barra->V[IDX1F(sistema->barrasPV[i])] < sistema->VfixadaPV[i] - folga) {
					
					if (global::output_processo_iterativo && !global::laconic_mode) {
						printf("Barra %d [%d]: PQ (limite inferior) -> PV\n", sistema->barrasPV[i], i);
					}
					iterativo->limQ[i] = 0;
					iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = barra->Qliq[IDX1F(sistema->barrasPV[i])];
					ans = 1;
					
				}
			}
		}
		
	}
	return ans;
}

void do_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
#pragma omp parallel for if (global::openmp)
	for (int i = 0; i < sistema->nPV; i++) { 
		if (!iterativo->limQ[i]) { 
			
			if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) < (sistema->limQinf[i] - global::tol)) {
				if (global::output_processo_iterativo && !global::laconic_mode) {
					printf("Barra %d [%d]: PV -> PQ (limite superior)\n", sistema->barrasPV[i], i);
				}

				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQinf[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
				iterativo->limQ[i] = 1;
			}
			else if ((iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]) > (sistema->limQsup[i] + global::tol)) {
				if (global::output_processo_iterativo && !global::laconic_mode) {
					printf("Barra %d [%d]: PV -> PQ (limite inferior)\n", sistema->barrasPV[i], i);
				}

				iterativo->limQ[i] = 1;
				iterativo->QliqLim[IDX1F(sistema->barrasPV[i])] = sistema->limQsup[i] - barra->Qload[IDX1F(sistema->barrasPV[i])];
			}
		}
	}
}

void calcResLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i < sistema->barraVO; i++) { 
			iterativo->gLim[IDX1F(i)] = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
		}

	#pragma omp parallel for if (global::openmp)
		for (int i = sistema->barraVO + 1; i <= sistema->nB; i++) { 
			float_type aux = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
			iterativo->gLim[IDX1F(i) - 1] = aux;
		}

	const int offset = sistema->nB - 1;
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < iterativo->nPQlim; i++) { 
			iterativo->gLim[i + offset] = iterativo->QliqLim[IDX1F(iterativo->barrasPQlim[i])] - iterativo->Qcalc[IDX1F(iterativo->barrasPQlim[i])];
		}
}

void initLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPQ; i++) {
			iterativo->barrasPQlim[i] = sistema->barrasPQ[i];
		}
	iterativo->nPQlim = sistema->nPQ;

	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nPV; i++) {
			iterativo->barrasPVlim[i] = sistema->barrasPV[i];
		}
	iterativo->nPVlim = sistema->nPV;

	iterativo->ngLim = sistema->nPQ + sistema->nPQ + sistema->nPV;

		for (int i = 0; i < sistema->nPV; i++) { 
			iterativo->limQ[i] = 0;
		}

		for (int i = 0; i < sistema->nB; i++) { 
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

void attVOlim(sistemaType* sistema, barraType* barra, iterativoType* iterativo) {
	#pragma omp parallel for if (global::openmp)
		for (int i = 1; i < sistema->barraVO; i++) {
			barra->theta[IDX1F(i)] += iterativo->gLim[IDX1F(i)];
		}
	#pragma omp parallel for if (global::openmp)
		for (int i = sistema->barraVO + 1; i <= sistema->nB; i++) {
			barra->theta[IDX1F(i)] += iterativo->gLim[IDX1F(i) - 1]; 
		}

	const int offset = sistema->nB - 1;
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < iterativo->nPQlim; i++) { 
			barra->V[IDX1F(iterativo->barrasPQlim[i])] += iterativo->gLim[i + offset]; 
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
