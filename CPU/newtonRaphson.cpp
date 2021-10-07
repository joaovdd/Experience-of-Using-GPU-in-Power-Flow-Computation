#include "newtonRaphson.h"

using namespace std;
using namespace Eigen;

float_type maxi(float_type* vet, int dim){
	float_type aux = 0;
	for (int i = 0; i < dim; i++){
		if (abs(vet[i]) > aux) {
			aux = abs(vet[i]);
			
		}
	}
	return aux;
}

void calcRes(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo){
	for (int i = 1; i < sistema->barraVO; i++){ 
		
		iterativo->g[IDX1F(i)] = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
		
	}

	for (int i = sistema->barraVO+1; i <= sistema->nB; i++){ 
		
		float_type aux = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
		iterativo->g[IDX1F(i) - 1] = aux;
		
	}

	const int offset = sistema->nB-1;
	for (int i = 0; i < sistema->nPQ; i++){ 
		
		iterativo->g[i+offset] = barra->Qliq[IDX1F(sistema->barrasPQ[i])] - iterativo->Qcalc[IDX1F(sistema->barrasPQ[i])];
		
	}
}

void attVO(sistemaType* sistema, barraType* barra, iterativoType* iterativo){
	for (int i = 1; i < sistema->barraVO; i++){ 
		barra->theta[IDX1F(i)] += iterativo->g[IDX1F(i)];  
		
	}
	for (int i = sistema->barraVO+1; i <= sistema->nB; i++) { 
		barra->theta[IDX1F(i)] += iterativo->g[IDX1F(i)-1];              
		
	}

	const int offset = sistema->nB - 1;
	for (int i = 0; i < sistema->nPQ; i++) {                    
		barra->V[IDX1F(sistema->barrasPQ[i])] += iterativo->g[i + offset]; 
		
	}

	if (global::verbose_mode) {
		printf("\nV^(%d):\n", iterativo->iteracao);
		showVec(barra->V, sistema->nB, 6);
		printf("\ntheta^(%d):\n", iterativo->iteracao);
		showVec(barra->theta, sistema->nB, 6);
	}

}

void evalLimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	chkLimQ(sistema, barra, ramo, iterativo);   
	geraVetPV(sistema, barra, ramo, iterativo); 
	geraVetPQ(sistema, barra, ramo, iterativo); 
	iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; 

#if JACOBIANO_ESPARSO_STENCIL
	if (global::metodo == esparso) {
		
		{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
			
			sistema->spJval.clear();
			sistema->spJsubmatType.clear();
			sistema->cooColIndSubMatJ.clear();
			sistema->cooRowIndSubMatJ.clear();
			sistema->cooColIndJ.clear();
			sistema->cooRowIndJ.clear();
			
			Jstencil0based(sistema, barra, ramo, iterativo);
		}
	}
#endif
}

void do_LimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	do_chkLimQ(sistema, barra, ramo, iterativo);   
	geraVetPV(sistema, barra, ramo, iterativo); 
	geraVetPQ(sistema, barra, ramo, iterativo); 
	iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; 

#if JACOBIANO_ESPARSO_STENCIL

	if (global::metodo == esparso) {
		
		{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
			
			sistema->spJval.clear();
			sistema->spJsubmatType.clear();
			sistema->cooColIndSubMatJ.clear();
			sistema->cooRowIndSubMatJ.clear();
			sistema->cooColIndJ.clear();
			sistema->cooRowIndJ.clear();
			
			Jstencil0based(sistema, barra, ramo, iterativo);
		}

	}
#endif
}

void undo_LimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	if (undo_chkLimQ(sistema, barra, ramo, iterativo)) {
		geraVetPV(sistema, barra, ramo, iterativo); 
		geraVetPQ(sistema, barra, ramo, iterativo); 
		iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; 
		
#if JACOBIANO_ESPARSO_STENCIL
		if (global::metodo == esparso) {
			
			{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
				
				sistema->spJval.clear();
				sistema->spJsubmatType.clear();
				sistema->cooColIndSubMatJ.clear();
				sistema->cooRowIndSubMatJ.clear();
				sistema->cooColIndJ.clear();
				sistema->cooRowIndJ.clear();
				
				Jstencil0based(sistema, barra, ramo, iterativo);
			}

		}
#endif
	}
}

bool simpleChkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	auto folga = global::tol_limInjReat; 

	for (int i = 0; i < sistema->nPV; i++) { 
		
		float_type aux = (iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]);
		if (aux < (sistema->limQinf[i] - folga)) {
			return false;
		}
		else if (aux > (sistema->limQsup[i] + folga)) {
			return false;
		}
	}
	return true;
}

bool hasTrue(bool* arr, int size) {
	for (int i = 0; i < size; i++) { 
		if (arr[i]) {
			return true;
		}
	}
	return false;
}

float_type accumulate(float_type* arr, int size) {
	float_type ans = 0;
	for (int i = 0; i < size; i++) { 
		ans += abs(arr[i]);
	}
	return ans;
}

void nR(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo){
	initLim(sistema, barra, ramo, iterativo);
	if (global::verbose_mode && global::lim_inj_reat) {
		printAlLim(*sistema, *barra, *ramo, *iterativo);
	}

	{ BENCHMARK_PROCESSOITERATIVO
		for (iterativo->iteracao = 0; iterativo->iteracao < global::no_max_iter; iterativo->iteracao++) {
			if (global::verbose_mode) {
				printf("\n****************************************************************** ITERACAO %d ******************************************************************\n", iterativo->iteracao);
			}

			{ BENCHMARK_CALCPQ
				switch (global::metodo)
				{
				case denso:
				case esparsoSimples:
					dnCalculePQ(sistema, barra, ramo, iterativo);
					break;
				case esparso:
					spCalcP(*sistema, *barra, *ramo, *iterativo);
					spCalcQ(*sistema, *barra, *ramo, *iterativo);
					break;
				default:
					printf("\n\n[ERRO] calcPQ: metodo inválido!\n\n");
					break;
				}
			}

			if ((iterativo->iteracao >= 1) && global::lim_inj_reat) {
				if (maxi(&(iterativo->gLim[0]), iterativo->ngLim) < 0.1) { 
					evalLimInjReat(sistema, barra, ramo, iterativo);
					
				}
				else if (hasTrue(iterativo->limQ, sistema->nPV) ) {
					
					undo_LimInjReat(sistema, barra, ramo, iterativo);
					
				}
			}

			calcResLim(sistema, barra, ramo, iterativo); 

			if (global::verbose_mode) {
				printf("Pcalc^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->Pcalc, sistema->nB, 4);
				printf("\nQcalc^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->Qcalc, sistema->nB, 4);
				printf("\ndeltaP^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->gLim, sistema->nB - 1, 4);
				printf("\ndeltaQ^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->gLim + sistema->nB - 1, iterativo->nPQlim, 4);

				printf("\nQg_barrasPV^(%d) =\n", iterativo->iteracao);
				cout << '[' << iterativo->Qcalc[IDX1F(sistema->barrasPV[0])] + barra->Qload[IDX1F(sistema->barrasPV[0])] << endl;
				for (int i = 1; i < (sistema->nPV - 1); i++) {
					cout << ' ' << iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])] << endl;
				}
				cout << ' ' << iterativo->Qcalc[IDX1F(sistema->barrasPV[sistema->nPV - 1])] + barra->Qload[IDX1F(sistema->barrasPV[sistema->nPV - 1])] << ']' << endl;

				cout << "\nlimQ:" << endl;
				showVec(iterativo->limQ, sistema->nPV, 4);
				cout << "\nbarrasPVlim:" << endl;
				showVec(iterativo->barrasPVlim, iterativo->nPVlim, 4);
				cout << "\nbarrasPQlim:" << endl;
				showVec(iterativo->barrasPQlim, iterativo->nPQlim, 4);
				cout << "\nlimQinf:" << endl;
				showVec(sistema->limQinf, sistema->nPV, 4);
				cout << "\nlimQsup:" << endl;
				showVec(sistema->limQsup, sistema->nPV, 4);
				cout << "\nQliqLim:" << endl;
				showVec(iterativo->QliqLim, sistema->nB, 4);
			}

			float_type maxG = maxi(iterativo->gLim, iterativo->ngLim); 
			if ((global::verbose_mode || global::output_processo_iterativo) && !global::laconic_mode) {
				printf("max(gLim) = %e; sum(gLim) = %e\n", maxG, accumulate(iterativo->gLim, iterativo->ngLim));
			}
			if (maxG <= global::tol) {
				if (global::lim_inj_reat && !simpleChkLimQ(sistema, barra, ramo, iterativo)) {
					
					evalLimInjReat(sistema, barra, ramo, iterativo);
					
				}
				else {
					if (global::verbose_mode) {
						printf("\n\nbreak!\nmax(g) = %e\n\n", maxi(iterativo->gLim, iterativo->ngLim));
					}
					return; 
				}

			}
			if (global::verbose_mode) {
				printf("     ...     continua:\n\n");
			}

			switch (global::metodo) {
				case esparso:
	#if JACOBIANO_ESPARSO_STENCIL
				{ BENCHMARK_JACOBIANOSTENCIL_FILL
					
					fillJstencil0based(sistema, barra, ramo, iterativo);
				}
	#else
				{ BENCHMARK_JACOBIANO
					
					spCalcJ_eficiente(sistema, barra, ramo, iterativo);
					
				}
	#endif
					break;
				default:
				{ BENCHMARK_JACOBIANO
					calcJacLim(sistema, barra, ramo, iterativo);
				}
					break;
			}

			{ BENCHMARK_SISTEMALINEAR
				if (solver(sistema, barra, ramo, iterativo)) {
					iterativo->iteracao = global::no_max_iter; 
					break; 
				}
			}

			if(global::verbose_mode){
				printf("\ndeltaX =\n");
				showVec(iterativo->gLim, sistema->nB - 1 + iterativo->nPQlim, 5);
			}
		
			attVOlim(sistema, barra, iterativo);

		}
	}

	if (maxi(&(iterativo->gLim[0]), iterativo->nPQlim + iterativo->nPVlim + iterativo->nPQlim) > global::tol) {
		
			printf("\n\nO processo iterativo nao convergiu!\nmax(g) = %e\n\n", maxi(iterativo->gLim, iterativo->nPQlim + iterativo->nPVlim + iterativo->nPQlim));
		
	}

}