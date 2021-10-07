#include "iterativo.h"

void initIter(sistemaType &sistema, iterativoType &iterativo){
	iterativo.Pcalc =  (float_type *)malloc(sistema.nB * sizeof(float_type));
	iterativo.Qcalc =  (float_type *)malloc(sistema.nB * sizeof(float_type));

	iterativo.iteracao = 0;
	iterativo.noMax = global::no_max_iter;
	for (int i = 0; i < sistema.nB-1; i++){
		iterativo.Pcalc[i] = 0.;
	}
	for (int i = 0; i < sistema.nPQ; i++){
		iterativo.Qcalc[i] = 0.;
	}

	iterativo.gLim = (float_type*)malloc((sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * sizeof(float_type));

	if (global::metodo == metodoType::denso || global::metodo == metodoType::esparsoSimples) {
		iterativo.Jlim = (float_type*)malloc((sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * sizeof(float_type));
		for (int i = 0; i < (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ); i++) {
			iterativo.Jlim[i] = 0.;
		}
	}

	iterativo.limQ = (bool*)malloc(sistema.nPV * sizeof(bool));
	memset(iterativo.limQ, 0, sistema.nPV * sizeof(bool));
	iterativo.barrasPVlim = (int*)malloc(sistema.nPV * sizeof(int));
	iterativo.barrasPQlim = (int*)malloc((sistema.nPV + sistema.nPQ) * sizeof(int));
	iterativo.QliqLim = (float_type*)malloc(sistema.nB * sizeof(float_type));
}

void finIter(iterativoType &iterativo) {
	if (iterativo.Pcalc == nullptr) { free(iterativo.Pcalc); }
	if (iterativo.Qcalc == nullptr) { free(iterativo.Qcalc); }
	if (iterativo.J == nullptr) { free(iterativo.J); }
	if (iterativo.g == nullptr) { free(iterativo.g); }

	if (iterativo.gLim == nullptr) { free(iterativo.gLim);}
	if (iterativo.Jlim == nullptr) { free(iterativo.Jlim); }
	if (iterativo.limQ == nullptr) { free(iterativo.limQ); }
	if (iterativo.barrasPVlim == nullptr) { free(iterativo.barrasPVlim); }
	if (iterativo.barrasPQlim == nullptr) { free(iterativo.barrasPQlim); }
	if (iterativo.QliqLim == nullptr) { free(iterativo.QliqLim); }
}
