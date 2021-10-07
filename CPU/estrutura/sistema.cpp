#include "sistema.h"

void initSistema(sistemaType &sistema){
	sistema.barraVO = 0;

	#ifndef FLUMEN_GPU 
		if (global::metodo == denso || global::metodo == esparsoSimples) {
			sistema.Y = (complex_type *)malloc(sistema.nB*sistema.nB * sizeof(complex_type));

			for (int i = 0; i < sistema.nB*sistema.nB; i++){
				sistema.Y[i] = _mkComplex(0., 0.);
			}
		}
	#else
		
		if (global::metodo == denso || global::metodo == hibridoB) {
			sistema.Y = (complex_type *)malloc(sistema.nB*sistema.nB * sizeof(complex_type));

			for (int i = 0; i < sistema.nB*sistema.nB; i++){
				sistema.Y[i] = _mkComplex
			}
		}
	#endif

	sistema.barrasPV = (int *)malloc(sistema.nPV * sizeof(int));
	sistema.barrasPQ = (int *)malloc(sistema.nPQ * sizeof(int));

	sistema.limQinf = (float_type*)malloc(sistema.nPV * sizeof(float_type));
	sistema.limQsup = (float_type*)malloc(sistema.nPV * sizeof(float_type));

	sistema.VfixadaPV = (float_type*)malloc(sistema.nPV * sizeof(float_type));

}

void finSistema(sistemaType &sistema) {
	if (sistema.Y == nullptr) { free(sistema.Y); }

	if (sistema.barrasPV == nullptr) { free(sistema.barrasPV); }
	if (sistema.barrasPQ == nullptr) { free(sistema.barrasPQ); }
	if (sistema.limQinf == nullptr) { free(sistema.limQinf); }
	if (sistema.limQsup == nullptr) { free(sistema.limQsup); }
	if (sistema.VfixadaPV == nullptr) { free(sistema.VfixadaPV); }

	if (sistema.spY == nullptr) { delete(sistema.spY); }
	if (sistema.spYval == nullptr) { free(sistema.spYval); }
	if (sistema.csrRowPtrY == nullptr) { free(sistema.csrRowPtrY); }
	if (sistema.csrColIndY == nullptr) { free(sistema.csrColIndY); }
}