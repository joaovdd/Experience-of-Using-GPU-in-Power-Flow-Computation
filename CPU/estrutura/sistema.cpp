#include "sistema.h"

void initSistema(sistemaType &sistema){
	//sistema.nPV = 0;
	//sistema.nPQ = 0;
	sistema.barraVO = 0;
	//sistema.nL = 0;
	//sistema.nB = 0;
	//sistema.baseMVA = 0;

	#ifndef FLUMEN_GPU // #if !defined FLUMEN_GPU || !defined __NVCC__ // https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#nvcc-identification-macro
		if (global::metodo == denso || global::metodo == esparsoSimples) {
			sistema.Y = (complex_type *)malloc(sistema.nB*sistema.nB * sizeof(complex_type));

			for (unsigned int i = 0; i < sistema.nB*sistema.nB; i++){
				sistema.Y[i] = _mkComplex(0., 0.);
			}
		}
	#else
		// ... para gpu
		if (global::metodo == denso || global::metodo == hibridoB) {
			sistema.Y = (complex_type *)malloc(sistema.nB*sistema.nB * sizeof(complex_type));

			for (unsigned int i = 0; i < sistema.nB*sistema.nB; i++){
				sistema.Y[i] = _mkComplex
			}
		}
	#endif

	sistema.barrasPV = (unsigned int *)malloc(sistema.nPV * sizeof(unsigned int));
	sistema.barrasPQ = (unsigned int *)malloc(sistema.nPQ * sizeof(unsigned int));

	//limite de injeção de reativos
	sistema.limQinf = (float_type*)malloc(sistema.nPV * sizeof(float_type));
	sistema.limQsup = (float_type*)malloc(sistema.nPV * sizeof(float_type));

	sistema.VfixadaPV = (float_type*)malloc(sistema.nPV * sizeof(float_type));
	//memset(sistema.VfixadaPV, 0, sistema.nPV * sizeof(float_type));
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