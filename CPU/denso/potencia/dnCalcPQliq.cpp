#include "dnCalcPQliq.h"

float_type calcP(const int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (int i = 0; i < sistema->nL; i++) { 
		float_type aux = 0;
		if (ramo->de[i] == k) { 
			aux = ramo->phi[i]; 
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->para[i])]; 
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * sin(aux);
			aux *= barra->V[IDX1F(ramo->para[i])];

			acc += aux;
			continue;
		}
		else if (ramo->para[i] == k) { 
			aux = -ramo->phi[i]; 
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->de[i])]; 
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * sin(aux);
			aux *= barra->V[IDX1F(ramo->de[i])];

			acc += aux;
			continue;
		}
	}

	float_type aux = 0;

	aux = _cuReal(sistema->Y[IDX2F(k, k, sistema->nB)]);
	aux = barra->V[IDX1F(k)] * aux;
	acc += aux;

	return (barra->V[IDX1F(k)] * acc);
}

float_type calcP2(const int k, sistemaType* sistema, barraType* barra, ramoType* ramo){
	float_type acc = 0;
	for (int i = 1; i <= sistema->nB; i++) { 
		float_type aux = 0;

		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(i)]; 
		aux = _cuReal(sistema->Y[IDX2F(k, i, sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, i, sistema->nB)]) * sin(aux);
		aux *= barra->V[IDX1F(i)];

		acc += aux;
	}
	return (barra->V[IDX1F(k)] * acc);
}

float_type calcQ(const int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (int i = 0; i < sistema->nL; i++) { 
		float_type aux = 0;
		if (ramo->de[i] == k) { 
			aux = ramo->phi[i]; 
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->para[i])]; 
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * cos(aux);
			aux *= barra->V[IDX1F(ramo->para[i])];

			acc += aux;
			continue;
		}
		else if (ramo->para[i] == k) { 
			aux = -ramo->phi[i]; 
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->de[i])]; 
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * cos(aux);
			aux *= barra->V[IDX1F(ramo->de[i])];

			acc += aux;
			continue;
		}
	}

	float_type aux = 0;

	aux = -_cuImag(sistema->Y[IDX2F(k, k, sistema->nB)]);
	aux *= barra->V[IDX1F(k)];
	acc += aux;

	return barra->V[IDX1F(k)] * acc;
}

float_type calcQ2(const int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (int i = 1; i <= sistema->nB; i++) { 
		float_type aux = 0;

		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(i)]; 
		aux = _cuReal(sistema->Y[IDX2F(k, i, sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, i, sistema->nB)]) * cos(aux);
		aux *= barra->V[IDX1F(i)];

		acc += aux;

	}
	return barra->V[IDX1F(k)] * acc;
}

void dnCalculePQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo){
	#pragma omp parallel for if (global::openmp)
		for(int i = 1; i <= sistema->nB; i++){
			iterativo->Pcalc[IDX1F(i)] = calcP2(i, sistema, barra, ramo);
		}

	#pragma omp parallel for if (global::openmp)
		for(int i = 1; i <= sistema->nB; i++){
    		iterativo->Qcalc[IDX1F(i)] = calcQ(i, sistema, barra, ramo);
		}

}

