#include "dnCalcPQliq.h"

// minimiza acesso a elementos
float_type calcP(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (unsigned int i = 0; i < sistema->nL; i++) { // i percorre ramos
		float_type aux = 0;
		if (ramo->de[i] == k) { // se k est� ligada � para[i] =>
			aux = ramo->phi[i]; //atan2(ramo->tap[i].y, ramo->tap[i].x); // defasagem do transformador
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->para[i])]; // theta_k para[i]
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * sin(aux);
			aux *= barra->V[IDX1F(ramo->para[i])];

			acc += aux;
			continue;
		}
		else if (ramo->para[i] == k) { // se k est� ligada � de[i] =>
			aux = -ramo->phi[i]; //-atan2(ramo->tap[i].y, ramo->tap[i].x); // defasagem do transformador da i-�sima LT
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->de[i])]; // theta_k de[i]
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * sin(aux);
			aux *= barra->V[IDX1F(ramo->de[i])];

			acc += aux;
			continue;
		}
	}

	float_type aux = 0;
	//aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(k)]; // theta_k para[i]
	//aux = _cuReal(sistema->Y[IDX2F(k, k, sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, k, sistema->nB)]) * sin(aux);
	aux = _cuReal(sistema->Y[IDX2F(k, k, sistema->nB)]);
	aux = barra->V[IDX1F(k)] * aux;
	acc += aux;

	return (barra->V[IDX1F(k)] * acc);
}

// denso
float_type calcP2(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo){
	float_type acc = 0;
	for (unsigned int i = 1; i <= sistema->nB; i++) { // i percorre colunas da k-�sima linha da Ybarra
		float_type aux = 0;

		// aux = dnPhif(k, i, sistema, ramo); // defasagem do transformador
		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(i)]; // theta_k para[i]
		aux = _cuReal(sistema->Y[IDX2F(k, i, sistema->nB)]) * cos(aux) + _cuImag(sistema->Y[IDX2F(k, i, sistema->nB)]) * sin(aux);
		aux *= barra->V[IDX1F(i)];

		acc += aux;
	}
	return (barra->V[IDX1F(k)] * acc);
}

// minimiza acesso a elementos
float_type calcQ(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (unsigned int i = 0; i < sistema->nL; i++) { // i percorre ramos
		float_type aux = 0;
		if (ramo->de[i] == k) { // se k est� ligada � para[i] =>
			aux = ramo->phi[i]; //atan2(ramo->tap[i].y, ramo->tap[i].x); // defasagem do transformador
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->para[i])]; // theta_k para[i]
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, ramo->para[i], sistema->nB)]) * cos(aux);
			aux *= barra->V[IDX1F(ramo->para[i])];

			acc += aux;
			continue;
		}
		else if (ramo->para[i] == k) { // se k est� ligada � de[i] =>
			aux = -ramo->phi[i]; //-atan2(ramo->tap[i].y, ramo->tap[i].x); // defasagem do transformador da i-�sima LT
			aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->de[i])]; // theta_k de[i]
			aux = _cuReal(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * cos(aux);
			aux *= barra->V[IDX1F(ramo->de[i])];

			acc += aux;
			continue;
		}
	}

	float_type aux = 0;
	//aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(k)]; // theta_k para[i]
	//aux = _cuReal(sistema->Y[IDX2F(k, k, sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, k, sistema->nB)]) * cos(aux);
	aux = -_cuImag(sistema->Y[IDX2F(k, k, sistema->nB)])/* * cos(aux)*/;
	aux *= barra->V[IDX1F(k)];
	acc += aux;

	return barra->V[IDX1F(k)] * acc;
}

// denso
float_type calcQ2(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type acc = 0;
	for (unsigned int i = 1; i <= sistema->nB; i++) { // i percorre colunas da k-�sima linha da Ybarra
		float_type aux = 0;

		//aux = dnPhif(k, i, sistema, ramo); // defasagem do transformador
		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(i)]; // theta_k para[i]
		aux = _cuReal(sistema->Y[IDX2F(k, i, sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, i, sistema->nB)]) * cos(aux);
		aux *= barra->V[IDX1F(i)];

		acc += aux;

		//aux = -ramo->phi[i]; //-atan2(ramo->tap[i].y, ramo->tap[i].x); // defasagem do transformador da i-�sima LT
		//aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(ramo->de[i])]; // theta_k de[i]
		//aux = _cuReal(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * sin(aux) - _cuImag(sistema->Y[IDX2F(k, ramo->de[i], sistema->nB)]) * cos(aux);
		//aux *= barra->V[IDX1F(ramo->de[i])];

		//acc += aux;
	}
	return barra->V[IDX1F(k)] * acc;
}

void dnCalculePQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo/*, const unsigned int* card*/){
	#pragma omp parallel for if (global::openmp)
		for(unsigned int i = 1; i <= sistema->nB; i++){
			iterativo->Pcalc[IDX1F(i)] = calcP2(i, sistema, barra, ramo);
		}

	//// calcula P para barras PQ e PV
	//#pragma omp parallel for if (global::openmp)
	//	for (int i = 1; i <= sistema->barraVO; i++) {
	//		iterativo->Pcalc[IDX1F(i)] = calcP2(i, sistema, barra, ramo);
	//	}
	//#pragma omp parallel for if (global::openmp)
	//	for (int i = sistema->barraVO + 1; i <= sistema->nB; i++) {
	//		iterativo->Pcalc[IDX1F(i)] = calcP2(i, sistema, barra, ramo);
	//	}

	 //Calcula os valores de Qk...
	#pragma omp parallel for if (global::openmp)
		for(unsigned int i = 1; i <= sistema->nB; i++){
    		iterativo->Qcalc[IDX1F(i)] = calcQ(i, sistema, barra, ramo);
		}

	//// calcula Q para barras PQ
	//#pragma omp parallel for if (global::openmp)
	//	for (int i = 0; i < iterativo->nPQlim; i++) {
	//		iterativo->Qcalc[IDX1F(sistema->barrasPQ[i])] = calcQ(sistema->barrasPQ[i], sistema, barra, ramo);
	//	}
}

//void dnCalculePQ_subsist2(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
//	// calcula P para barra VO
//	iterativo->Pcalc[IDX1F(sistema->barraVO)] = calcP2(sistema->barraVO, sistema, barra, ramo);
//
//	// calcula Q para barras PV
//	#pragma omp parallel for if (global::openmp)
//		for (int i = 0; i < sistema->nPV; i++) {
//			iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] = calcQ(sistema->barrasPV[i], sistema, barra, ramo);
//		}
//	// e referência
//	iterativo->Qcalc[IDX1F(sistema->barraVO)] = calcQ(sistema->barraVO, sistema, barra, ramo);
//}