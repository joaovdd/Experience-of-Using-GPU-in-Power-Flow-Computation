#include "spCalcPQliq.h"

void spCalcP(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo){
#pragma omp parallel for if (global::openmp)
	for (int idx = 1; idx <= sistema.nB; idx++){
		
		float_type aux = 0;
		float_type acc = 0;

		
		
		for(int i = sistema.csrRowPtrY[IDX1F(idx)]; i < sistema.csrRowPtrY[IDX1F(idx) + 1]; i++){
			
			aux = barra.theta[sistema.cooRowIndY[i]] - barra.theta[sistema.csrColIndY[i]]; 
			aux = real(sistema.spYvalE[i]) * cos(aux) + imag(sistema.spYvalE[i]) * sin(aux);
			aux *= barra.V[sistema.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barra.V[IDX1F(idx)];
		iterativo.Pcalc[IDX1F(idx)] = acc; 
	}
}

void spCalcQ(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo){
#pragma omp parallel for if (global::openmp)
	for (int idx = 1; idx <= sistema.nB; idx++){
 		
 		float_type aux = 0;
		float_type acc = 0;

		
		
		for(int i = sistema.csrRowPtrY[IDX1F(idx)]; i < sistema.csrRowPtrY[IDX1F(idx) + 1]; i++){
			
			aux = barra.theta[sistema.cooRowIndY[i]] - barra.theta[sistema.csrColIndY[i]]; 
			aux = real(sistema.spYvalE[i]) * sin(aux) - imag(sistema.spYvalE[i]) * cos(aux);
			aux *= barra.V[sistema.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
 		acc *= barra.V[IDX1F(idx)];
		iterativo.Qcalc[IDX1F(idx)] = acc; 
 	}
}

