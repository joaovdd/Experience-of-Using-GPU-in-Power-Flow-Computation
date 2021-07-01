#include "spCalcPQliq.h"

// Calcula valores de P
void spCalcP(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo){
	//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
#pragma omp parallel for if (global::openmp)
	for (unsigned int idx = 1; idx <= sistema.nB; idx++){
		//calculo de P_idx	
		float_type aux = 0;
		float_type acc = 0;

		// idx itera sistema.cooRowIndY
		// i <- posições de spYval com os eltos da idx-ésima linha de Y
		for(int i = sistema.csrRowPtrY[IDX1F(idx)]; i < sistema.csrRowPtrY[IDX1F(idx) + 1]; i++){
			//aux = phif(idx, sistema.csrColIndY[i] + 1, &sistema, &ramo); // defasagem do transformador
			aux = barra.theta[sistema.cooRowIndY[i]] - barra.theta[sistema.csrColIndY[i]]; // theta_k para[i]
			aux = real(sistema.spYvalE[i]) * cos(aux) + imag(sistema.spYvalE[i]) * sin(aux);
			aux *= barra.V[sistema.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barra.V[IDX1F(idx)];
		iterativo.Pcalc[IDX1F(idx)] = acc; // salva o valor
	}
}

void spCalcQ(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo){
	//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
#pragma omp parallel for if (global::openmp)
	for (unsigned int idx = 1; idx <= sistema.nB; idx++){
 		//calculo de P_idx	
 		float_type aux = 0;
		float_type acc = 0;

		// idx itera sistema.cooRowIndY
		// i é a coluna // i   itera sistema.csrColIndY
		for(int i = sistema.csrRowPtrY[IDX1F(idx)]; i < sistema.csrRowPtrY[IDX1F(idx) + 1]; i++){
			// aux = phif(idx, sistema.csrColIndY[i] + 1, &sistema, &ramo); // defasagem do transformador
			aux = barra.theta[sistema.cooRowIndY[i]] - barra.theta[sistema.csrColIndY[i]]; // theta_k para[i]
			aux = real(sistema.spYvalE[i]) * sin(aux) - imag(sistema.spYvalE[i]) * cos(aux);
			aux *= barra.V[sistema.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
 		acc *= barra.V[IDX1F(idx)];
		iterativo.Qcalc[IDX1F(idx)] = acc; // salva o valor
 	}
}

