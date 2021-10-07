#include "dnCalcJac.h"

void calcHlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (int idx = 1; idx <= sistema->nB - 1; idx++) {
			
			bool offm = (idx >= sistema->barraVO); 
			idx += offm; 

			for (int i = 1; i <= sistema->nB - 1; i++) {
				bool offk = (i >= sistema->barraVO);
				i += offk; 
				if (i != idx) {
					
					float_type aux = barra->theta[IDX1F(i)] - barra->theta[IDX1F(idx)];
					
					iterativo->Jlim[IDX2F(i - offk, idx - offm, szJ)] = barra->V[IDX1F(i)] * barra->V[IDX1F(idx)] * (sistema->Y[IDX2F(i, idx, sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(i, idx, sistema->nB)].y * cos(aux));
				}
				i -= offk;
			}
			
			iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -iterativo->Qcalc[IDX1F(idx)] - barra->V[IDX1F(idx)] * barra->V[IDX1F(idx)] * sistema->Y[IDX2F(idx, idx, sistema->nB)].y;

			idx -= offm;
		}
}

void calcLlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (int idx = 1; idx <= iterativo->nPQlim; idx++) {
			
			int offset = iterativo->nPVlim + iterativo->nPQlim; 

			for (int i = 1; i <= iterativo->nPQlim; i++) {
				if (iterativo->barrasPQlim[IDX1F(i)] != iterativo->barrasPQlim[IDX1F(idx)]) {
					
					float_type aux = barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] - barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];

					iterativo->Jlim[IDX2F(i + offset, idx + offset, szJ)] = barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] * (sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y * cos(aux));
				}
			}
			
			iterativo->Jlim[IDX2F(idx + offset, idx + offset, szJ)] = (iterativo->Qcalc[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] - barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(idx)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y) / barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];
		}
}

void calcMlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (int idx = 1; idx <= (sistema->nB - 1); idx++) {
			
			bool offm = (idx >= sistema->barraVO); 
			idx += offm;

			int offset = sistema->nB - 1; 

			for (int i = 1; i <= iterativo->nPQlim; i++) {
				if (iterativo->barrasPQlim[IDX1F(i)] != idx) {
					
					float_type aux = barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] - barra->theta[IDX1F(idx)];

					iterativo->Jlim[IDX2F(i + offset, idx - offm, szJ)] = -barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] * barra->V[IDX1F(idx)] * (sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], idx, sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], idx, sistema->nB)].y * sin(aux));
				}
				else {
					
					iterativo->Jlim[IDX2F(i + offset, idx - offm, szJ)] = iterativo->Pcalc[IDX1F(idx)] - barra->V[IDX1F(idx)] * barra->V[IDX1F(idx)] * sistema->Y[IDX2F(idx, idx, sistema->nB)].x;
				}
			}
			idx -= offm;
		}
}

void calcNlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (int idx = 1; idx <= iterativo->nPQlim; idx++) {
			
			int offset = sistema->nB - 1; 

			for (int i = 1; i <= sistema->nPQ + sistema->nPV; i++) {
				bool offk = (i >= sistema->barraVO);
				i += offk;
				if (i != iterativo->barrasPQlim[IDX1F(idx)]) {
					
					float_type aux = barra->theta[IDX1F(i)] - barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];

					iterativo->Jlim[IDX2F(i - offk, idx + offset, szJ)] = barra->V[IDX1F(i)] * (sistema->Y[IDX2F(i, iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(i, iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y * sin(aux));
				}
				else {
					
					iterativo->Jlim[IDX2F(i - offk, idx + offset, szJ)] = (iterativo->Pcalc[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] + barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(idx)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x) / barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];
				}
				i -= offk;
			}
		}
}

void calcJacLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	calcHlimf(sistema, barra, ramo, iterativo);

	calcLlimf(sistema, barra, ramo, iterativo);

	calcMlimf(sistema, barra, ramo, iterativo);

	calcNlimf(sistema, barra, ramo, iterativo);

	if (global::verbose_mode) {
		printf("Jacobiano =\n");
		showMat(iterativo->Jlim, sistema->nB - 1 + iterativo->nPQlim);
	}
}