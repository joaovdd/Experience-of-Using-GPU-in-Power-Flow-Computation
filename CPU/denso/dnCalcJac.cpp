#include "dnCalcJac.h"

// Jacobiano com lim inj reat

void calcHlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	unsigned int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (unsigned int idx = 1; idx <= sistema->nB - 1; idx++) {
			//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
			bool offm = (idx >= sistema->barraVO); // barra idx, posicao idx+offswing
			idx += offm; // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2

			// idx percorre colunas de H diretamente
			// idx=1:nPQ
			// i percorre linhas de H diretamente
			// i=1:nPQ+nPV

			for (unsigned int i = 1; i <= sistema->nB - 1; i++) {
				bool offk = (i >= sistema->barraVO);
				i += offk; // pula barra swing
				if (i != idx) {
					//fórmula pra km
					//float_type aux = dnPhif(i, idx, sistema, ramo);
					float_type aux = barra->theta[IDX1F(i)] - barra->theta[IDX1F(idx)];
					// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
					iterativo->Jlim[IDX2F(i - offk, idx - offm, szJ)] = barra->V[IDX1F(i)] * barra->V[IDX1F(idx)] * (sistema->Y[IDX2F(i, idx, sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(i, idx, sistema->nB)].y * cos(aux));
				}
				i -= offk;
			}
			// Hkk = -Qk - V^2_k*Bkk
			iterativo->Jlim[IDX2F(idx - offm, idx - offm, szJ)] = -iterativo->Qcalc[IDX1F(idx)] - barra->V[IDX1F(idx)] * barra->V[IDX1F(idx)] * sistema->Y[IDX2F(idx, idx, sistema->nB)].y;

			idx -= offm;
		}
}

void calcLlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	unsigned int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (unsigned int idx = 1; idx <= iterativo->nPQlim; idx++) {
			//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1
			unsigned int offset = iterativo->nPVlim + iterativo->nPQlim; // deslocamento para armazenar elementos na matriz J

			// i   é iterador do vetor barrasPQ, percorre linhas  de L
			// idx é iterador do vetor barrasPQ, percorre colunas de L

			//bool offm = (barrasPQ[IDX1F(idx)] > swing); // se idx>swing => armazenara elto uma posição a menos em J

			for (unsigned int i = 1; i <= iterativo->nPQlim; i++) {
				if (iterativo->barrasPQlim[IDX1F(i)] != iterativo->barrasPQlim[IDX1F(idx)]) {
					//fórmula pra km
					//float_type aux = dnPhif(iterativo->barrasPQlim[IDX1F(i)], iterativo->barrasPQlim[IDX1F(idx)], sistema, ramo);
					float_type aux = barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] - barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];

					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
					iterativo->Jlim[IDX2F(i + offset, idx + offset, szJ)] = barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] * (sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y * cos(aux));
				}
			}
			// Lkk= (Qk - V^2_k*Bkk)/Vk
			iterativo->Jlim[IDX2F(idx + offset, idx + offset, szJ)] = (iterativo->Qcalc[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] - barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(idx)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y) / barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];
		}
}

void calcMlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	unsigned int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (unsigned int idx = 1; idx <= (sistema->nB - 1); idx++) {
			//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1 // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2
			bool offm = (idx >= sistema->barraVO); // barra idx, posicao idx+offswing
			idx += offm;

			unsigned int offset = sistema->nB - 1; // =nPV+nPQ deslocamento para armazenar elementos na matriz J

			// i é iterador do vetor barrasPQ, percorre linhas de M
			// i=1:nPQ
			// idx percorre colunas de M diretamente
			// idx=1:nB - swing

			for (unsigned int i = 1; i <= iterativo->nPQlim; i++) {
				if (iterativo->barrasPQlim[IDX1F(i)] != idx) {
					//fórmula pra km
					//float_type aux = dnPhif(iterativo->barrasPQlim[IDX1F(i)], idx, sistema, ramo);
					float_type aux = barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] - barra->theta[IDX1F(idx)];

					// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
					iterativo->Jlim[IDX2F(i + offset, idx - offm, szJ)] = -barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(i)])] * barra->V[IDX1F(idx)] * (sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], idx, sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(i)], idx, sistema->nB)].y * sin(aux));
				}
				else {
					// Mkk= (Pk - V^2_k*Gkk)
					iterativo->Jlim[IDX2F(i + offset, idx - offm, szJ)] = iterativo->Pcalc[IDX1F(idx)] - barra->V[IDX1F(idx)] * barra->V[IDX1F(idx)] * sistema->Y[IDX2F(idx, idx, sistema->nB)].x;
				}
			}
			idx -= offm;
		}
}

void calcNlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	unsigned int szJ = iterativo->nPVlim + iterativo->nPQlim + iterativo->nPQlim;
	#pragma omp parallel for if (global::openmp)
		for (unsigned int idx = 1; idx <= iterativo->nPQlim; idx++) {
			//unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x + 1; //começa de 1 // if idx >= swing : idx = idx = threadIdx.x + blockDim.x * blockIdx.x + 2

			unsigned int offset = sistema->nB - 1; // =nPV+nPQ deslocamento para armazenar elementos na matriz J (colunas)

			// idx é iterador do vetor barrasPQ, percorre colunas de N
			// idx=1:nPQ
			// i percorre linhas de N diretamente
			// i=1:nPQ+nPV

			for (unsigned int i = 1; i <= sistema->nPQ + sistema->nPV; i++) {
				bool offk = (i >= sistema->barraVO);
				i += offk;
				if (i != iterativo->barrasPQlim[IDX1F(idx)]) {
					//fórmula pra km=iidx
					//float_type aux = dnPhif(i, iterativo->barrasPQlim[IDX1F(idx)], sistema, ramo);
					float_type aux = barra->theta[IDX1F(i)] - barra->theta[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];

					// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
					iterativo->Jlim[IDX2F(i - offk, idx + offset, szJ)] = barra->V[IDX1F(i)] * (sistema->Y[IDX2F(i, iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(i, iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].y * sin(aux));
				}
				else {
					// Nkk= (Qk - V^2_k*Bkk)/Vkk
					iterativo->Jlim[IDX2F(i - offk, idx + offset, szJ)] = (iterativo->Pcalc[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] + barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])] * sistema->Y[IDX2F(iterativo->barrasPQlim[IDX1F(idx)], iterativo->barrasPQlim[IDX1F(idx)], sistema->nB)].x) / barra->V[IDX1F(iterativo->barrasPQlim[IDX1F(idx)])];
				}
				i -= offk;
			}
		}
}

void calcJacLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int sz = (sistema->nB - 1 + sistema->nPQ) * (sistema->nB - 1 + sistema->nPQ);

	//printf("Jacobiano =\n");
	//showMat(iterativo->J, sistema->nB - 1 + sistema->nPQ);
	calcHlimf(sistema, barra, ramo, iterativo);

	//printf("Jacobiano =\n");
	//showMat(iterativo->J, sistema->nB - 1 + sistema->nPQ);

	calcLlimf(sistema, barra, ramo, iterativo);

	//printf("Jacobiano =\n");
	//showMat(iterativo->J, sistema->nB - 1 + sistema->nPQ);

	calcMlimf(sistema, barra, ramo, iterativo);

	//printf("Jacobiano =\n");
	//showMat(iterativo->J, sistema->nB - 1 + sistema->nPQ);

	calcNlimf(sistema, barra, ramo, iterativo);

	if (global::verbose_mode) {
		printf("Jacobiano =\n");
		showMat(iterativo->Jlim, sistema->nB - 1 + iterativo->nPQlim);
	}
}