#pragma once
//#include "../sistInfo.h"
//#include "../sistema.h"
//#include "../dim.h"
//#include "../PcalcQcalc.h"
//#include "../jacobiano.h"

#include <algorithm>
#include <vector>

//#include "lapacke.h"

#include <stdlib.h>
#include <stdio.h>

#include "Eigen/Sparse"
using namespace std;
using namespace Eigen;

std::pair<int, int> HtoJ(int row, int col, sistema* sistema) {
	//std::pair<int, int> ans = 
	return { row - (row >= sistema->barraVO), col - (col >= sistema->barraVO) };
}

int HtoJrow(int row, sistema* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int HtoJcol(int col, sistema* sistema) {
	return col+1 - (col+1 >= sistema->barraVO);
}

std::pair<int, int> LtoJ(int row, int col, sistema* sistema) {
	//std::pair<int, int> ans = 
	return { row + sistema->nPV + sistema->nPQ, col + sistema->nPV + sistema->nPQ };
}

int LtoJrow(int row, sistema* sistema) {
	return row + sistema->nB - 1;
}

int LtoJcol(int col, sistema* sistema) {
	return col + sistema->nB - 1;
}

std::pair<int, int> MtoJ(int row, int col, sistema* sistema) {
	//std::pair<int, int> ans = 
	return { row + sistema->nPV + sistema->nPQ, col - (col >= sistema->barraVO) };
}

int MtoJrow(int row, sistema* sistema) {
	return row + sistema->nB - 1;
}

int MtoJcol(int col, sistema* sistema) {
	return col - (col >= sistema->barraVO);
}

std::pair<int, int> NtoJ(int row, int col, sistema* sistema) {
	//std::pair<int, int> ans = 
	return { row - (row >= sistema->barraVO), col + sistema->nPV + sistema->nPQ };
}

int NtoJrow(int row, sistema* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int NtoJcol(int col, sistema* sistema) {
	return col+1 + sistema->nB - 1;
}

// agora com HLMNpos
void Jstencil0based(sistema* sistema, h_sparse* h_sparse, barra* barra/*, ramo* ramo*/, iterativo* iterativo) {

	int nnzJ = 0/*, cntH = 0, cntL = 0, cntM = 0, cntN = 0*/;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(H);
			h_sparse->Hpos.push_back(nnzJ);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			h_sparse->spJsubmatType.push_back(N);
			h_sparse->Npos.push_back(nnzJ);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(H);
			h_sparse->Hpos.push_back(nnzJ);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			h_sparse->spJsubmatType.push_back(N);
			h_sparse->Npos.push_back(nnzJ);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
	}

	// M L]

	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(M);
			h_sparse->Mpos.push_back(nnzJ);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { // se coluna for PQ
				// linha de L
				h_sparse->spJsubmatType.push_back(L);
				h_sparse->Lpos.push_back(nnzJ);
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = pntRowPQ - iterativo->barrasPQlim;

				h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);
				nnzJ++;
			}
		}
	}
	h_sparse->nnzJ = nnzJ;
	h_sparse->spJval.reserve(nnzJ);
	h_sparse->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha /*+1...*/) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

// coloca primeiro os eltos da diagonal principal para minimizar a divergência
void Jstencil(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {

	int nnzJ = 0/*, cntH = 0, cntL = 0, cntM = 0, cntN = 0*/;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);
			
			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);
			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
	}

	// M L]

	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			h_sparse->spJsubmatType.push_back(M);

			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Mpos[nDiagM] = nnzJ;
				nDiagM++;
			}
			else {
				h_sparse->Mpos.push_back(nnzJ);
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { // se coluna for PQ
				// linha de L
				h_sparse->spJsubmatType.push_back(L);
				h_sparse->Lpos.push_back(nnzJ);
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = pntRowPQ - iterativo->barrasPQlim;

				h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
					h_sparse->Lpos[nDiagL] = nnzJ;
					nDiagL++;
				}
				else {
					h_sparse->Lpos.push_back(nnzJ);
				}

				nnzJ++;
			}
		}
	}
	h_sparse->nnzJ = nnzJ;
	h_sparse->spJval.reserve(nnzJ);
	h_sparse->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha /*+1...*/) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

void fillJstencil0based(sistema* sistema, h_sparse* h_sparse, barra* barra, ramo* ramo, iterativo* iterativo) {
	//h_sparse->spJval.clear();
	//h_sparse->spJval.resize(h_sparse->nnzJ);
	#pragma omp parallel for if (openmp)
		for (int i = 0; i < h_sparse->nnzJ; i++) {
			switch (h_sparse->spJsubmatType[i]) {
			case H:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					//fórmula pra km
					// float_type aux = phif(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema, ramo);
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];
					// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
					// Y denso
					// h_sparse->spJval[i] = barra->V[IDX1F(h_sparse->cooRowIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
					// Y esparso
					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					// Hkk = -Qk - V^2_k*Bkk
					// Y denso
					// h_sparse->spJval[i] = -iterativo->Qcalc[IDX1F(h_sparse->cooColIndSubMatJ[i])] - barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y;
					// Y esparso
					h_sparse->spJval[i] = -iterativo->Qcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag();
				}
				break;
			case L:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					//fórmula pra km
					//float_type aux = phif(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema, ramo);
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
					// Y denso
					// h_sparse->spJval[i] = barra->V[IDX1F(h_sparse->cooRowIndSubMatJ[i])] * (sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
					// Y esparso
					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					// Lkk= (Qk - V^2_k*Bkk)/Vk
					// Y denso
					// h_sparse->spJval[i] = (iterativo->Qcalc[IDX1F(h_sparse->cooColIndSubMatJ[i])] - barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y) / barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])];
					// Y esparso
					h_sparse->spJval[i] = (iterativo->Qcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag()) / barra->V[h_sparse->cooColIndSubMatJ[i]];
				}
				break;
			case M:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					//fórmula pra km
					// float_type aux = phif(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema, ramo);
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
					// Y denso
					// h_sparse->spJval[i] = -barra->V[IDX1F(h_sparse->cooRowIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
					// Y esparso
					h_sparse->spJval[i] = -barra->V[h_sparse->cooRowIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					// Mkk= (Pk - V^2_k*Gkk)
					// Y denso
					// h_sparse->spJval[i] = iterativo->Pcalc[IDX1F(h_sparse->cooColIndSubMatJ[i])] - barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x;
					// Y esparso
					h_sparse->spJval[i] = iterativo->Pcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real();
				}
				break;
			case N:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					//fórmula pra km=iidx
					//float_type aux = phif(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema, ramo);
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
					// Y denso
					// h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * (sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
					// Y esparso
					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					// Nkk= (Qk - V^2_k*Bkk)/Vkk
					// Y denso
					// h_sparse->spJval[i] = (iterativo->Pcalc[IDX1F(h_sparse->cooColIndSubMatJ[i])] + barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i], sistema->nB)].x) / barra->V[IDX1F(h_sparse->cooColIndSubMatJ[i])];
					// Y esparso
					h_sparse->spJval[i] = (iterativo->Pcalc[h_sparse->cooColIndSubMatJ[i]] + barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real()) / barra->V[h_sparse->cooColIndSubMatJ[i]];
				}
				break;
			default:
				printf("ERRO [fillJstencil] tipo inválido de H");
				break;
			}
		}
}

// Rotinas para reordenação do cálculo do jacobiano

void reordJac(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int szHdp = 0, szLdp = 0, szMdp = 0, szNdp = 0, aux;

	// H
	for (int i = 0; i < h_sparse->Hpos.size(); i++) {
		// se elto está na diagonal principal da submatriz, então colóca-o no início da
		// fila para evitar divergência do kernel
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Hpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Hpos[i]]) {
			aux = h_sparse->Hpos[i];
			// right shift
			for (int j = i - 1; j >= szHdp; j--) {
				h_sparse->Hpos[j + 1] = h_sparse->Hpos[j];
			}
			h_sparse->Hpos[szHdp] = aux;
			szHdp++;
		}
	}
	// L
	for (int i = 0; i < h_sparse->Lpos.size(); i++) {
		// se elto está na diagonal principal da submatriz, então colóca-o no início da
		// fila para evitar divergência do kernel
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Lpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Lpos[i]]) {
			aux = h_sparse->Lpos[i];
			// right shift
			for (int j = i - 1; j >= szLdp; j--) {
				h_sparse->Lpos[j + 1] = h_sparse->Lpos[j];
			}
			h_sparse->Lpos[szLdp] = aux;
			szLdp++;
		}
	}
	// M
	for (int i = 0; i < h_sparse->Mpos.size(); i++) {
		// se elto está na diagonal principal da submatriz, então colóca-o no início da
		// fila para evitar divergência do kernel
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Mpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Mpos[i]]) {
			aux = h_sparse->Mpos[i];
			// right shift
			for (int j = i - 1; j >= szMdp; j--) {
				h_sparse->Mpos[j + 1] = h_sparse->Mpos[j];
			}
			h_sparse->Mpos[szMdp] = aux;
			szMdp++;
		}
	}
	// N
	for (int i = 0; i < h_sparse->Npos.size(); i++) {
		// se elto está na diagonal principal da submatriz, então colóca-o no início da
		// fila para evitar divergência do kernel
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Npos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Npos[i]]) {
			aux = h_sparse->Npos[i];
			// right shift
			for (int j = i - 1; j >= szNdp; j--) {
				h_sparse->Npos[j + 1] = h_sparse->Npos[j];
			}
			h_sparse->Npos[szNdp] = aux;
			szNdp++;
		}
	}
}

int GLOBALauxIsPQ = 0;
bool isPQ(int barra, sistema* sistema, iterativo* iterativo) {
	while (iterativo->barrasPQlim[GLOBALauxIsPQ] < barra) {
		GLOBALauxIsPQ++;
		if (GLOBALauxIsPQ == sistema->nPQ) {
			return false;
		}
	}
	return iterativo->barrasPQlim[GLOBALauxIsPQ] == barra;
}

// coloca primeiro os eltos da diagonal principal para minimizar a divergência
void Jstencil_eficiente(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {

	int nnzJ = 0/*, cntH = 0, cntL = 0, cntM = 0, cntN = 0*/;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	// sistema->cooRowIndY[j] ==> i/**/

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
			//	continue;
			//}
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
		GLOBALauxIsPQ = 0;
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);
			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
			//	continue;
			//}
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
		GLOBALauxIsPQ = 0;
	}

	// M L]

	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(M);

			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
				h_sparse->Mpos[nDiagM] = nnzJ;
				nDiagM++;
			}
			else {
				h_sparse->Mpos.push_back(nnzJ);
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			//int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) /*pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim*/) { // se coluna for PQ
				// linha de L
				//h_sparse->spJsubmatType.push_back(L);
				//h_sparse->Lpos.push_back(nnzJ);
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = GLOBALauxIsPQ; //pntRowPQ - iterativo->barrasPQlim;

				h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
					h_sparse->Lpos[nDiagL] = nnzJ;
					nDiagL++;
				}
				else {
					h_sparse->Lpos.push_back(nnzJ);
				}

				nnzJ++;
			}
			GLOBALauxIsPQ = 0;
		}
	}
	h_sparse->nnzJ = nnzJ;
	h_sparse->spJval.reserve(nnzJ);
	h_sparse->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha /*+1...*/) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

// coloca primeiro os eltos da diagonal principal para minimizar a divergência
// void Jstencil_eficienteTeste(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {

// 	int nnzJ = 0/*, cntH = 0, cntL = 0, cntM = 0, cntN = 0*/;

// 	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
// 	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
// 	// loops for não são paralelizáveis (dependem do valor de nnzJ)

// 	h_sparse->Hpos.resize(sistema->nB - 1);
// 	int nDiagH = 0;
// 	h_sparse->Lpos.resize(sistema->nPQ);
// 	int nDiagL = 0;
// 	h_sparse->Mpos.resize(sistema->nPQ);
// 	int nDiagM = 0;
// 	h_sparse->Npos.resize(sistema->nPQ);
// 	int nDiagN = 0;

// 	int nnzH = 0, nnzL = 0, nnzM = 0, nnzN = 0;

// 	h_sparse->cooColIndSubMatJ.resize(137197);
// 	h_sparse->cooRowIndSubMatJ.resize(137197);

// 	h_sparse->Hpos.resize(70000);
// 	h_sparse->Lpos.resize(70000);
// 	h_sparse->Mpos.resize(70000);
// 	h_sparse->Npos.resize(70000);



// 	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			//h_sparse->spJsubmatType.push_back(H);
// 			h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);
// 			h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);

// 			//h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
// 			//h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

// 			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 				h_sparse->Hpos[nDiagH] = nnzJ;
// 				nDiagH++;
// 			}
// 			else {
// 				h_sparse->Hpos[sistema->nB - 1 + nnzH] = (nnzJ);
// 				nnzH++;
// 			}

// 			nnzJ++;
// 		}
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
// 			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
// 			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
// 			//	continue;
// 			//}
// 			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
// 				continue;
// 			}
// 			//h_sparse->spJsubmatType.push_back(N);
// 			h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);
// 			h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);

// 			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

// 			//h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
// 			//h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

// 			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 				h_sparse->Npos[nDiagN] = nnzJ;
// 				nDiagN++;
// 			}
// 			else {
// 				h_sparse->Npos[sistema->nPQ + nnzN] = (nnzJ);
// 				nnzN++;
// 			}

// 			nnzJ++;
// 		}
// 		GLOBALauxIsPQ = 0;
// 	}
// 	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			//h_sparse->spJsubmatType.push_back(H);
// 			h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);
// 			h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);

// 			//h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
// 			//h_sparse->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);
// 			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 				h_sparse->Hpos[nDiagH] = nnzJ;
// 				nDiagH++;
// 			}
// 			else {
// 				h_sparse->Hpos[sistema->nB - 1 + nnzH] = (nnzJ);
// 				nnzH++;
// 			}

// 			nnzJ++;
// 		}
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
// 			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
// 			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
// 			//	continue;
// 			//}
// 			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
// 				continue;
// 			}
// 			//h_sparse->spJsubmatType.push_back(N);
// 			h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);
// 			h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);

// 			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

// 			//h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
// 			//h_sparse->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

// 			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 				h_sparse->Npos[nDiagN] = nnzJ;
// 				nDiagN++;
// 			}
// 			else {
// 				h_sparse->Npos[sistema->nPQ + nnzN] = (nnzJ);
// 				nnzN++;
// 			}

// 			nnzJ++;
// 		}
// 		GLOBALauxIsPQ = 0;
// 	}

// 	// M L]

// 	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
// 		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

// 			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			//h_sparse->spJsubmatType.push_back(M);

// 			h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);
// 			h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);

// 			//h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
// 			//h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

// 			if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 				h_sparse->Mpos[nDiagM] = nnzJ;
// 				nDiagM++;
// 			}
// 			else {
// 				h_sparse->Mpos[sistema->nPQ + nnzM] = (nnzJ);
// 				nnzM++;
// 			}

// 			nnzJ++;
// 		}

// 		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
// 			//int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
// 			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) /*pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim*/) { // se coluna for PQ
// 				// linha de L
// 				//h_sparse->spJsubmatType.push_back(L);
// 				//h_sparse->Lpos.push_back(nnzJ);
// 				h_sparse->cooColIndSubMatJ[nnzJ] = (sistema->csrColIndY[j]);
// 				h_sparse->cooRowIndSubMatJ[nnzJ] = (sistema->cooRowIndY[j]);

// 				int col = GLOBALauxIsPQ; //pntRowPQ - iterativo->barrasPQlim;

// 				//h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
// 				//h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

// 				if (sistema->csrColIndY[j] == sistema->cooRowIndY[j]) {
// 					h_sparse->Lpos[nDiagL] = nnzJ;
// 					nDiagL++;
// 				}
// 				else {
// 					h_sparse->Lpos[sistema->nPQ + nnzL] = (nnzJ);
// 					nnzL++;
// 				}

// 				nnzJ++;
// 			}
// 			GLOBALauxIsPQ = 0;
// 		}
// 	}
// 	h_sparse->nnzJ = nnzJ;
// 	h_sparse->spJval.reserve(nnzJ);
// 	h_sparse->spJval.resize(nnzJ);

// 	// coo to csr
// 	int /*acc = 0,*/ linha = 0;
// 	h_sparse->csrRowPtrJ.clear();
// 	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
// 	h_sparse->csrRowPtrJ[0] = 0;
// 	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;
// 	//linha++; // se é 1 based
// 	//for (int j = 1; j < h_sparse->nnzJ; j++) {
// 	//	if (h_sparse->cooRowIndJ[j] != linha /*+1...*/) {
// 	//		h_sparse->csrRowPtrJ[linha + 1] = j;
// 	//		linha++;
// 	//	}
// 	//}
// }

// coloca primeiro os eltos da diagonal principal para minimizar a divergência
void Jstencil_eficiente_(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {

	int nnzJ = 0/*, cntH = 0, cntL = 0, cntM = 0, cntN = 0*/;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); // i <- row index

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(i, sistema) - 1); // i <- row index

			if (sistema->csrColIndY[j] == i) { // i <- row index
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
			//	continue;
			//}
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); // i <- row index

			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(i, sistema) - 1); // i <- row index

			if (sistema->csrColIndY[j] == i) { // i <- row index
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
		GLOBALauxIsPQ = 0;
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(H);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); // i <- row index

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(i, sistema) - 1); // i <- row index
			if (sistema->csrColIndY[j] == i) { // i <- row index
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			//int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			//if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
			//	continue;
			//}
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(N);
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); // i <- row index

			int col = GLOBALauxIsPQ; //pntColPQ - iterativo->barrasPQlim;

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(i, sistema) - 1); // i <- row index

			if (sistema->csrColIndY[j] == i) { // i <- row index
				h_sparse->Npos[nDiagN] = nnzJ;
				nDiagN++;
			}
			else {
				h_sparse->Npos.push_back(nnzJ);
			}

			nnzJ++;
		}
		GLOBALauxIsPQ = 0;
	}

	// M L]

	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			//h_sparse->spJsubmatType.push_back(M);

			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(IDX1F(iterativo->barrasPQlim[i])); // IDX1F(iterativo->barrasPQlim[i]) <- row index

			h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->csrColIndY[j] == i) { // i <- row index
				h_sparse->Mpos[nDiagM] = nnzJ;
				nDiagM++;
			}
			else {
				h_sparse->Mpos.push_back(nnzJ);
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			//int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) /*pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim*/) { // se coluna for PQ
				// linha de L
				//h_sparse->spJsubmatType.push_back(L);
				//h_sparse->Lpos.push_back(nnzJ);
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(IDX1F(iterativo->barrasPQlim[i])); // IDX1F(iterativo->barrasPQlim[i]) <- row index

				int col = GLOBALauxIsPQ; //pntRowPQ - iterativo->barrasPQlim;

				h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->csrColIndY[j] == i) { // i <- row index
					h_sparse->Lpos[nDiagL] = nnzJ;
					nDiagL++;
				}
				else {
					h_sparse->Lpos.push_back(nnzJ);
				}

				nnzJ++;
			}
			GLOBALauxIsPQ = 0;
		}
	}
	h_sparse->nnzJ = nnzJ;
	h_sparse->spJval.reserve(nnzJ);
	h_sparse->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha /*+1...*/) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}
