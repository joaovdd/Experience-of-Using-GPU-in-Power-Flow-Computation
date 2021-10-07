#pragma once

#include <algorithm>
#include <vector>

#include <stdlib.h>
#include <stdio.h>

#include "Eigen/Sparse"
using namespace std;
using namespace Eigen;

std::pair<int, int> HtoJ(int row, int col, sistema* sistema) {
	return { row - (row >= sistema->barraVO), col - (col >= sistema->barraVO) };
}

int HtoJrow(int row, sistema* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int HtoJcol(int col, sistema* sistema) {
	return col+1 - (col+1 >= sistema->barraVO);
}

std::pair<int, int> LtoJ(int row, int col, sistema* sistema) {
	return { row + sistema->nPV + sistema->nPQ, col + sistema->nPV + sistema->nPQ };
}

int LtoJrow(int row, sistema* sistema) {
	return row + sistema->nB - 1;
}

int LtoJcol(int col, sistema* sistema) {
	return col + sistema->nB - 1;
}

std::pair<int, int> MtoJ(int row, int col, sistema* sistema) {
	return { row + sistema->nPV + sistema->nPQ, col - (col >= sistema->barraVO) };
}

int MtoJrow(int row, sistema* sistema) {
	return row + sistema->nB - 1;
}

int MtoJcol(int col, sistema* sistema) {
	return col - (col >= sistema->barraVO);
}

std::pair<int, int> NtoJ(int row, int col, sistema* sistema) {
	return { row - (row >= sistema->barraVO), col + sistema->nPV + sistema->nPQ };
}

int NtoJrow(int row, sistema* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int NtoJcol(int col, sistema* sistema) {
	return col+1 + sistema->nB - 1;
}

void Jstencil0based(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int nnzJ = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

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

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { 
				
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

	int  linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;

	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha ) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

void Jstencil(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int nnzJ = 0;

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

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

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { 
				
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

	int  linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;

	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha ) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

void fillJstencil0based(sistema* sistema, h_sparse* h_sparse, barra* barra, ramo* ramo, iterativo* iterativo) {
	#pragma omp parallel for if (openmp)
		for (int i = 0; i < h_sparse->nnzJ; i++) {
			switch (h_sparse->spJsubmatType[i]) {
			case H:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];
					
					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					
					h_sparse->spJval[i] = -iterativo->Qcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag();
				}
				break;
			case L:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					
					h_sparse->spJval[i] = (iterativo->Qcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag()) / barra->V[h_sparse->cooColIndSubMatJ[i]];
				}
				break;
			case M:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					h_sparse->spJval[i] = -barra->V[h_sparse->cooRowIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					
					h_sparse->spJval[i] = iterativo->Pcalc[h_sparse->cooColIndSubMatJ[i]] - barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real();
				}
				break;
			case N:
				if (h_sparse->cooRowIndSubMatJ[i] != h_sparse->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[h_sparse->cooRowIndSubMatJ[i]] - barra->theta[h_sparse->cooColIndSubMatJ[i]];

					h_sparse->spJval[i] = barra->V[h_sparse->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(h_sparse->cooRowIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					
					h_sparse->spJval[i] = (iterativo->Pcalc[h_sparse->cooColIndSubMatJ[i]] + barra->V[h_sparse->cooColIndSubMatJ[i]] * barra->V[h_sparse->cooColIndSubMatJ[i]] * sistema->spY->coeff(h_sparse->cooColIndSubMatJ[i], h_sparse->cooColIndSubMatJ[i]).real()) / barra->V[h_sparse->cooColIndSubMatJ[i]];
				}
				break;
			default:
				printf("ERRO [fillJstencil] tipo inválido de H");
				break;
			}
		}
}

void reordJac(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int szHdp = 0, szLdp = 0, szMdp = 0, szNdp = 0, aux;

	for (int i = 0; i < h_sparse->Hpos.size(); i++) {
		
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Hpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Hpos[i]]) {
			aux = h_sparse->Hpos[i];
			
			for (int j = i - 1; j >= szHdp; j--) {
				h_sparse->Hpos[j + 1] = h_sparse->Hpos[j];
			}
			h_sparse->Hpos[szHdp] = aux;
			szHdp++;
		}
	}

	for (int i = 0; i < h_sparse->Lpos.size(); i++) {
		
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Lpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Lpos[i]]) {
			aux = h_sparse->Lpos[i];
			
			for (int j = i - 1; j >= szLdp; j--) {
				h_sparse->Lpos[j + 1] = h_sparse->Lpos[j];
			}
			h_sparse->Lpos[szLdp] = aux;
			szLdp++;
		}
	}

	for (int i = 0; i < h_sparse->Mpos.size(); i++) {
		
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Mpos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Mpos[i]]) {
			aux = h_sparse->Mpos[i];
			
			for (int j = i - 1; j >= szMdp; j--) {
				h_sparse->Mpos[j + 1] = h_sparse->Mpos[j];
			}
			h_sparse->Mpos[szMdp] = aux;
			szMdp++;
		}
	}

	for (int i = 0; i < h_sparse->Npos.size(); i++) {
		
		if (h_sparse->cooRowIndSubMatJ[h_sparse->Npos[i]] == h_sparse->cooColIndSubMatJ[h_sparse->Npos[i]]) {
			aux = h_sparse->Npos[i];
			
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

void Jstencil_eficiente(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int nnzJ = 0;

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; 

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
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; 

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

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
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

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) ) { 
				
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = GLOBALauxIsPQ; 

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

	int  linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;

	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha ) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

void Jstencil_eficiente_(sistema* sistema, h_sparse* h_sparse, barra* barra, iterativo* iterativo) {
	int nnzJ = 0;

	h_sparse->Hpos.resize(sistema->nB - 1);
	int nDiagH = 0;
	h_sparse->Lpos.resize(sistema->nPQ);
	int nDiagL = 0;
	h_sparse->Mpos.resize(sistema->nPQ);
	int nDiagM = 0;
	h_sparse->Npos.resize(sistema->nPQ);
	int nDiagN = 0;

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); 

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(i, sistema) - 1); 

			if (sistema->csrColIndY[j] == i) { 
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); 

			int col = GLOBALauxIsPQ; 

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(i, sistema) - 1); 

			if (sistema->csrColIndY[j] == i) { 
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
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); 

			h_sparse->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			h_sparse->cooRowIndJ.push_back(HtoJrow(i, sistema) - 1); 
			if (sistema->csrColIndY[j] == i) { 
				h_sparse->Hpos[nDiagH] = nnzJ;
				nDiagH++;
			}
			else {
				h_sparse->Hpos.push_back(nnzJ);
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(i); 

			int col = GLOBALauxIsPQ; 

			h_sparse->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(NtoJrow(i, sistema) - 1); 

			if (sistema->csrColIndY[j] == i) { 
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

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
			h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			h_sparse->cooRowIndSubMatJ.push_back(IDX1F(iterativo->barrasPQlim[i])); 

			h_sparse->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			h_sparse->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->csrColIndY[j] == i) { 
				h_sparse->Mpos[nDiagM] = nnzJ;
				nDiagM++;
			}
			else {
				h_sparse->Mpos.push_back(nnzJ);
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) ) { 
				
				h_sparse->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				h_sparse->cooRowIndSubMatJ.push_back(IDX1F(iterativo->barrasPQlim[i])); 

				int col = GLOBALauxIsPQ; 

				h_sparse->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				h_sparse->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->csrColIndY[j] == i) { 
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

	int  linha = 0;
	h_sparse->csrRowPtrJ.clear();
	h_sparse->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	h_sparse->csrRowPtrJ[0] = 0;
	h_sparse->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = h_sparse->nnzJ;

	for (int j = 1; j < h_sparse->nnzJ; j++) {
		if (h_sparse->cooRowIndJ[j] != linha ) {
			h_sparse->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}
