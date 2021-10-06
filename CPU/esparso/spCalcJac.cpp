#include "spCalcJac.h"

using namespace std;
using namespace Eigen;

std::pair<int, int> HtoJ(int row, int col, sistemaType* sistema) {
	//std::pair<int, int> ans = 
	return { row - (row >= sistema->barraVO), col - (col >= sistema->barraVO) };
}

int HtoJrow(int row, sistemaType* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int HtoJcol(int col, sistemaType* sistema) {
	return col+1 - (col+1 >= sistema->barraVO);
}

std::pair<int, int> LtoJ(int row, int col, sistemaType* sistema) {
	//std::pair<int, int> ans = 
	return { row + sistema->nPV + sistema->nPQ, col + sistema->nPV + sistema->nPQ };
}

int LtoJrow(int row, sistemaType* sistema) {
	return row + sistema->nB - 1;
}

int LtoJcol(int col, sistemaType* sistema) {
	return col + sistema->nB - 1;
}

std::pair<int, int> MtoJ(int row, int col, sistemaType* sistema) {
	//std::pair<int, int> ans = 
	return { row + sistema->nPV + sistema->nPQ, col - (col >= sistema->barraVO) };
}

int MtoJrow(int row, sistemaType* sistema) {
	return row + sistema->nB - 1;
}

int MtoJcol(int col, sistemaType* sistema) {
	return col - (col >= sistema->barraVO);
}

std::pair<int, int> NtoJ(int row, int col, sistemaType* sistema) {
	//std::pair<int, int> ans = 
	return { row - (row >= sistema->barraVO), col + sistema->nPV + sistema->nPQ };
}

int NtoJrow(int row, sistemaType* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int NtoJcol(int col, sistemaType* sistema) {
	return col+1 + sistema->nB - 1;
}

// void Jstencil1based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {

// 	int nnzJ = 0, cntH = 0, cntL = 0, cntM = 0, cntN = 0;

// 	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	
// 	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			sistema->spJsubmatType.push_back(H);
// 			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
// 			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
// 			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema));
// 			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema));

// 			nnzJ++;
// 		}
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
// 			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
// 			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
// 				continue;
// 			}
// 			sistema->spJsubmatType.push_back(N);
// 			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j] + 1);
// 			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j] + 1);

// 			int col = pntColPQ - iterativo->barrasPQlim;

// 			sistema->cooColIndJ.push_back(NtoJcol(col, sistema));
// 			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema));

// 			nnzJ++;
// 		}
// 	}
// 	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			sistema->spJsubmatType.push_back(H);
// 			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j] + 1);
// 			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j] + 1);
// 			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema));
// 			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema));

// 			nnzJ++;
// 		}
// 		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
// 			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
// 			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
// 				continue;
// 			}
// 			sistema->spJsubmatType.push_back(N);
// 			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j] + 1);
// 			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j] + 1);

// 			int col = pntColPQ - iterativo->barrasPQlim;

// 			sistema->cooColIndJ.push_back(NtoJcol(col, sistema));
// 			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema));

// 			nnzJ++;
// 		}
// 	}

// 	// M L]

// 	for (int i = 0; i < iterativo->nPQlim; i++) { // i percorre as barras PQ (linhas do jacobiano) para montar M e L
// 		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J

// 			// linha de M (percorre linhas PQ de Y e pega todas as colunas)
// 			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
// 				continue;
// 			}
// 			sistema->spJsubmatType.push_back(M);
// 			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j] + 1);
// 			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j] + 1);
// 			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema));
// 			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema));

// 			nnzJ++;
// 		}

// 		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
// 			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
// 			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { // se coluna for PQ
// 				// linha de L
// 				sistema->spJsubmatType.push_back(L);
// 				sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j] + 1);
// 				sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j] + 1);

// 				int col = pntRowPQ - iterativo->barrasPQlim;

// 				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema));
// 				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema));
// 				nnzJ++;
// 			}
// 		}
// 	}
// 	sistema->nnzJ = nnzJ;
// 	sistema->spJval.reserve(nnzJ);
// 	sistema->spJval.resize(nnzJ);

// 	//criar cco to csr
// }

void Jstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {

	int nnzJ = 0; // , cntH = 0, cntL = 0, cntM = 0, cntN = 0;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; // pula coluna VO
			}
			sistema->spJsubmatType.push_back(H);
			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			sistema->spJsubmatType.push_back(N);
			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			sistema->spJsubmatType.push_back(H);
			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			sistema->spJsubmatType.push_back(N);
			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

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
			sistema->spJsubmatType.push_back(M);
			sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { // se coluna for PQ
				// linha de L
				sistema->spJsubmatType.push_back(L);
				sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = pntRowPQ - iterativo->barrasPQlim;

				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);
				nnzJ++;
			}
		}
	}
	sistema->nnzJ = nnzJ;
	sistema->spJval.reserve(nnzJ);
	sistema->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha /*+1...*/) {
			sistema->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

// void fillJstencil1based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
// 	//sistema->spJval.clear();
// 	//sistema->spJval.resize(sistema->nnzJ);
// 	for (int i = 0; i < sistema->nnzJ; i++) {
// 		switch (sistema->spJsubmatType[i]) {
// 			case H:
// 				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
// 					//fórmula pra km
// 					float_type aux = phif(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i]), sistema, ramo);
// 					aux += barra->theta[IDX1F(sistema->cooRowIndSubMatJ[i])] - barra->theta[IDX1F(sistema->cooColIndSubMatJ[i])];
// 					// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
// 					// Y denso
// 					// sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
// 					// Y esparso
// 					sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real() * sin(aux) - sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag() * cos(aux));
// 				}
// 				else {
// 					// Hkk = -Qk - V^2_k*Bkk
// 					// Y denso
// 					// sistema->spJval[i] = -iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y;
// 					// Y esparso
// 					sistema->spJval[i] = -iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->spY->coeff(IDX1F(sistema->cooColIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag();
// 				}
// 				break;
// 			case L:
// 				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
// 					//fórmula pra km
// 					float_type aux = phif(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i]), sistema, ramo);
// 					aux += barra->theta[IDX1F(sistema->cooRowIndSubMatJ[i])] - barra->theta[IDX1F(sistema->cooColIndSubMatJ[i])];

// 					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
// 					// Y denso
// 					// sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
// 					// Y esparso
// 					sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * (sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real() * sin(aux) - sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag() * cos(aux));
// 				}
// 				else {
// 					// Lkk= (Qk - V^2_k*Bkk)/Vk
// 					// Y denso
// 					// sistema->spJval[i] = (iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
// 					// Y esparso
// 					sistema->spJval[i] = (iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->spY->coeff(IDX1F(sistema->cooColIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag()) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
// 				}
// 				break;
// 			case M:
// 				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
// 					//fórmula pra km
// 					float_type aux = phif(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i]), sistema, ramo);
// 					aux += barra->theta[IDX1F(sistema->cooRowIndSubMatJ[i])] - barra->theta[IDX1F(sistema->cooColIndSubMatJ[i])];

// 					// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
// 					// Y denso
// 					// sistema->spJval[i] = -barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
// 					// Y esparso
// 					sistema->spJval[i] = -barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real() * cos(aux) + sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag() * sin(aux));
// 				}
// 				else {
// 					// Mkk= (Pk - V^2_k*Gkk)
// 					// Y denso
// 					// sistema->spJval[i] = iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x;
// 					// Y esparso
// 					sistema->spJval[i] = iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->spY->coeff(IDX1F(sistema->cooColIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real();
// 				}
// 				break;
// 			case N:
// 				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
// 					//fórmula pra km=iidx
// 					float_type aux = phif(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i]), sistema, ramo);
// 					aux += barra->theta[IDX1F(sistema->cooRowIndSubMatJ[i])] - barra->theta[IDX1F(sistema->cooColIndSubMatJ[i])];

// 					// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
// 					// Y denso
// 					// sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
// 					// Y esparso
// 					sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * (sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real() * cos(aux) + sistema->spY->coeff(IDX1F(sistema->cooRowIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).imag() * sin(aux));
// 				}
// 				else {
// 					// Nkk= (Qk - V^2_k*Bkk)/Vkk
// 					// Y denso
// 					// sistema->spJval[i] = (iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] + barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
// 					// Y esparso
// 					sistema->spJval[i] = (iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] + barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->spY->coeff(IDX1F(sistema->cooColIndSubMatJ[i]), IDX1F(sistema->cooColIndSubMatJ[i])).real()) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
// 				}
// 				break;
// 			default:
// 				printf("ERRO [fillJstencil] tipo inválido de H");
// 				break;
// 		}
// 	}
// }

void fillJstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	//sistema->spJval.clear();
	//sistema->spJval.resize(sistema->nnzJ);
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nnzJ; i++) {
			switch (sistema->spJsubmatType[i]) {
			case H:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					//fórmula pra km
					//float_type aux = phif(sistema->cooRowIndSubMatJ[i] + 1, sistema->cooColIndSubMatJ[i] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];
					// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
					// Y denso
					// sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
					// Y esparso
					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					// Hkk = -Qk - V^2_k*Bkk
					// Y denso
					// sistema->spJval[i] = -iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y;
					// Y esparso
					sistema->spJval[i] = -iterativo->Qcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag();
				}
				break;
			case L:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					//fórmula pra km
					//float_type aux = phif(sistema->cooRowIndSubMatJ[i] + 1, sistema->cooColIndSubMatJ[i] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
					// Y denso
					// sistema->spJval[i] = barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * sin(aux) - sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * cos(aux));
					// Y esparso
					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					// Lkk= (Qk - V^2_k*Bkk)/Vk
					// Y denso
					// sistema->spJval[i] = (iterativo->Qcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
					// Y esparso
					sistema->spJval[i] = (iterativo->Qcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag()) / barra->V[sistema->cooColIndSubMatJ[i]];
				}
				break;
			case M:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					//fórmula pra km
					// float_type aux = phif(sistema->cooRowIndSubMatJ[i] + 1, sistema->cooColIndSubMatJ[i] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
					// Y denso
					// sistema->spJval[i] = -barra->V[IDX1F(sistema->cooRowIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
					// Y esparso
					sistema->spJval[i] = -barra->V[sistema->cooRowIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					// Mkk= (Pk - V^2_k*Gkk)
					// Y denso
					// sistema->spJval[i] = iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] - barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x;
					// Y esparso
					sistema->spJval[i] = iterativo->Pcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real();
				}
				break;
			case N:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					//fórmula pra km=iidx
					// float_type aux = phif(sistema->cooRowIndSubMatJ[i] + 1, sistema->cooColIndSubMatJ[i] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
					// Y denso
					//sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * (sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x * cos(aux) + sistema->Y[IDX2F(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].y * sin(aux));
					// Y esparso
					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					// Nkk= (Qk - V^2_k*Bkk)/Vkk
					// Y denso
					// sistema->spJval[i] = (iterativo->Pcalc[IDX1F(sistema->cooColIndSubMatJ[i])] + barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * barra->V[IDX1F(sistema->cooColIndSubMatJ[i])] * sistema->Y[IDX2F(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i], sistema->nB)].x) / barra->V[IDX1F(sistema->cooColIndSubMatJ[i])];
					// Y esparso
					sistema->spJval[i] = (iterativo->Pcalc[sistema->cooColIndSubMatJ[i]] + barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real()) / barra->V[sistema->cooColIndSubMatJ[i]];
				}
				break;
			default:
				printf("ERRO [fillJstencil] tipo inválido de H");
				break;
			}
		}
}

void spCalcJ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {

	sistema->spJval.clear();
	sistema->cooColIndJ.clear();
	sistema->cooRowIndJ.clear();

	int nnzJ = 0; // , cntH = 0, cntL = 0, cntM = 0, cntN = 0;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; // pula coluna VO
			}
			//sistema->spJsubmatType.push_back(H);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				// Y esparso
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			//sistema->spJsubmatType.push_back(N);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km=iidx
				//float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Nkk= (Qk - V^2_k*Bkk)/Vkk
				// Y esparso
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
			}

			nnzJ++;
		}
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; // pula coluna VO
			}
			//sistema->spJsubmatType.push_back(H);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				// Y esparso
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de N
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procure sistema->csrColIndY[cntN] em iterativo->barrasPQlim
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { // apenas prossegue se é coluna PQ de Y
				continue;
			}
			//sistema->spJsubmatType.push_back(N);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km=iidx
				//float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Nkk= (Qk - V^2_k*Bkk)/Vkk
				// Y esparso
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
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
			// sistema->spJsubmatType.push_back(M);
			// sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			// sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
				// Y esparso
				sistema->spJval.push_back(-barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Mkk= (Pk - V^2_k*Gkk)
				// Y esparso
				sistema->spJval.push_back(iterativo->Pcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real());
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { // se coluna for PQ
				// linha de L
				//sistema->spJsubmatType.push_back(L);
				//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = pntRowPQ - iterativo->barrasPQlim;

				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
					//fórmula pra km
					// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
					// Y esparso
					sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
				}
				else {
					// Lkk= (Qk - V^2_k*Bkk)/Vk
					// Y esparso
					sistema->spJval.push_back((iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag()) / barra->V[sistema->csrColIndY[j]]);
				}
				nnzJ++;
			}
		}
	}
	sistema->nnzJ = nnzJ;
	// sistema->spJval.reserve(nnzJ);
	// sistema->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha /*+1...*/) {
			sistema->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

int GLOBALauxIsPQ = 0;
bool isPQ(int barra, sistemaType* sistema, iterativoType* iterativo) {
	while (iterativo->barrasPQlim[GLOBALauxIsPQ] < barra) {
		GLOBALauxIsPQ++;
		if (GLOBALauxIsPQ == sistema->nPQ) {
			return false;
		}
	}
	return iterativo->barrasPQlim[GLOBALauxIsPQ] == barra;
}

void spCalcJ_eficiente(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {

	sistema->spJval.clear();
	sistema->cooColIndJ.clear();
	sistema->cooRowIndJ.clear();

	int nnzJ = 0; //, cntH = 0, cntL = 0, cntM = 0, cntN = 0;

	// cabe otimização da parte de N com buscas no vetor coo de Y pelos elementos PQ de forma direta ?
	// segundo profiler do visual studio: parte mais onerosa são os push_backs. Contribuição insignificante para o resultado final (até ieee118x8).
	// loops for não são paralelizáveis (dependem do valor de nnzJ)

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { // i percorre as barras diretamente e pula a barra swing (linhas do jacobiano) para montar H e N
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; // pula coluna VO
			}
			//sistema->spJsubmatType.push_back(H);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				//float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				// Y esparso
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
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
			//sistema->spJsubmatType.push_back(N);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; // pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km=iidx
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Nkk= (Qk - V^2_k*Bkk)/Vkk
				// Y esparso
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
			}

			nnzJ++;
		}
		GLOBALauxIsPQ = 0;
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { // i percorre as barras diretamente e pula a linha d Y referente à barra swing
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { // linha de H
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; // pula coluna VO
			}
			//sistema->spJsubmatType.push_back(H);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				// Hkm = Vk*Vm * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				// Hkk = -Qk - V^2_k*Bkk
				// Y esparso
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
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
			//sistema->spJsubmatType.push_back(N);
			//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

			int col = GLOBALauxIsPQ; // pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km=iidx
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// Vk * (Gkm*sen(theta) - Bkm*cos(theta))
				// Y esparso
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Nkk= (Qk - V^2_k*Bkk)/Vkk
				// Y esparso
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
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
			// sistema->spJsubmatType.push_back(M);
			// sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
			// sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);
			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				//fórmula pra km
				// float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				// -Vk*Vm * (Gkm*cos(theta) + Bkm*sen(theta))
				// Y esparso
				sistema->spJval.push_back(-barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				// Mkk= (Pk - V^2_k*Gkk)
				// Y esparso
				sistema->spJval.push_back(iterativo->Pcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real());
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { // j percorre i-ésima linha de Y para preencher uma linha de J
			//int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); // procura sistema->csrColIndY[j] em iterativo->barrasPQlim i.e., elto da coluna também é PQ?
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) /*pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim*/) { // se coluna for PQ
				// linha de L
				//sistema->spJsubmatType.push_back(L);
				//sistema->cooColIndSubMatJ.push_back(sistema->csrColIndY[j]);
				//sistema->cooRowIndSubMatJ.push_back(sistema->cooRowIndY[j]);

				int col = GLOBALauxIsPQ; //pntRowPQ - iterativo->barrasPQlim;

				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
					//fórmula pra km
					//float_type aux = phif(sistema->cooRowIndY[j] + 1, sistema->csrColIndY[j] + 1, sistema, ramo);
					float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

					// Vk * (Gkm*sen(theta) - Bkmcos(theta))
					// Y esparso
					sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
				}
				else {
					// Lkk= (Qk - V^2_k*Bkk)/Vk
					// Y esparso
					sistema->spJval.push_back((iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag()) / barra->V[sistema->csrColIndY[j]]);
				}
				nnzJ++;
			}
			GLOBALauxIsPQ = 0;
		}
	}
	sistema->nnzJ = nnzJ;
	// sistema->spJval.reserve(nnzJ);
	// sistema->spJval.resize(nnzJ);

	// coo to csr
	int /*acc = 0,*/ linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;
	//linha++; // se é 1 based
	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha /*+1...*/) {
			sistema->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}