#include "spCalcJac.h"

using namespace std;
using namespace Eigen;

std::pair<int, int> HtoJ(int row, int col, sistemaType* sistema) {
	return { row - (row >= sistema->barraVO), col - (col >= sistema->barraVO) };
}

int HtoJrow(int row, sistemaType* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int HtoJcol(int col, sistemaType* sistema) {
	return col+1 - (col+1 >= sistema->barraVO);
}

std::pair<int, int> LtoJ(int row, int col, sistemaType* sistema) {
	return { row + sistema->nPV + sistema->nPQ, col + sistema->nPV + sistema->nPQ };
}

int LtoJrow(int row, sistemaType* sistema) {
	return row + sistema->nB - 1;
}

int LtoJcol(int col, sistemaType* sistema) {
	return col + sistema->nB - 1;
}

std::pair<int, int> MtoJ(int row, int col, sistemaType* sistema) {
	return { row + sistema->nPV + sistema->nPQ, col - (col >= sistema->barraVO) };
}

int MtoJrow(int row, sistemaType* sistema) {
	return row + sistema->nB - 1;
}

int MtoJcol(int col, sistemaType* sistema) {
	return col - (col >= sistema->barraVO);
}

std::pair<int, int> NtoJ(int row, int col, sistemaType* sistema) {
	return { row - (row >= sistema->barraVO), col + sistema->nPV + sistema->nPQ };
}

int NtoJrow(int row, sistemaType* sistema) {
	return row+1 - (row+1 >= sistema->barraVO);
}

int NtoJcol(int col, sistemaType* sistema) {
	return col+1 + sistema->nB - 1;
}

void Jstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	int nnzJ = 0; 

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
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
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
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

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

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

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { 
				
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

	int  linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;

	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha ) {
			sistema->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}

void fillJstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	#pragma omp parallel for if (global::openmp)
		for (int i = 0; i < sistema->nnzJ; i++) {
			switch (sistema->spJsubmatType[i]) {
			case H:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];
					
					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					
					sistema->spJval[i] = -iterativo->Qcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag();
				}
				break;
			case L:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * cos(aux));
				}
				else {
					
					sistema->spJval[i] = (iterativo->Qcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag()) / barra->V[sistema->cooColIndSubMatJ[i]];
				}
				break;
			case M:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					sistema->spJval[i] = -barra->V[sistema->cooRowIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					
					sistema->spJval[i] = iterativo->Pcalc[sistema->cooColIndSubMatJ[i]] - barra->V[sistema->cooColIndSubMatJ[i]] * barra->V[sistema->cooColIndSubMatJ[i]] * sistema->spY->coeff(sistema->cooColIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real();
				}
				break;
			case N:
				if (sistema->cooRowIndSubMatJ[i] != sistema->cooColIndSubMatJ[i]) {
					
					float_type aux = barra->theta[sistema->cooRowIndSubMatJ[i]] - barra->theta[sistema->cooColIndSubMatJ[i]];

					sistema->spJval[i] = barra->V[sistema->cooRowIndSubMatJ[i]] * (sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndSubMatJ[i], sistema->cooColIndSubMatJ[i]).imag() * sin(aux));
				}
				else {
					
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

	int nnzJ = 0; 

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; 
			}
			
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
				continue;
			}
			
			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
			}

			nnzJ++;
		}
	}
	for (int i = IDX1F(sistema->barraVO) + 1; i < sistema->nB; i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; 
			}
			
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			int* pntColPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntColPQ == iterativo->barrasPQlim + iterativo->nPQlim) { 
				continue;
			}
			
			int col = pntColPQ - iterativo->barrasPQlim;

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
			}

			nnzJ++;
		}
	}

	for (int i = 0; i < iterativo->nPQlim; i++) { 
		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 

			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue;
			}
			
			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(-barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back(iterativo->Pcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real());
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			int* pntRowPQ = std::find(iterativo->barrasPQlim, iterativo->barrasPQlim + iterativo->nPQlim, sistema->csrColIndY[j] + 1); 
			if (pntRowPQ != iterativo->barrasPQlim + iterativo->nPQlim) { 
				
				int col = pntRowPQ - iterativo->barrasPQlim;

				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
					
					float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

					sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
				}
				else {
					
					sistema->spJval.push_back((iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag()) / barra->V[sistema->csrColIndY[j]]);
				}
				nnzJ++;
			}
		}
	}
	sistema->nnzJ = nnzJ;

	int  linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;

	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha ) {
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

	int nnzJ = 0; 

	for (int i = 0; i < IDX1F(sistema->barraVO); i++) { 
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			if (sistema->csrColIndY[j] == IDX1F(sistema->barraVO)) {
				continue; 
			}
			
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			int col = GLOBALauxIsPQ; 

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
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
			
			sistema->cooColIndJ.push_back(HtoJcol(sistema->csrColIndY[j], sistema) - 1);
			sistema->cooRowIndJ.push_back(HtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];
				
				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
			}
			else {
				
				sistema->spJval.push_back(-iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag());
			}

			nnzJ++;
		}
		for (int j = sistema->csrRowPtrY[i]; j < sistema->csrRowPtrY[i + 1]; j++) { 
			
			if (!isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo)) {
				continue;
			}
			
			int col = GLOBALauxIsPQ; 

			sistema->cooColIndJ.push_back(NtoJcol(col, sistema) - 1);
			sistema->cooRowIndJ.push_back(NtoJrow(sistema->cooRowIndY[j], sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back((iterativo->Pcalc[sistema->csrColIndY[j]] + barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real()) / barra->V[sistema->csrColIndY[j]]);
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
			
			sistema->cooColIndJ.push_back(MtoJcol(sistema->csrColIndY[j] + 1, sistema) - 1);
			sistema->cooRowIndJ.push_back(MtoJrow(i + 1, sistema) - 1);

			if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
				
				float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

				sistema->spJval.push_back(-barra->V[sistema->cooRowIndY[j]] * barra->V[sistema->csrColIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * cos(aux) + sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * sin(aux)));
			}
			else {
				
				sistema->spJval.push_back(iterativo->Pcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).real());
			}

			nnzJ++;
		}

		for (int j = sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i])]; j < sistema->csrRowPtrY[IDX1F(iterativo->barrasPQlim[i]) + 1]; j++) { 
			
			if (isPQ(sistema->csrColIndY[j] + 1, sistema, iterativo) ) { 
				
				int col = GLOBALauxIsPQ; 

				sistema->cooColIndJ.push_back(LtoJcol(col + 1, sistema) - 1);
				sistema->cooRowIndJ.push_back(LtoJrow(i + 1, sistema) - 1);

				if (sistema->cooRowIndY[j] != sistema->csrColIndY[j]) {
					
					float_type aux = barra->theta[sistema->cooRowIndY[j]] - barra->theta[sistema->csrColIndY[j]];

					sistema->spJval.push_back(barra->V[sistema->cooRowIndY[j]] * (sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).real() * sin(aux) - sistema->spY->coeff(sistema->cooRowIndY[j], sistema->csrColIndY[j]).imag() * cos(aux)));
				}
				else {
					
					sistema->spJval.push_back((iterativo->Qcalc[sistema->csrColIndY[j]] - barra->V[sistema->csrColIndY[j]] * barra->V[sistema->csrColIndY[j]] * sistema->spY->coeff(sistema->csrColIndY[j], sistema->csrColIndY[j]).imag()) / barra->V[sistema->csrColIndY[j]]);
				}
				nnzJ++;
			}
			GLOBALauxIsPQ = 0;
		}
	}
	sistema->nnzJ = nnzJ;

	int  linha = 0;
	sistema->csrRowPtrJ.clear();
	sistema->csrRowPtrJ.resize(sistema->nB - 1 + iterativo->nPQlim + 1);
	sistema->csrRowPtrJ[0] = 0;
	sistema->csrRowPtrJ[sistema->nB - 1 + iterativo->nPQlim] = sistema->nnzJ;

	for (int j = 1; j < sistema->nnzJ; j++) {
		if (sistema->cooRowIndJ[j] != linha ) {
			sistema->csrRowPtrJ[linha + 1] = j;
			linha++;
		}
	}
}