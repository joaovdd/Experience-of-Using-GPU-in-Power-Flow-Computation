#include "phifBshf.h"
#include "idx.h"

// criar matriz esparsa para phif e bshf via triplets (ramo->de[i], ramo->para[i], ramo->phi[i])
float_type dnPhif(int a, int b, sistemaType* sistema, ramoType* ramo) {
	for (int i = 0; i < sistema->nL; i++) { // i percorre ramos
		if (ramo->de[i] == a) { // se k est� ligada � para[i] =>
			if (ramo->para[i] == b) { // se k est� ligada � para[i] =>
				return ((ramo->phi[i]));
			}
		}
		else {
			if (ramo->para[i] == a) { // se k est� ligada � para[i] =>
				if (ramo->de[i] == b) { // se k est� ligada � para[i] =>
					return (-(ramo->phi[i]));
				}
			}
		}
	}
	return(0.);
}

// retorna a posição do elemento de índices (lin, col) no vetor de valores da matriz
// armazenada da forma csr (overload para indexadores int)
int coeffPos(int lin, int col, int* csrRowPtrY, int* csrColIndY, int nnzY) {
	for (int i = csrRowPtrY[lin]; i < csrRowPtrY[lin + 1]; i++) {
		if (csrColIndY[i] == col) {
			return i;
		}
	}
	return -1;
}

// retorna valor da defasagem angular entre as barras a e b
// pode se tornar mais eficiente (será?) ao se utilizar a função coeff (O(log(nnz_j)).
float_type phif(int row, int col, sistemaType* sistema, ramoType* ramo) {
	//return ramo->eigen_phi.coeff(IDX1F(row), IDX1F(col));
	auto aux = coeffPos(IDX1F(row), IDX1F(col), ramo->eigen_phi.outerIndexPtr(), ramo->eigen_phi.innerIndexPtr(), ramo->eigen_phi.nonZeros());
	if (aux == -1) {
		aux = coeffPos(IDX1F(col), IDX1F(row), ramo->eigen_phi.outerIndexPtr(), ramo->eigen_phi.innerIndexPtr(), ramo->eigen_phi.nonZeros());
		if (aux == -1) {
			return 0.;
		}
		else {
			return -ramo->eigen_phi.valuePtr()[aux];
		}
	}
	else {
		return ramo->eigen_phi.valuePtr()[aux];
	}
}

float_type bshf(int a, int b, sistemaType* sistema, ramoType* ramo) {
	for (int i = 0; i < sistema->nL; i++) { // i percorre ramos
		if (ramo->de[i] == a) { // se k est� ligada � para[i] =>
			if (ramo->para[i] == b) { // se k est� ligada � para[i] =>
				return (ramo->bsh[i]);
			}
		}
		else {
			if (ramo->para[i] == a) { // se k est� ligada � para[i] =>
				if (ramo->de[i] == b) { // se k est� ligada � para[i] =>
					return (ramo->bsh[i]);
				}
			}
		}
	}
	return(0.);
}
