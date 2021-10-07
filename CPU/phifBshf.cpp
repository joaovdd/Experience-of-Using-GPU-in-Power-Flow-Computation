#include "phifBshf.h"
#include "idx.h"

float_type dnPhif(int a, int b, sistemaType* sistema, ramoType* ramo) {
	for (int i = 0; i < sistema->nL; i++) { 
		if (ramo->de[i] == a) { 
			if (ramo->para[i] == b) { 
				return ((ramo->phi[i]));
			}
		}
		else {
			if (ramo->para[i] == a) { 
				if (ramo->de[i] == b) { 
					return (-(ramo->phi[i]));
				}
			}
		}
	}
	return(0.);
}

int coeffPos(int lin, int col, int* csrRowPtrY, int* csrColIndY, int nnzY) {
	for (int i = csrRowPtrY[lin]; i < csrRowPtrY[lin + 1]; i++) {
		if (csrColIndY[i] == col) {
			return i;
		}
	}
	return -1;
}

float_type phif(int row, int col, sistemaType* sistema, ramoType* ramo) {
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
	for (int i = 0; i < sistema->nL; i++) { 
		if (ramo->de[i] == a) { 
			if (ramo->para[i] == b) { 
				return (ramo->bsh[i]);
			}
		}
		else {
			if (ramo->para[i] == a) { 
				if (ramo->de[i] == b) { 
					return (ramo->bsh[i]);
				}
			}
		}
	}
	return(0.);
}
