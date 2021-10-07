#include "dim.h"

void showVec(bool* show, int dim, int prec) {
	std::cout << '[' << show[0] << std::endl;

	for (int i = 1; i < dim - 1; i++) {
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

void showVec(float_type* show, int dim, int prec) {
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (int i = 1; i < dim - 1; i++) {
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

void showVec(int* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

void showVecR(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x << std::endl;

	for (int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].x << std::endl;
	}

	std::cout << ' ' << show[dim - 1].x << ']' << std::endl;
}

void showVecI(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].y << std::endl;

	for (int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].y << std::endl;
	}

	std::cout << ' ' << show[dim - 1].y << ']' << std::endl;
}

void showVec(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x << " + j*" << show[0].y << std::endl;

	for (int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].x << " + j*" << show[i].y << std::endl;
	}

	std::cout << ' ' << show[dim - 1].x << " + j*" << show[dim - 1].y << ']' << std::endl;
}

void showMatI(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].y;

	for (int j = 1; j < dim; j++){
    	for (int i = 1; i < dim; i++){
	    	
			printf(" %8.5f", show[j + i * dim].y);
	    }
		std::cout << std::endl;
	}
	std::cout << ']' << std::endl;
}

void showMatR(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x;

	for (int j = 1; j < dim; j++){
    	for (int i = 1; i < dim; i++){
	    	
			printf(" %10.5f", show[j + i * dim].x);
	    }
		std::cout << std::endl;
	}
	std::cout<< ']' << std::endl;
}

void showMat(const complex_type* show, const int dim){
	printf("\n[");

	for (int j = 0; j < dim; j++){
    	for (int i = 0; i < dim; i++){
    		printf("(%5.1f + j*%5.1f) ", show[j+i*dim].x ,show[j+i*dim].y);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

void showMat(const float_type* show, const int dim) {
	printf("\n[");

	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < dim; i++) {
			printf("%5.1f ", show[j + i * dim]);
		}
		printf("\n ");
	}
	printf("]\n\n");
}

void showMatRI(const complex_type* show, const int dim) {
	printf("\nRe=\n[");

	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < dim; i++) {
			printf("%8.4f ", show[j + i * dim].x);
		}
		printf("\n ");
	}
	printf("]\n\n");

	printf("\nIm=\n[");

	for (int j = 0; j < dim; j++) {
		for (int i = 0; i < dim; i++) {
			printf("%8.4f ", show[j + i * dim].y);
		}
		printf("\n ");
	}
	printf("]\n\n");
}

void printAll(sistemaType& sistema, barraType &barra, ramoType &ramo){
	std::cout<< "nB  = " << sistema.nB << std::endl;
	std::cout<< "nPQ = " << sistema.nPQ << std::endl;
	std::cout<< "nPV = " << sistema.nPV << std::endl;
	std::cout<< "nl  = " << sistema.nL << std::endl;

	showVec(barra.Qload, sistema.nB, 4);
	std::cout<< "Pg:" << std::endl;
	showVec(barra.Pg, sistema.nB, 4);
	std::cout<< "Qg:" << std::endl;
	showVec(barra.Qg, sistema.nB, 4);
	std::cout<< "Pliq:" << std::endl;
	showVec(barra.Pliq, sistema.nB, 4);
	std::cout<< "Qliq:" << std::endl;
	showVec(barra.Qliq, sistema.nB, 4);
	std::cout<< "VBase:" << std::endl;
	showVec(barra.Vbase, sistema.nB, 4);
	std::cout<< "V:" << std::endl;
	showVec(barra.V, sistema.nB, 4);
	std::cout<< "theta:" << std::endl;
	showVec(barra.theta, sistema.nB, 4);
	std::cout<< "bsh(barra):" << std::endl;
	showVec(barra.bsh, sistema.nB, 4);
	std::cout<< "gsh(barra):" << std::endl;
	showVec(barra.gsh, sistema.nB, 4);

	std::cout<< "de:" << std::endl;
	showVec(ramo.de, sistema.nL, 1);
	std::cout<< "para:" << std::endl;
	showVec(ramo.para, sistema.nL, 1);
	std::cout<< "bsh(ramo):" << std::endl;
	showVec(ramo.bsh, sistema.nL, 4);
	std::cout<< "phi(trafo):" << std::endl;
	showVec(ramo.phi, sistema.nL, 4);

	std::cout<< "z:" << std::endl;	
	showVec(ramo.z, sistema.nL, 4);

	std::cout<< "barrasPV :" << std::endl;
	showVec(sistema.barrasPV, sistema.nPV, 1);
	std::cout<< "barrasPQ :" << std::endl;
	showVec(sistema.barrasPQ, sistema.nPQ, 1);
	std::cout<< "barraVO : " << sistema.barraVO << std::endl;

	std::cout<< "Ybus:" << std::endl;
	if (global::metodo != esparso) {
		showMat(sistema.Y, sistema.nPQ+sistema.nPV+1);

		showMatRI(sistema.Y, sistema.nPQ + sistema.nPV + 1);
	}
	else{
		std::cout<< *sistema.spY << std::endl;
	}
}
