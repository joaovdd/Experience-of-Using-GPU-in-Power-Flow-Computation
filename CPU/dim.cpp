#include "dim.h"

//      o que mostrar, dimens�o, no de casas decimais
void showVec(bool* show, unsigned int dim, unsigned int prec) {
	//std::cout.unsetf(std::ios::floatfield);
	//std::cout.setf(std::ios::fixed, std::ios::floatfield);
	//std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++) {
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(float_type* show, unsigned int dim, unsigned int prec) {
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++) {
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(unsigned int* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(int* show, int dim, unsigned int prec) {
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0] << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++) {
		std::cout << ' ' << show[i] << std::endl;
	}

	std::cout << ' ' << show[dim - 1] << ']' << std::endl;
}


//      o que mostrar, dimens�o, no de casas decimais
void showVecR(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].x << std::endl;
	}

	std::cout << ' ' << show[dim - 1].x << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVecI(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].y << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].y << std::endl;
	}

	std::cout << ' ' << show[dim - 1].y << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x << " + j*" << show[0].y << std::endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		std::cout << ' ' << show[i].x << " + j*" << show[i].y << std::endl;
	}

	std::cout << ' ' << show[dim - 1].x << " + j*" << show[dim - 1].y << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMatI(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].y;

	for (unsigned int j = 1; j < dim; j++){
    	for (unsigned int i = 1; i < dim; i++){
	    	//cout << ' ' << show[j+i*dim].y;
			printf(" %8.5f", show[j + i * dim].y);
	    }
		std::cout << std::endl;
	}
	std::cout << ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMatR(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	std::cout << '[' << show[0].x;

	for (unsigned int j = 1; j < dim; j++){
    	for (unsigned int i = 1; i < dim; i++){
	    	//cout << ' ' << show[j+i*dim].x;
			printf(" %10.5f", show[j + i * dim].x);
	    }
		std::cout << std::endl;
	}
	std::cout<< ']' << std::endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMat(const complex_type* show, const unsigned int dim){
	printf("\n[");

	for (unsigned int j = 0; j < dim; j++){
    	for (unsigned int i = 0; i < dim; i++){
    		printf("(%5.1f + j*%5.1f) ", show[j+i*dim].x ,show[j+i*dim].y);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

//      o que mostrar, dimens�o, no de casas decimais
void showMat(const float_type* show, const unsigned int dim) {
	printf("\n[");

	for (unsigned int j = 0; j < dim; j++) {
		for (unsigned int i = 0; i < dim; i++) {
			printf("%5.1f ", show[j + i * dim]);
		}
		printf("\n ");
	}
	printf("]\n\n");
}

//// o que mostrar, dimens�o, no de casas decimais
//void showMat(const float* show, const unsigned int dim) {
//	printf("\n[");
//
//	for (unsigned int j = 0; j < dim; j++) {
//		for (unsigned int i = 0; i < dim; i++) {
//			printf("%5.1f ", show[j + i * dim]);
//		}
//		printf("\n ");
//	}
//	printf("]\n\n");
//}

//      o que mostrar, dimens�o, no de casas decimais
void showMat(const bool* show, const unsigned int dim) {
	printf("\n[");

	for (unsigned int j = 0; j < dim; j++) {
		for (unsigned int i = 0; i < dim; i++) {
			printf("%d ", show[j + i * dim]);
		}
		printf("\n ");
	}
	printf("]\n\n");
}

void showMatRI(const complex_type* show, const unsigned int dim) {
	printf("\nRe=\n[");

	for (unsigned int j = 0; j < dim; j++) {
		for (unsigned int i = 0; i < dim; i++) {
			printf("%8.4f ", show[j + i * dim].x);
		}
		printf("\n ");
	}
	printf("]\n\n");

	printf("\nIm=\n[");

	for (unsigned int j = 0; j < dim; j++) {
		for (unsigned int i = 0; i < dim; i++) {
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
	//std::cout<< "bseMVA  = " << sistema.baseMVA << std::endl;

	std::cout<< "Pload:" << std::endl;
	showVec(barra.Pload, sistema.nB, 4);
	std::cout<< "Qload:" << std::endl;
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
	//std::cout<< "r:" << std::endl;
	//showVecR(ramo.z, nl, 4);
	//std::cout<< "x:" << std::endl;
	//showVecI(ramo.z, nl, 4);
	//std::cout<< "tap (Re):" << std::endl;
	//showVecR(ramo.tap, nl, 1);
	//std::cout<< "tap (Im):" << std::endl;
	//showVecI(ramo.tap, nl, 1);
	std::cout<< "tap :" << std::endl;
	showVec(ramo.tap, sistema.nL, 4);
	std::cout<< "barrasPV :" << std::endl;
	showVec(sistema.barrasPV, sistema.nPV, 1);
	std::cout<< "barrasPQ :" << std::endl;
	showVec(sistema.barrasPQ, sistema.nPQ, 1);
	std::cout<< "barraVO : " << sistema.barraVO << std::endl;


//	std::cout<< "Ybus (Re):" << std::endl;
//	showMatR(sistema.Y, sistema.nB, 1);
//
//	std::cout<< "Ybus (Im):" << std::endl;
//	showMatI(sistema.Y, sistema.nB, 1);

	std::cout<< "Ybus:" << std::endl;
	if (global::metodo != esparso) {
	showMat(sistema.Y, sistema.nPQ+sistema.nPV+1);

	showMatRI(sistema.Y, sistema.nPQ + sistema.nPV + 1);
	}
	else{
		std::cout<< *sistema.spY << std::endl;
	}
}
