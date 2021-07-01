#pragma once
#include <cuComplex.h>

//#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include "sistInfo.h"

using namespace std;

//      o que mostrar, dimens�o, no de casas decimais
void showVec(float_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(unsigned int* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(unsigned short* show, unsigned short dim, unsigned short prec) {
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (unsigned short i = 1; i < dim - 1; i++) {
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVec(bool* show, unsigned short dim, unsigned short prec) {
	//std::cout.unsetf(std::ios::floatfield);
	//std::cout.setf(std::ios::fixed, std::ios::floatfield);
	//std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (unsigned short i = 1; i < dim - 1; i++) {
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVecR(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].x << endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].x << endl;
	}

	cout << ' ' << show[dim - 1].x << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showVecI(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y << endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].y << endl;
	}

	cout << ' ' << show[dim - 1].y << ']' << endl;
}

void showMatRI(const complex_type* show, const unsigned short dim) {
	printf("\nRe=\n[");

	for (unsigned short j = 0; j < dim; j++) {
		for (unsigned short i = 0; i < dim; i++) {
			printf("%8.4f ", show[j + i * dim].x);
		}
		printf("\n ");
	}
	printf("]\n\n");

	printf("\nIm=\n[");

	for (unsigned short j = 0; j < dim; j++) {
		for (unsigned short i = 0; i < dim; i++) {
			printf("%8.4f ", show[j + i * dim].y);
		}
		printf("\n ");
	}
	printf("]\n\n");
}


//      o que mostrar, dimens�o, no de casas decimais
void showVec(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].x << " + j*" << show[0].y << endl;

	for (unsigned int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].x << " + j*" << show[i].y << endl;
	}

	cout << ' ' << show[dim - 1].x << " + j*" << show[dim - 1].y << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMatI(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y;

	for (unsigned int j = 1; j < dim; j++){
    	for (unsigned int i = 1; i < dim; i++){
	    	cout << ' ' << show[j+i*dim].y;
	    }
		cout << endl;
	}
	cout << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMatR(complex_type* show, unsigned int dim, unsigned int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y;

	for (unsigned int j = 1; j < dim; j++){
    	for (unsigned int i = 1; i < dim; i++){
	    	cout << ' ' << show[j+i*dim].x;
	    }
		cout << endl;
	}
	cout << ']' << endl;
}

//      o que mostrar, dimens�o, no de casas decimais
void showMat(const complex_type* show, const unsigned int dim){
	printf("\n[");

	for (unsigned int j = 0; j < dim; j++){
    	for (unsigned int i = 0; i < dim; i++){
    		printf("(%7.3f + j*%7.3f) ", show[j+i*dim].x ,show[j+i*dim].y);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

//      o que mostrar, dimens�o, no de casas decimais
void showMat(const float_type* show, const unsigned int dim){
	printf("\n[");

	for (unsigned int j = 0; j < dim; j++){
    	for (unsigned int i = 0; i < dim; i++){
    		printf("%5.1f ", show[j+i*dim]);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

void printAll(sistema &sistema, barra &barra, ramo &ramo){
	cout << "nB  = " << sistema.nB << endl;
	cout << "nPQ = " << sistema.nPQ << endl;
	cout << "nPV = " << sistema.nPV << endl;
	cout << "nl  = " << sistema.nL << endl;

	cout << "Pload:" << endl;
	showVec(barra.Pload, sistema.nB, 1);
	cout << "Qload:" << endl;
	showVec(barra.Qload, sistema.nB, 1);
	cout << "Pg:" << endl;
	showVec(barra.Pg, sistema.nB, 1);
	cout << "Qg:" << endl;
	showVec(barra.Qg, sistema.nB, 1);
	cout << "VBase:" << endl;
	showVec(barra.Vbase, sistema.nB, 1);
	cout << "Vb:" << endl;
	showVec(barra.V, sistema.nB, 5);
	cout << "bsh:" << endl;
	showVec(barra.bsh, sistema.nB, 1);
	cout << "gsh:" << endl;
	showVec(barra.gsh, sistema.nB, 1);
	
	cout << "de:" << endl;
	showVec(ramo.de, sistema.nL, 1);
	cout << "para:" << endl;
	showVec(ramo.para, sistema.nL, 1);
	cout << "bsh:" << endl;
	showVec(ramo.bsh, sistema.nL, 1);
	
	cout << "z:" << endl;	
	showVec(ramo.z, sistema.nL, 4);
	//cout << "r:" << endl;
	//showVecR(ramo.z, nl, 4);
	//cout << "x:" << endl;
	//showVecI(ramo.z, nl, 4);
	//cout << "tap (Re):" << endl;
	//showVecR(ramo.tap, nl, 1);
	//cout << "tap (Im):" << endl;
	//showVecI(ramo.tap, nl, 1);
	cout << "tap :" << endl;
	showVec(ramo.tap, sistema.nL, 1);
	cout << "barrasPV :" << endl;
	showVec(sistema.barrasPV, sistema.nPV, 1);
	cout << "barrasPQ :" << endl;
	showVec(sistema.barrasPQ, sistema.nPQ, 1);
	cout << "barraVO : " << sistema.barraVO << endl;


//	cout << "Ybus (Re):" << endl;
//	showMatR(sistema.Y, sistema.nB, 1);
//
//	cout << "Ybus (Im):" << endl;
//	showMatI(sistema.Y, sistema.nB, 1);

	//cout << "Ybus:" << endl;
	//showMat(sistema.Y, sistema.nPQ+sistema.nPV+1);

	std::cout << "Ybus:" << std::endl;
	if (global::metodo != metodo::esparso && global::metodo != metodo::hibridoA && global::metodo != metodo::hibridoB) {
		showMat(sistema.Y, sistema.nPQ + sistema.nPV + 1);

		showMatRI(sistema.Y, sistema.nPQ + sistema.nPV + 1);
	}
	else {
		std::cout << *sistema.spY << std::endl;
	}
}
