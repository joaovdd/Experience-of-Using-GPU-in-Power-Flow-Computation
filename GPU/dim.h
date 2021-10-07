#pragma once
#include <cuComplex.h>

#include <fstream>
#include <iostream>
#include <math.h>
#include "sistInfo.h"

using namespace std;

void showVec(float_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

void showVec(int* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0] << endl;

	for (int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

void showVec(bool* show, int dim, int prec) {
	cout << '[' << show[0] << endl;

	for (int i = 1; i < dim - 1; i++) {
		cout << ' ' << show[i] << endl;
	}

	cout << ' ' << show[dim - 1] << ']' << endl;
}

void showVecR(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].x << endl;

	for (int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].x << endl;
	}

	cout << ' ' << show[dim - 1].x << ']' << endl;
}

void showVecI(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y << endl;

	for (int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].y << endl;
	}

	cout << ' ' << show[dim - 1].y << ']' << endl;
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

void showVec(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].x << " + j*" << show[0].y << endl;

	for (int i = 1; i < dim - 1; i++){
		cout << ' ' << show[i].x << " + j*" << show[i].y << endl;
	}

	cout << ' ' << show[dim - 1].x << " + j*" << show[dim - 1].y << ']' << endl;
}

void showMatI(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y;

	for (int j = 1; j < dim; j++){
    	for (int i = 1; i < dim; i++){
	    	cout << ' ' << show[j+i*dim].y;
	    }
		cout << endl;
	}
	cout << ']' << endl;
}

void showMatR(complex_type* show, int dim, int prec){
	std::cout.unsetf(std::ios::floatfield);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(prec);

	cout << '[' << show[0].y;

	for (int j = 1; j < dim; j++){
    	for (int i = 1; i < dim; i++){
	    	cout << ' ' << show[j+i*dim].x;
	    }
		cout << endl;
	}
	cout << ']' << endl;
}

void showMat(const complex_type* show, const int dim){
	printf("\n[");

	for (int j = 0; j < dim; j++){
    	for (int i = 0; i < dim; i++){
    		printf("(%7.3f + j*%7.3f) ", show[j+i*dim].x ,show[j+i*dim].y);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

void showMat(const float_type* show, const int dim){
	printf("\n[");

	for (int j = 0; j < dim; j++){
    	for (int i = 0; i < dim; i++){
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

	std::cout << "Pload:" << std::endl;
	showVec(barra.Pload, sistema.nB, 4);
	std::cout << "Qload:" << std::endl;
	showVec(barra.Qload, sistema.nB, 4);
	std::cout << "Pg:" << std::endl;
	showVec(barra.Pg, sistema.nB, 4);
	std::cout << "Qg:" << std::endl;
	showVec(barra.Qg, sistema.nB, 4);
	std::cout<< "Pliq:" << std::endl;
	showVec(barra.Pliq, sistema.nB, 4);
	std::cout<< "Qliq:" << std::endl;
	showVec(barra.Qliq, sistema.nB, 4);
	std::cout << "VBase:" << std::endl;
	showVec(barra.Vbase, sistema.nB, 4);
	std::cout << "V:" << std::endl;
	showVec(barra.V, sistema.nB, 4);
	std::cout<< "theta:" << std::endl;
	showVec(barra.theta, sistema.nB, 4);
	std::cout << "bsh:" << std::endl;
	showVec(barra.bsh, sistema.nB, 4);
	std::cout << "gsh:" << std::endl;
	showVec(barra.gsh, sistema.nB, 4);

	std::cout << "de:" << std::endl;
	showVec(ramo.de, sistema.nL, 4);
	std::cout << "para:" << std::endl;
	showVec(ramo.para, sistema.nL, 4);
	std::cout << "bsh(ramo):" << std::endl;
	showVec(ramo.bsh, sistema.nL, 4);
	std::cout<< "phi(trafo):" << std::endl;
	showVec(ramo.phi, sistema.nL, 4);

	std::cout << "z:" << std::endl;	
	showVec(ramo.z, sistema.nL, 4);

	std::cout << "tap :" << std::endl;
	showVec(ramo.tap, sistema.nL, 4);
	std::cout << "barrasPV :" << std::endl;
	showVec(sistema.barrasPV, sistema.nPV, 4);
	std::cout << "barrasPQ :" << std::endl;
	showVec(sistema.barrasPQ, sistema.nPQ, 4);
	std::cout << "barraVO : " << sistema.barraVO << std::endl;

	std::cout << "Ybus:" << std::endl;
	if (global::metodo != metodo::esparso && global::metodo != metodo::hibridoA && global::metodo != metodo::hibridoB) {
		showMat(sistema.Y, sistema.nPQ + sistema.nPV + 1);

		showMatRI(sistema.Y, sistema.nPQ + sistema.nPV + 1);
	}
	else {
		std::cout << *sistema.spY << std::endl;
	}
}
