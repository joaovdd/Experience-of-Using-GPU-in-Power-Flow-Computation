#include "dnCalcYbus.h"

void calcYbus(sistemaType &sistema, barraType &barra, ramoType &ramo) {
	for (int i = 1; i <= sistema.nL; i++) {
		complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); 
		aux = _cuDiv(_mkComplex(-1., 0.), aux); 
		
		complex_type aux2 = _cuDiv(_mkComplex(-1., 0.), ramo.tap[IDX1F(i)]); 
		aux2 = _cuDiv(aux2, ramo.z[IDX1F(i)]); 
		
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = _cuDiv(aux, ramo.z[IDX1F(i)]); 
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux2;

		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); 
		dAux = 1 / (dAux * dAux); 

		aux = _mkComplex(dAux, 0.); 
		aux = _cuDiv(aux, ramo.z[IDX1F(i)]); 
		aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); 
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

		aux = _cuDiv(_mkComplex(1., 0.), ramo.z[IDX1F(i)]);
		aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); 
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
	}

	for (int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])); 
	}
}