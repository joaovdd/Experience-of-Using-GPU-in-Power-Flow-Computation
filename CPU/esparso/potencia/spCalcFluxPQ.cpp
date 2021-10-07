#include "spCalcFluxPQ.h"

float_type SpP(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = phif(k, m, sistema, ramo); 
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; 
	return (barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * real(sistema->spY->coeff(IDX1F(k), IDX1F(m)))  - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * real(sistema->spY->coeff(IDX1F(k), IDX1F(m)))  * cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * imag(sistema->spY->coeff(IDX1F(k), IDX1F(m)))  * sin(aux));
}

float_type SpQ(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = phif(k, m, sistema, ramo); 
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; 
	return (-barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * ( imag(sistema->spY->coeff(IDX1F(k), IDX1F(m))) + bshf(k, m, sistema, ramo)) + barra->V[IDX1F(k)] * barra->V[IDX1F(m)] *  imag(sistema->spY->coeff(IDX1F(k), IDX1F(m))) * cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] *  real(sistema->spY->coeff(IDX1F(k), IDX1F(m))) * sin(aux));
}

void SpCalcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo) {
#pragma omp parallel for if (global::openmp)
	for (int i = 1; i <= sistema->nL; i++) {
		ramo->Pdp[IDX1F(i)] = SpP(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
	}

#pragma omp parallel for if (global::openmp)
	for (int i = 1; i <= sistema->nL; i++) {
		ramo->Qdp[IDX1F(i)] = SpQ(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
	}
#pragma omp parallel for if (global::openmp)
	for (int i = 1; i <= sistema->nL; i++) {
		ramo->Ppd[IDX1F(i)] = SpP(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
	}

#pragma omp parallel for if (global::openmp)
	for (int i = 1; i <= sistema->nL; i++) {
		ramo->Qpd[IDX1F(i)] = SpQ(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
	}
}