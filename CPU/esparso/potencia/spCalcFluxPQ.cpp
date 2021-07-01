#include "spCalcFluxPQ.h"

float_type SpP(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = phif(k, m, sistema, ramo); // defasagem do transformador
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; // theta_k para[i]
	return (barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * real(sistema->spY->coeff(IDX1F(k), IDX1F(m))) /*sistema->Y[IDX2F(k, m, sistema->nB)].x*/ - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * real(sistema->spY->coeff(IDX1F(k), IDX1F(m))) /*sistema->Y[IDX2F(k, m, sistema->nB)].x*/ * cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * imag(sistema->spY->coeff(IDX1F(k), IDX1F(m))) /*sistema->Y[IDX2F(k, m, sistema->nB)].y*/ * sin(aux));
}

float_type SpQ(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = phif(k, m, sistema, ramo); // defasagem do transformador
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; // theta_k para[i]
	return (-barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * (/*sistema->Y[IDX2F(k, m, sistema->nB)].y*/ imag(sistema->spY->coeff(IDX1F(k), IDX1F(m))) + bshf(k, m, sistema, ramo)) + barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * /*sistema->Y[IDX2F(k, m, sistema->nB)].y*/ imag(sistema->spY->coeff(IDX1F(k), IDX1F(m))) * cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * /*sistema->Y[IDX2F(k, m, sistema->nB)].x*/ real(sistema->spY->coeff(IDX1F(k), IDX1F(m))) * sin(aux));
}

void SpCalcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo) {
	// Calcula os valores de Pkm...
#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 1; i <= sistema->nL; i++) {
		ramo->Pdp[IDX1F(i)] = SpP(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
	}
	// Calcula os valores de Qkm...
#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 1; i <= sistema->nL; i++) {
		ramo->Qdp[IDX1F(i)] = SpQ(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
	}
#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 1; i <= sistema->nL; i++) {
		ramo->Ppd[IDX1F(i)] = SpP(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
	}
	// Calcula os valores de Qkm...
#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 1; i <= sistema->nL; i++) {
		ramo->Qpd[IDX1F(i)] = SpQ(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
	}
}