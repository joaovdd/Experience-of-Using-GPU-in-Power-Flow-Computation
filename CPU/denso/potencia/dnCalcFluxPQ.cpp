#include "dnCalcFluxPQ.h"

float_type P(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = dnPhif(k, m, sistema, ramo); // defasagem do transformador
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; // theta_k para[i]
	return (barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * sistema->Y[IDX2F(k, m, sistema->nB)].x - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * sistema->Y[IDX2F(k, m, sistema->nB)].x*cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * sistema->Y[IDX2F(k, m, sistema->nB)].y*sin(aux));
}

float_type Q(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo) {
	float_type aux = 0;
	aux = dnPhif(k, m, sistema, ramo); // defasagem do transformador
	aux += barra->theta[IDX1F(k)] - barra->theta[IDX1F(m)]; // theta_k para[i]
	return (-barra->V[IDX1F(k)] * barra->V[IDX1F(k)] * (sistema->Y[IDX2F(k, m, sistema->nB)].y + bshf(k, m, sistema, ramo)) + barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * sistema->Y[IDX2F(k, m, sistema->nB)].y*cos(aux) - barra->V[IDX1F(k)] * barra->V[IDX1F(m)] * sistema->Y[IDX2F(k, m, sistema->nB)].x*sin(aux));
}

void calcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo) {
	// Calcula os valores de Pkm...
	#pragma omp parallel for if (global::openmp)
		for (unsigned int i = 1; i <= sistema->nL; i++) {
			ramo->Pdp[IDX1F(i)] = P(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
		}
	// Calcula os valores de Qkm...
	#pragma omp parallel for if (global::openmp)
		for (unsigned int i = 1; i <= sistema->nL; i++) {
			ramo->Qdp[IDX1F(i)] = Q(ramo->de[IDX1F(i)], ramo->para[IDX1F(i)], sistema, barra, ramo);
		}
	#pragma omp parallel for if (global::openmp)
		for (unsigned int i = 1; i <= sistema->nL; i++) {
			ramo->Ppd[IDX1F(i)] = P(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
		}
	// Calcula os valores de Qkm...
	#pragma omp parallel for if (global::openmp)
		for (unsigned int i = 1; i <= sistema->nL; i++) {
			ramo->Qpd[IDX1F(i)] = Q(ramo->para[IDX1F(i)], ramo->de[IDX1F(i)], sistema, barra, ramo);
		}
}
