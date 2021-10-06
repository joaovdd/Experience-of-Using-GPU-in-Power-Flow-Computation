#pragma once
#include "../global.h"
#include "sistema.h"

struct iterativoType {
	float_type* Pcalc; //[NB];
	float_type* Qcalc; //[NB];

	float_type* J; //[(NB-1+NPQ)*(NB-1+NPQ)];
	int iteracao, noMax;

	float_type* g;

	//float_type tol;

	// limite de injeção de reativos
	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	//limite de injeção de reativos
	bool* limQ; // FUTURO?:dois bits: o 1o é 0 se o limite atingido é o superior, 1 se o inferior. O segundo é
	            // 0 se o limite de injeção de reativos não foi atingido, 1 caso contrário
				// nPV entradas
	int* barrasPVlim; //[NPV]; // vetor contendo as barras de tipo PV após a análise de limite de reativos
	int nPVlim; // número de barras de tipo PV após a análise de limite de reativos
	int* barrasPQlim; //[NPV]; // vetor contendo as barras de tipo PV após a análise de limite de reativos
	int nPQlim; // número de barras de tipo PV após a análise de limite de reativos

	float_type* QliqLim;       // atualização do limite do valor de Qliq para as novas barras de tipo PQ
};

void initIter(sistemaType &sistema, iterativoType &iterativo);
void finIter(iterativoType &iterativo);