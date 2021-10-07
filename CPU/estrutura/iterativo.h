#pragma once
#include "../global.h"
#include "sistema.h"

struct iterativoType {
	float_type* Pcalc; 
	float_type* Qcalc; 

	float_type* J; 
	int iteracao, noMax;

	float_type* g;

	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	bool* limQ; 
	            
	int* barrasPVlim; 
	int nPVlim; 
	int* barrasPQlim; 
	int nPQlim; 

	float_type* QliqLim;       
};

void initIter(sistemaType &sistema, iterativoType &iterativo);
void finIter(iterativoType &iterativo);