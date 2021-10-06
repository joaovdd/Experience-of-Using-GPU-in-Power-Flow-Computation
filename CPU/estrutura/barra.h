#pragma once
#include "../global.h"
#include "sistema.h"

// indexador é o numero da barra (via IDX1F)
struct barraType {
	int* id; // bus identifier number
	float_type* V;//[NB];
	float_type* theta;//[NB];

	float_type* Pliq;//[NB];
	float_type* Qliq;//[NB];
	float_type* Pload;//[NB];
	float_type* Qload;//[NB];
	float_type* Pg;//[NB];
	float_type* Qg;//[NB];

	float_type* Vbase;//[NB];

	float_type* gsh;//[NB];
	float_type* bsh;//[NB];

	//float_type* phi;
};

void initBus(sistemaType &sistema, barraType &barra);
void finBus(barraType &barra);