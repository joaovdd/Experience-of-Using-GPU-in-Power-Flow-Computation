#pragma once
#include "../global.h"
#include "sistema.h"

struct barraType {
	int* id; 
	float_type* V;
	float_type* theta;

	float_type* Pliq;
	float_type* Qliq;
	float_type* Pload;
	float_type* Qload;
	float_type* Pg;
	float_type* Qg;

	float_type* Vbase;

	float_type* gsh;
	float_type* bsh;

};

void initBus(sistemaType &sistema, barraType &barra);
void finBus(barraType &barra);