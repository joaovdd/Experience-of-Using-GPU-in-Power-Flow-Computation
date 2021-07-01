#pragma once
#include "../../global.h"
#include "../../estruturas.h"
#include "../../idx.h"
#include "../../phifBshf.h"

float_type SpP(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type SpQ(unsigned int k, unsigned int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

void SpCalcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo);