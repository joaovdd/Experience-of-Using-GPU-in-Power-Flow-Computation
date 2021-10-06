#pragma once
#include "../../global.h"
#include "../../estruturas.h"
#include "../../idx.h"
#include "../../phifBshf.h"

float_type SpP(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type SpQ(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

void SpCalcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo);