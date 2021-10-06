#pragma once
#include "../../global.h"
#include "../../estruturas.h"
#include "../../idx.h"
#include "../../phifBshf.h"

float_type P(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type Q(int k, int m, sistemaType* sistema, barraType* barra, ramoType* ramo);

void calcFlux(sistemaType* sistema, barraType* barra, ramoType* ramo);