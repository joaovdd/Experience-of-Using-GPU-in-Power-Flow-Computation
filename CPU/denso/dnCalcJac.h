#pragma once
#include "../estruturas.h"
#include "../dim.h"
#include "../global.h"
#include "../phifBshf.h"

void calcHlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void calcLlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void calcMlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void calcNlimf(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void calcJacLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);