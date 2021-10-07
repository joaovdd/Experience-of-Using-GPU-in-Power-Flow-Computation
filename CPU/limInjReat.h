#pragma once
#include "estruturas.h"
#include "dim.h"

#include "esparso.h"

#include "benchmarks.h"

void geraVetPV(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void geraVetPQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void do_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

bool undo_chkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void calcResLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void initLim(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void attVOlim(sistemaType* sistema, barraType* barra, iterativoType* iterativo);

void printAlLim(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo);