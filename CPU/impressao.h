#pragma once
#include <stdio.h>
#include <fstream>
#include <chrono>

#include "global.h"
#include "idx.h"
#include "estruturas.h"

#include "benchmarks.h"

#include <iterator>


void impressao(sistemaType &sistema, barraType &barra, ramoType &ramo, iterativoType &iterativo);

void impressao2(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo);

void benchmarkModePrint(iterativoType& iterPon, /*std::chrono::duration<double, std::milli>*/double& duracao);

void benchmarksPrint(iterativoType& iterPon);

void benchmarksPrintFile(iterativoType& iterPon);