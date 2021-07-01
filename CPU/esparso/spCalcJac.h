#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include "../leitura.h"
#include "../estruturas.h"
#include "../dim.h"
#include "../phifBshf.h"

//#include "lapacke.h"
#ifndef FLUMEN_GPU
    #include "../externo/Eigen/Sparse"
#endif
void Jstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void fillJstencil0based(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void spCalcJ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void spCalcJ_eficiente(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);