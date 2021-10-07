#pragma once
#include <math.h>

#include "leitura.h"
#include "estruturas.h"
#include "dim.h"

#include "limInjReat.h"
#include "denso/dnPotencia.h"
#include "esparso/spPotencia.h"
#include "denso/dnCalcJac.h"
#include "esparso/spCalcJac.h"
#include "solver.h"

#include <stdlib.h>
#include <stdio.h>

#ifndef FLUMEN_GPU
    #include "externo/Eigen/Dense"
    #include "externo/Eigen/Sparse"
#endif

#include "benchmarks.h"

void nR(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);