#pragma once
#include "global.h"
#include "estruturas.h"

#include "denso/dnSolver.h"
#include "esparso/spSolver.h"

bool solver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);