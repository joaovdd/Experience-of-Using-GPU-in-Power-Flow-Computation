#pragma once
#include "../idx.h"
#include "../estrutura/sistema.h"
#include "../estrutura/barra.h"
#include "../estrutura/ramo.h"

void calcYbusSp(sistemaType &sistema, barraType &barra, ramoType &ramo);

void calcYbusSp_eficinte(sistemaType& sistema, barraType& barra, ramoType& ramo);

void calcYbusSp_Matpower(sistemaType& sistema, barraType& barra, ramoType& ramo);