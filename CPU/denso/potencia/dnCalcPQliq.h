#pragma once
#include "../../global.h"
#include "../../estruturas.h"
#include "../../idx.h"
#include "../../phifBshf.h"

float_type calcP(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type calcP2(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type calcQ(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo);

float_type calcQ2(const unsigned int k, sistemaType* sistema, barraType* barra, ramoType* ramo);

void dnCalculePQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);

void dnCalculePQ_subsist2(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo);