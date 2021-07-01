#pragma once
#include <iostream>
#include <fstream>

#include "../idx.h"
#include "../estruturas.h"
#include "../global.h"
#include "id2i.h"

#include "config.h"

bool lerTamanhos(std::string cdfFile, sistemaType &sistema);

bool lerTamanhosCDFX(std::string cdfFile, sistemaType& sistema);

bool readCDF(std::string cdfFile, sistemaType& sistema, barraType& barra, ramoType& ramo);

bool readCDFX(std::string cdfFile, sistemaType& sistema, barraType& barra, ramoType& ramo);