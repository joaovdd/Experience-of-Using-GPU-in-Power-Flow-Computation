#pragma once
//#include "externo/cuComplex.h"
#include <cuComplex.h>

// se desligado, o calculo do jacobiano ser� esparso ordin�rio.
// se ligado, ser�o utilizadas rotinas Stencil, otimizadas para GPU.
#define JACOBIANO_ESPARSO_STENCIL false // ativar apenas um dos modos de c�lculo do jacobiano

// defina como verdadeiro para que seja usado double,
// caso contr�rio ser� utilizado float 
#define DOUBLE_MODE true

// antigo

//// defini��o do tipo de ponto flutuante
//typedef double/*double*/ float_type;
//
//// defini��o do tipo de n�mero complexo
//typedef /*cuFloatComplex*/cuDoubleComplex complex_type;
//
////
//#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
////#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)

// *****************************************************************************************
// N�o alterar a partir daqui***************************************************************
// *****************************************************************************************

// defini��o do tipo de ponto flutuante
#if DOUBLE_MODE
typedef double float_type;
#else
typedef float float_type;
#endif

// defini��o do tipo de n�mero complexo
#if DOUBLE_MODE
typedef cuDoubleComplex complex_type;
#else
typedef cuFloatComplex complex_type;
#endif

// defini��o do constritor complexo
#if DOUBLE_MODE
#define MK_COMPLEX_FUNCTION(x, y) make_cuDoubleComplex(x, y) 
#else
#define MK_COMPLEX_FUNCTION(x, y) make_cuFloatComplex(x, y)
#endif