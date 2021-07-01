#pragma once
#include <string>
#include <type_traits>

//#include "externo/cuComplex.h"
#include <cuComplex.h>

#include "estrutura/metodo.h"

#include "opcoesDeCompilacao.h"

namespace global {
	extern float_type v_inicial;
	extern float_type theta_inicial;
	extern unsigned int no_max_iter;
	extern float_type tol;
	extern float_type tol_limInjReat;
	extern std::string arq_entrada;
	extern bool verbose_mode;
	extern bool lim_inj_reat;
	extern bool laconic_mode;
	extern bool openmp;
	extern metodoType metodo;
	extern bool temporizador;
	extern output_benchmarkType output_benchmark; // output_benchmarkType definition in metodo.h
	extern bool output_ans;
	extern bool output_processo_iterativo;
}

static inline cuDoubleComplex _cuAdd(cuDoubleComplex x, cuDoubleComplex y) { return cuCadd(x, y); }
static inline cuDoubleComplex _cuSub(cuDoubleComplex x, cuDoubleComplex y) { return cuCsub(x, y); }
static inline cuDoubleComplex _cuDiv(cuDoubleComplex x, cuDoubleComplex y) { return cuCdiv(x, y); }
static inline cuDoubleComplex _cuMul(cuDoubleComplex x, cuDoubleComplex y) { return cuCmul(x, y); }
static inline cuDoubleComplex _cuCon(cuDoubleComplex x) { return cuConj(x); }
static inline double _cuAbs(cuDoubleComplex x) { return cuCabs(x); }
static inline double _cuReal(cuDoubleComplex x) { return cuCreal(x); }
static inline double _cuImag(cuDoubleComplex x) { return cuCimag(x); }

// static inline cuDoubleComplex _mkComplex(double x, double y) { return make_cuDoubleComplex(x, y); }

static inline cuFloatComplex _cuAdd(cuFloatComplex x, cuFloatComplex y) { return cuCaddf(x, y); }
static inline cuFloatComplex _cuSub(cuFloatComplex x, cuFloatComplex y) { return cuCsubf(x, y); }
static inline cuFloatComplex _cuDiv(cuFloatComplex x, cuFloatComplex y) { return cuCdivf(x, y); }
static inline cuFloatComplex _cuMul(cuFloatComplex x, cuFloatComplex y) { return cuCmulf(x, y); }
static inline cuFloatComplex _cuCon(cuFloatComplex x) { return cuConjf(x); }
static inline float _cuAbs(cuFloatComplex x) { return cuCabsf(x); }
static inline float _cuReal(cuFloatComplex x) { return cuCrealf(x); }
static inline float _cuImag(cuFloatComplex x) { return cuCimagf(x); }

static inline complex_type _mkComplex(float_type x, float_type y) { return MK_COMPLEX_FUNCTION(x, y); }
