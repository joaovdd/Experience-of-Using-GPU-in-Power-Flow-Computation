#include <string>

#include "estrutura/metodo.h"

#include "opcoesDeCompilacao.h"

namespace global {
	float_type v_inicial = 1.0, theta_inicial = 0.0;
	int no_max_iter = 100;
	float_type tol = 0.0001;
	float_type tol_limInjReat = 0.0001; 
	std::string arq_entrada = "ieee14.cdf";
	bool verbose_mode = 0;
	bool lim_inj_reat = 0;
	bool laconic_mode = 0;
	bool openmp = 1;
	metodoType metodo;
	bool temporizador = 1;
	output_benchmarkType output_benchmark = output_benchmarkType::all; 
	bool output_ans = 1;
	bool output_processo_iterativo = 0;
}