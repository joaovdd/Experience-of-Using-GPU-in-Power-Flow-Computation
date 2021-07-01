#include <string>

#include "estrutura/metodo.h"

#include "opcoesDeCompilacao.h"

// leitura dos valores iniciais de magnitude e �ngulo da tens�o nodal,
// No m�ximo de itera��es e nome do arquivo a ser aberto

namespace global {
	float_type v_inicial = 1.0, theta_inicial = 0.0;
	unsigned int no_max_iter = 100;
	float_type tol = 0.0001;
	float_type tol_limInjReat = 0.0001; // (p.u.; 0.01 MVAr se sistema.baseMVA = 100 MVAr)
	std::string arq_entrada = "ieee14.cdf";
	bool verbose_mode = 0;
	bool lim_inj_reat = 0;
	bool laconic_mode = 0;
	bool openmp = 1;
	metodoType metodo;
	bool temporizador = 1;
	output_benchmarkType output_benchmark = output_benchmarkType::all; // output_benchmarkType definition in metodo.h
	bool output_ans = 1;
	bool output_processo_iterativo = 0;
}