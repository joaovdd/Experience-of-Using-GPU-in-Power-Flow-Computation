#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <iostream>

#include <chrono>
#include <omp.h>

#include "global.h"

#include "leitura/config.h"
#include "estruturas.h"
#include "dim.h"

#include "impressao.h"
#include "leitura.h"

#include "newtonRaphson.h"
#include "denso/dnPotencia.h"
#include "denso/dnCalcYbus.h"

#include "benchmarks.h"

// ATENÇÃO: limite de sistemas de 65.535 barras pelo uso de int
// Bug: barra swing do método esparso não pode ser a última!

//#define JACOBIANO_ESPARSO_STENCIL true // ativar apenas um dos modos de cálculo do jacobiano -> no arquivo global.h
//#define JACOBIANO_ESPARSO_PADRAO false

int main()
{
	//timerType timer(processoIterativo);
	loadFile(); // config
	sistemaType sistema;
	barraType barra;
	ramoType ramo;
	iterativoType iterativo;
	
	lerArquivoEAlocarMemoria(sistema, barra, ramo, iterativo); // data sistema

	if(global::openmp /*&& (global::metodo == denso)*/){Eigen::setNbThreads(0);}
	else{Eigen::setNbThreads(1);}

	// auto inicio = std::chrono::high_resolution_clock::now();

	{ timerType timer(geral);

		// calcula Ybus
		switch (global::metodo)
		{
		case denso:
		case esparsoSimples:
		{ BENCHMARK_ADMITANCIA
			calcYbus(sistema, barra, ramo);
			break;
		}

		case esparso:
		{ BENCHMARK_ADMITANCIA
			//calcYbusSp(sistema, barra, ramo);
			//calcYbusSp_eficinte(sistema, barra, ramo);
			calcYbusSp_Matpower(sistema, barra, ramo);
			//std::cout << *sistema.spY << std::endl;
			break;
		}

		default:
			printf("\n\n[ERRO] calcYbus: metodo inválido!\n\n");
			return -1;
			break;
		}

		if (global::verbose_mode) {
			printAll(sistema, barra, ramo);
		}

		switch (global::metodo) {
		case denso:
			if (global::laconic_mode) { printf("[DENSO]\n"); }
			else { printf("\n\n[DENSO]\n"); }
			break;
		case denso_LAPACKE:
			if (global::laconic_mode) { printf("[DENSO_LAPACK]\n"); }
			else { printf("\n\n[DENSO_LAPACK]\n"); }
			// printf("\n\n[DENSO_LAPACK]\n\n");
			break;
		case esparso:
			if (global::laconic_mode) { printf("[ESPARSO]\n"); }
			else { printf("\n\n[ESPARSO]\n"); }
			//printf("\n\n[ESPARSO]\n\n");
			break;
		case esparsoSimples:
			if (global::laconic_mode) { printf("[ESPARSO SIMPLES]\n"); }
			else { printf("\n\n[ESPARSO SIMPLES]\n"); }
			//printf("\n\n[ESPARSO SIMPLES]\n\n");
			break;
		default:
			std::cout << "[METODO] METODO INVALIDO LIDO DO ARQUIVO!\n" << std::endl;
		}

		nR(&sistema, &barra, &ramo, &iterativo);

		if (global::verbose_mode) {
			std::cout << "V:" << std::endl;
			showVec(barra.V, sistema.nB, 4);
			std::cout << "theta:" << std::endl;
			showVec(barra.theta, sistema.nB, 4);
		}

		{ BENCHMARK_FLUXO
			switch (global::metodo) {
			case denso:
			case denso_LAPACKE:
			case esparsoSimples:
				calcFlux(&sistema, &barra, &ramo);
				break;
			case esparso:
				SpCalcFlux(&sistema, &barra, &ramo);
				break;
			case hibridoA:
			case hibridoB:
			case nda:
			default:
				break;
			}
		}
	}

	//if (global::laconic_mode) {
	//	printf("%d iteracoes.\n%f ms.\n", iterativo.iteracao, std::get<2>(global::tracker.benchmarkTable[0]));
	//}
	//else {
	//	impressao2(sistema, barra, ramo, iterativo); 
	//}
	if (!global::laconic_mode &&
	    global::output_ans) {
		impressao2(sistema, barra, ramo, iterativo);
	}
	if (global::laconic_mode || 
		global::output_benchmark == output_benchmarkType::all || 
		global::output_benchmark == output_benchmarkType::file) {
		benchmarksPrintFile(iterativo);
	}
	if (!global::laconic_mode && 
		(global::output_benchmark == output_benchmarkType::screen ||
		 global::output_benchmark == output_benchmarkType::all    )) {
		benchmarksPrint(iterativo);
	}

 	finBranch(ramo);
	finSistema(sistema);
	finBus(barra);
	finIter(iterativo);

	if (!global::laconic_mode) {
		printf("\nExecucao concluida."); // Pressione qualquer tecla para sair.
	}

	return 0;
}
