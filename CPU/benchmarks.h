// Copyright (c) 2020, Jo�o Daibes

#define BENCHMARK_MODE true // compila programa com as rotinas para medi��o de tempo de execu��o

#if BENCHMARK_MODE

#define BENCHMARK_GERAL timerType timer(benchmarkType::geral);
#define BENCHMARK_ADMITANCIA timerType timer(benchmarkType::admitancia);
#define BENCHMARK_PROCESSOITERATIVO timerType timer(benchmarkType::processoIterativo);
#define BENCHMARK_CALCPQ timerType timer(benchmarkType::calcPQ, iterativo->iteracao);
#define BENCHMARK_JACOBIANO timerType timer(benchmarkType::jacobiano, iterativo->iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_FILL timerType timer(benchmarkType::jacobianoStencil_fill, iterativo->iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_REBUILD timerType timer(benchmarkType::jacobianoStencil_rebuild, iterativo->iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_BUILD timerType timer(benchmarkType::jacobianoStencil_build, iterativo->iteracao);
#define BENCHMARK_SISTEMALINEAR timerType timer(benchmarkType::sistemaLinear, iterativo->iteracao);
#define BENCHMARK_FLUXO timerType timer(benchmarkType::fluxo);

#else

#define BENCHMARK_GERAL  
#define BENCHMARK_ADMITANCIA  
#define BENCHMARK_PROCESSOITERATIVO  
#define BENCHMARK_CALCPQ  
#define BENCHMARK_JACOBIANO  
#define BENCHMARK_JACOBIANOSTENCIL_FILL  
#define BENCHMARK_JACOBIANOSTENCIL_REBUILD  
#define BENCHMARK_JACOBIANOSTENCIL_BUILD  
#define BENCHMARK_SISTEMALINEAR  
#define BENCHMARK_FLUXO

#endif

#ifndef BENCHMARKS
#define BENCHMARKS

#include <iostream>
#include <vector>
#include <tuple>
#include <chrono>

enum benchmarkType {
	geral,
	admitancia,
	processoIterativo,
	calcPQ,
	jacobiano,
	jacobianoStencil_fill,
	jacobianoStencil_rebuild,
	jacobianoStencil_build,
	sistemaLinear,
	fluxo
};

class trackerType {
public:
	//                 tipo de medida; itera��o; valor
	std::vector<std::tuple<benchmarkType, int, double>> benchmarkTable;

	void log(benchmarkType tipo, int iteracao, double duracao) {
		benchmarkTable.push_back(std::make_tuple(tipo, iteracao, duracao));
	}
};

namespace global {
	extern trackerType tracker;
}

class timerType {
public:
	timerType(benchmarkType tipo, int iteracao) {
		m_inicioTimepoint = std::chrono::high_resolution_clock::now();
		m_tipo = tipo;
		m_iteracao = iteracao;
	}
	timerType(benchmarkType tipo) {
		m_inicioTimepoint = std::chrono::high_resolution_clock::now();
		m_tipo = tipo;
	}
	~timerType() {
		stop();
	}
	void stop() {
		auto fimTimepoint = std::chrono::high_resolution_clock::now();

		std::chrono::duration<double, std::nano> duration = fimTimepoint - m_inicioTimepoint;
		global::tracker.log(m_tipo, m_iteracao, duration.count() / 1000000); 
	}
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> m_inicioTimepoint;
	benchmarkType m_tipo;
	int m_iteracao;
};
#endif