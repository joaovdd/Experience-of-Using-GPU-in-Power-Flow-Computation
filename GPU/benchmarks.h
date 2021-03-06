

#define BENCHMARK_MODE true 

#if BENCHMARK_MODE

#define BENCHMARK_GERAL timerType timer(benchmarkType::geral);
#define BENCHMARK_ADMITANCIA timerType timer(benchmarkType::admitancia);
#define BENCHMARK_INITGPU checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::init_GPU);
#define BENCHMARK_PROCESSOITERATIVO checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::processoIterativo);
#define BENCHMARK_INIT_LIB checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::init_lib, iterPon.iteracao);
#define BENCHMARK_CUDAMEMCPY checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::cudaMemcpy, iterPon.iteracao);
#define BENCHMARK_CALCPQ checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::calcPQ, iterPon.iteracao);
#define BENCHMARK_JACOBIANO checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::jacobiano, iterPon.iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_FILL checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::jacobianoStencil_fill, iterPon.iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_REBUILD checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::jacobianoStencil_rebuild, iterPon.iteracao);
#define BENCHMARK_JACOBIANOSTENCIL_BUILD checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::jacobianoStencil_build, iterPon.iteracao);
#define BENCHMARK_SISTEMALINEAR checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::sistemaLinear, iterPon.iteracao);
#define BENCHMARK_FLUXO checkCudaErrors(cudaDeviceSynchronize()); timerType timer(benchmarkType::fluxo);

#define AVISO_ERRO_TEMPO_TOTAL " "
#define AVISO_ERRO_TEMPO_TOTAL_FILE " "

#define BENCHMARK_SYNC checkCudaErrors(cudaDeviceSynchronize());

#else

#define BENCHMARK_GERAL timerType timer(benchmarkType::geral);  
#define BENCHMARK_ADMITANCIA
#define BENCHMARK_INITGPU
#define BENCHMARK_PROCESSOITERATIVO
#define BENCHMARK_INIT_LIB
#define BENCHMARK_CUDAMEMCPY
#define BENCHMARK_CALCPQ
#define BENCHMARK_JACOBIANO
#define BENCHMARK_JACOBIANOSTENCIL_FILL
#define BENCHMARK_JACOBIANOSTENCIL_REBUILD
#define BENCHMARK_JACOBIANOSTENCIL_BUILD
#define BENCHMARK_SISTEMALINEAR
#define BENCHMARK_FLUXO

#define AVISO_ERRO_TEMPO_TOTAL ""
#define AVISO_ERRO_TEMPO_TOTAL_FILE ""

#define BENCHMARK_SYNC

#endif

#ifndef BENCHMARKS
#define BENCHMARKS

#include <iostream>
#include <vector>
#include <tuple>
#include <chrono>

#include <limits> 

enum class benchmarkType {
	geral,
	admitancia,
	processoIterativo,
	init_GPU,
	init_lib,
	cudaMemcpy,
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

	std::vector<std::tuple<benchmarkType, int, double>> benchmarkTable;

	void log(benchmarkType tipo, int iteracao, double duracao) {
		benchmarkTable.push_back(std::make_tuple(tipo, iteracao, duracao));
	}
};

namespace global {
	trackerType tracker;
	cudaEvent_t start, stop;
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