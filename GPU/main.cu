#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "opcoesDeCompilacao_2.h"

#include "helper.cuh"

#include <cuComplex.h>
#include <math.h>

#include <chrono>

#include "benchmarks.h"

#include "sistInfo.h"
#include "sistema.cuh"

#include "PQcalc.cuh"
#include "Jacobiano.cuh"
#include "dim.cuh"
#include "newtonRaphson.cuh"
#include "impressao.h"
#include "arquivo.h"

int main(void)
{
	loadFile();
	sistema h_sistema, sistPon, * d_sistema = NULL;
	barra h_barra, barraPon, * d_barra = NULL;
	ramo h_ramo, ramoPon, * d_ramo = NULL;
	iterativo h_iterativo, iterPon, * d_iterativo = NULL;

	lerArquivoEAlocarMemoria(h_sistema, h_barra, h_ramo, h_iterativo);

	switch (global::metodo) {
	case metodo::denso:
		if (global::laconic_mode) { printf("[DENSO]\n\n"); }
		else { printf("\n\n[DENSO]\n\n"); }
		break;
	case metodo::hibridoA:
		if (global::laconic_mode) { printf("[HIBRIDO A]\n\n"); }
		else { printf("\n\n[HIBRIDO A]\n\n"); }
		break;
	case metodo::hibridoB:
		if (global::CPUsolverFlg) {
			if (global::laconic_mode) { printf("[HIBRIDO C]\n\n"); }
			else { printf("\n\n[HIBRIDO C]\n\n"); }
		}
		else {
			if (global::laconic_mode) { printf("[HIBRIDO B]\n\n"); }
			else { printf("\n\n[HIBRIDO B]\n\n"); }
		}
		break;
	case metodo::esparso:
		if (global::laconic_mode) { printf("[ESPARSO]\n\n"); }
		else { printf("\n\n[ESPARSO]\n\n"); }
		break;
	default:
		std::cout << "[METODO] METODO INVALIDO LIDO DO ARQUIVO!\n" << std::endl;
	}

	constexpr int nStreams = 4;
	cudaStream_t streams[nStreams] = { nullptr };

	cudaDeviceProp deviceProp;
	deviceProp = initGPU(); 

	cudaEventCreate(&global::start);
	cudaEventCreate(&global::stop);
	{ cudaEventRecord(global::start); 
		{ BENCHMARK_ADMITANCIA
			
			switch (global::metodoDeCalculoDeYbus) {
			case metodoDeCalculoDeYbus::dnCPU:
				calcYbus(h_sistema, h_barra, h_ramo);
				break;
			case metodoDeCalculoDeYbus::spCPU:
				if (global::isMatpower)
					calcYbusSp_Matpower(h_sistema, h_barra, h_ramo);
				else
					calcYbusSp_eficinte(h_sistema, h_barra, h_ramo);
				break;
			default:
				printf("\n\n[ERRO] calcYbus: metodo inv??lido!\n\n");
				return -1;
				break;
			}
		}

		if (global::verbose_mode) {
			printAll(h_sistema, h_barra, h_ramo);
		}

		{ BENCHMARK_INITGPU
			
			checkCudaErrors(cudaMalloc(&d_sistema, sizeof(sistema)));
			checkCudaErrors(cudaMalloc(&d_barra, sizeof(barra)));
			checkCudaErrors(cudaMalloc(&d_ramo, sizeof(ramo)));
			checkCudaErrors(cudaMalloc(&d_iterativo, sizeof(iterativo)));

			if (global::streams) {
				if (global::metodo == metodo::esparso || global::metodo == metodo::hibridoB) {
					for (int i = 0; i < nStreams; i++)
					{
						checkCudaErrors(cudaStreamCreate(&streams[i]));
					}
				}
				else if (global::metodo == metodo::denso) {
					
				}
			}

			checkCudaErrors(cudaGetLastError());

			sistPon = d_initSistema(h_sistema, d_sistema);
			ramoPon = d_initRamo(h_sistema, h_ramo, d_ramo);
			barraPon = d_initBarra(h_sistema, d_barra);
			iterPon = d_initIter(h_sistema, d_iterativo);
			
			checkCudaErrors(cudaGetLastError());

			checkCudaErrors(cudaDeviceSynchronize());
		}

		{ BENCHMARK_CUDAMEMCPY
			
			sistemacpyH2D(h_sistema, d_sistema, sistPon);
			ramocpyH2D(h_sistema, h_ramo, d_ramo, ramoPon);
			barracpyH2D(h_sistema, h_barra, d_barra, barraPon);
			itercpyH2D(h_sistema, h_iterativo, d_iterativo, iterPon);
		}

		if (global::verbose_mode) {
			printf("V =\n");
			d_showVecf(barraPon.V, h_sistema.nB);
		}

		switch (global::metodo) {
		case metodo::hibridoA:
			d_criarYesparso(sistPon);
			
			break;
		default:
			break;
		}

		mNR(h_sistema, h_barra, h_ramo, h_iterativo, d_sistema, d_barra, d_ramo, d_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp, streams);

		{ BENCHMARK_CUDAMEMCPY
			cudaMemcpy(h_barra.V, barraPon.V, sizeof(float_type) * (sistPon.nB), cudaMemcpyDeviceToHost);

			cudaMemcpy(h_barra.theta, barraPon.theta, sizeof(float_type) * (sistPon.nB), cudaMemcpyDeviceToHost);
			BENCHMARK_SYNC
		}

		{ BENCHMARK_FLUXO
			switch (global::metodo) {
			case metodo::denso:
				
				calcFluxf_ef(h_sistema, h_barra, h_ramo, sistPon, barraPon, ramoPon, deviceProp); 
				break;
			case metodo::hibridoA:
			case metodo::hibridoB:
			case metodo::esparso:
				calcFluxf_Eficiente_Sp(h_sistema, h_barra, h_ramo, h_iterativo, sistPon, barraPon, ramoPon, iterPon, deviceProp, streams);
				break;
			default:
				std::cout << "[calcFluxf] METODO INVALIDO LIDO DO ARQUIVO!\n" << std::endl;
			}
			BENCHMARK_SYNC
		}

		{ BENCHMARK_CUDAMEMCPY
			checkCudaErrors(cudaMemcpy(h_ramo.Ppd, ramoPon.Ppd, sistPon.nL * sizeof(float_type), cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_ramo.Pdp, ramoPon.Pdp, sistPon.nL * sizeof(float_type), cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_ramo.Qpd, ramoPon.Qpd, sistPon.nL * sizeof(float_type), cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_ramo.Qdp, ramoPon.Qdp, sistPon.nL * sizeof(float_type), cudaMemcpyDeviceToHost));

			checkCudaErrors(cudaMemcpy(h_iterativo.Pcalc, iterPon.Pcalc, sistPon.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(h_iterativo.Qcalc, iterPon.Qcalc, sistPon.nB * sizeof(float_type), cudaMemcpyDeviceToHost));

			checkCudaErrors(cudaDeviceSynchronize()); 
		}

		cudaEventRecord(global::stop);
	}

	if (!global::laconic_mode &&
		global::output_ans) {
		impressao2(h_sistema, h_barra, h_ramo, h_iterativo);
	}
	if (global::laconic_mode ||
		global::output_benchmark == output_benchmarkType::all ||
		global::output_benchmark == output_benchmarkType::file) {
		benchmarksPrintFile(iterPon);
	}
	if (!global::laconic_mode &&
		(global::output_benchmark == output_benchmarkType::screen ||
			global::output_benchmark == output_benchmarkType::all)) {
		benchmarksPrint(iterPon);
	}

	for (int i = 0; i < nStreams; i++) {
		if (streams[i]) { checkCudaErrors(cudaStreamDestroy(streams[i])); }
	}

	d_finSistema (d_sistema, sistPon);
	d_finBarra (d_barra, barraPon);
	d_finRamo (d_ramo, ramoPon);
	d_finIter (d_iterativo, iterPon);

	finSistema(h_sistema);
	finBranch(h_ramo);
	finBus(h_barra);
	finIter(h_iterativo);
}
