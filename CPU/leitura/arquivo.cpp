#include "arquivo.h"

void lerArquivoEAlocarMemoria(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	auto posPonto = global::arq_entrada.find_last_of('.');
	std::string ext = global::arq_entrada.substr(posPonto + 1, global::arq_entrada.size() - posPonto);

	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	if (ext == "cdf") {
		
		lerTamanhos(global::arq_entrada, sistema); 

		initBranch(sistema, ramo);
		initSistema(sistema);
		initBus(sistema, barra);
		initIter(sistema, iterativo);

		if (readCDF(global::arq_entrada, sistema, barra, ramo)) { 
			printf("Erro ao abrir o arquivo...");
		}
		InitCsrPhi(sistema, ramo);
	}
	else if (ext == "cdfx") {
		lerTamanhosCDFX(global::arq_entrada, sistema); 

		initBranch(sistema, ramo);
		initSistema(sistema);
		initBus(sistema, barra);
		initIter(sistema, iterativo);

		if (readCDFX(global::arq_entrada, sistema, barra, ramo)) { 
			printf("Erro ao abrir o arquivo...");
		}
		InitCsrPhi(sistema, ramo);
	}
	else if (ext == "pwf") {
		printf("PWF == %s", ext.c_str());
		lerPWFEAlocarMemoria(global::arq_entrada, sistema, barra, ramo, iterativo);
		InitCsrPhi(sistema, ramo);
	}
	else if (ext == "m") {
		
		matPowerDataType mpData = lerMatPowerEAlocarMemoria(global::arq_entrada, sistema, barra, ramo, iterativo);
		
		InitCsrPhi(sistema, ramo);
	}

	else {
		printf("nda == %s", ext.c_str());
	}
}
