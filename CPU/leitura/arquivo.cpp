#include "arquivo.h"

void lerArquivoEAlocarMemoria(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	//std::string ext = global::arq_entrada.substr(global::arq_entrada.size() - 3, 3);
	
	auto posPonto = global::arq_entrada.find_last_of('.');
	std::string ext = global::arq_entrada.substr(posPonto + 1, global::arq_entrada.size() - posPonto);

	std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

	if (ext == "cdf") {
		//printf("CDF == %s", ext.c_str());
		lerTamanhos(global::arq_entrada, sistema); // para aloca��o din�mica das vari�veis

		initBranch(sistema, ramo);
		initSistema(sistema);
		initBus(sistema, barra);
		initIter(sistema, iterativo);

		if (readCDF(global::arq_entrada, sistema, barra, ramo)) { // l� dados do arquivo .CDF
			printf("Deu ruim...");
		}
		InitCsrPhi(sistema, ramo);
	}
	else if (ext == "cdfx") {
		lerTamanhosCDFX(global::arq_entrada, sistema); // para aloca��o din�mica das vari�veis

		initBranch(sistema, ramo);
		initSistema(sistema);
		initBus(sistema, barra);
		initIter(sistema, iterativo);

		if (readCDFX(global::arq_entrada, sistema, barra, ramo)) { // l� dados do arquivo .CDF
			printf("Deu ruim...");
		}
		InitCsrPhi(sistema, ramo);
	}
	else if (ext == "pwf") {
		printf("PWF == %s", ext.c_str());
		lerPWFEAlocarMemoria(global::arq_entrada, sistema, barra, ramo, iterativo);
		InitCsrPhi(sistema, ramo);
	}
	else {
		printf("nda == %s", ext.c_str());
	}
}
