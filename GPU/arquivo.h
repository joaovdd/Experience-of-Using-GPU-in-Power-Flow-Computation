#pragma once
#include "m.h"

bool lerTamanhosCDFX(std::string cdfFile, sistema& sistema) {
	std::string line;
	std::fstream CDF;
	CDF.open(cdfFile.c_str(), std::ios::in);
	if (CDF.is_open())
	{
		sistema.nB = 0;
		sistema.nPQ = 0;
		sistema.nPV = 0;
		sistema.barraVO = 0;
		sistema.baseMVA = 0;

		std::getline(CDF, line); std::getline(CDF, line); 
		while (std::getline(CDF, line)) 
		{
			if (line.find("-999") == std::string::npos) {
				sistema.nB++;	
				switch (atoi(line.substr(25, 2).c_str())) { 
				case 0: 
					sistema.nPQ++;
					break;
				case 1: 
					sistema.nPQ++;
					break;
				case 2: 
					sistema.nPV++;
					break;
				case 3:
					if (sistema.barraVO == 0) {
					}
					else {
						std::cout << "ERRO: mais de uma barrra swing definida!" << std::endl;
						return 1;
					}
					break;
				default:
					std::cout << "ERRO: valor de tipo de barra inválido!" << std::endl;
					return 1;
				}
			}
			else
				break; 
		}
		std::getline(CDF, line); 
		sistema.nL = 0;
		std::vector<int> de, para;

		while (std::getline(CDF, line)) 
		{
			if (line.find("-999") == std::string::npos) {
				int auxde = atoi(line.substr(0, 5).c_str()); 
				int auxpara = atoi(line.substr(6, 5).c_str()); 
				bool flgRamoNovo = 1;

				for (int i = 0; i < sistema.nL; i++) {
					if ((de[i] == auxde) && (para[i] == auxpara)) {
						
						flgRamoNovo = 0;
						break;
					}
				}
				if (flgRamoNovo) { 
					de.push_back(auxde);
					para.push_back(auxpara);
					sistema.nL++; 
				}
			}
			else
				break; 
		}
		if (sistema.nPQ + sistema.nPV + 1 != sistema.nB) { printf("não fechou!!"); }
		CDF.close();
	}
	else {
		std::cout << "ERRO: arquivo .cdf n�o pode ser aberto!" << std::endl;
		return 1;
	}
	return 0;
}

bool readCDFX(std::string cdfFile, sistema& sistema, barra& barra, ramo& ramo) {
	std::string line;
	std::fstream CDF;
	CDF.open(cdfFile.c_str(), std::ios::in);
	if (CDF.is_open())
	{
		sistema.nB = 0;
		sistema.nPQ = 0;
		sistema.nPV = 0;
		sistema.barraVO = 0;
		std::getline(CDF, line);

		sistema.baseMVA = atoi(line.substr(31, 6).c_str());

		int i = 0; 

		std::getline(CDF, line); 
		while (std::getline(CDF, line)) 
		{
			if (line.find("-999") == std::string::npos) {
				i++; 

				barra.id[IDX1F(i)] = atoi(line.substr(0, 5).c_str());

				barra.Vbase[IDX1F(i)] = atof(line.substr(77, 7).c_str());  

				switch (atoi(line.substr(25, 2).c_str())) { 
				case 0: 
					sistema.barrasPQ[sistema.nPQ] = i;
					sistema.nPQ++;
					break;
				case 1: 
					sistema.barrasPQ[sistema.nPQ] = i;
					sistema.nPQ++;
					break;
				case 2: 
					sistema.barrasPV[sistema.nPV] = i;
					sistema.nPV++;

					barra.V[IDX1F(i)] = atof(line.substr(28, 6).c_str());
					sistema.VfixadaPV[IDX1F(sistema.nPV)] = barra.V[IDX1F(i)]; 

					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(91, 7).c_str()) / sistema.baseMVA;  
					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(99, 7).c_str()) / sistema.baseMVA;  

					break;
				case 3:
					if (sistema.barraVO == 0) {
						sistema.barraVO = i;

						barra.V[IDX1F(i)] = atof(line.substr(28, 6).c_str());

						barra.theta[IDX1F(i)] = atof(line.substr(34, 7).c_str()) * 3.14159265358979323846 / 180.;  
					}
					else {
						std::cout << "ERRO: mais de uma barrra swing definida!" << std::endl;
						return 1;
					}
					break;
				default:
					std::cout << "ERRO: valor de tipo de barra inválido!" << std::endl;
					return 1;
				}
				
				barra.Pload[IDX1F(i)] = atof(line.substr(41, 9).c_str()) / sistema.baseMVA;  
				barra.Qload[IDX1F(i)] = atof(line.substr(50, 10).c_str()) / sistema.baseMVA; 
				barra.Pg[IDX1F(i)] = atof(line.substr(60, 8).c_str()) / sistema.baseMVA;  
				barra.Qg[IDX1F(i)] = atof(line.substr(68, 8).c_str()) / sistema.baseMVA;  
				barra.Vbase[IDX1F(i)] = atof(line.substr(77, 7).c_str());  
				barra.bsh[IDX1F(i)] = atof(line.substr(115, 8).c_str()); 
				barra.gsh[IDX1F(i)] = atof(line.substr(107, 8).c_str()); 

				barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  
				barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  

			}
			else
				break; 
		}

		sistema.nB = sistema.nPQ + sistema.nPV + 1;

		std::getline(CDF, line); 
		sistema.nL = 0;
		int nRamosDuplicatas = 0;
		while (std::getline(CDF, line)) 
		{
			if (line.find("-999") == std::string::npos) {
				
				int auxde = id2i(atoi(line.substr(0, 5).c_str()), sistema, barra); 
				int auxpara = id2i(atoi(line.substr(6, 5).c_str()), sistema, barra); 
				bool flgRamoNovo = 1;
				int i;

				for (i = 0; i < sistema.nL; i++) {
					if ((ramo.de[i] == auxde) && (ramo.para[i] == auxpara)) {
						
						flgRamoNovo = 0;
						nRamosDuplicatas++;
						break;
					}
				}
				if (flgRamoNovo) { 
					sistema.nL++;
					ramo.de[IDX1F(sistema.nL)] = auxde;
					ramo.para[IDX1F(sistema.nL)] = auxpara;

					ramo.z[IDX1F(sistema.nL)] = _mkComplex(atof(line.substr(21, 10).c_str()),  
						atof(line.substr(31, 11).c_str())); 

					float_type auxA = atof(line.substr(78, 6).c_str()); 
					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(85, 7).c_str()) * 3.14159265358979323846 / 180.; 
					if (auxA != 0.) {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); 
					}
					else {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
					}

					ramo.bsh[IDX1F(sistema.nL)] = atof(line.substr(42, 10).c_str()); 
				}
				else { 
					if (!global::laconic_mode) {
						printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", ramo.de[IDX1F(sistema.nL)], ramo.para[IDX1F(sistema.nL)]);
					}

					float_type auxA = atof(line.substr(78, 6).c_str()); 
					float_type auxPhi = atof(line.substr(85, 7).c_str()) * 3.14159265358979323846 / 180.;
					complex_type aux;
					if ((auxA != 0.)) {
						aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); 
					}
					else {
						aux = _mkComplex(cos(auxPhi), sin(auxPhi));
					}
					
					if ((aux.x != ramo.tap[i].x) || (aux.y != ramo.tap[i].y))
						printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);

					aux = _mkComplex(.0, .0);

					aux = _mkComplex(atof(line.substr(21, 10).c_str()),  
						atof(line.substr(31, 11).c_str())); 

					ramo.z[i] = _cuDiv(_cuMul(aux, ramo.z[i]), _cuAdd(ramo.z[i], aux)); 

					float_type aux2 = atof(line.substr(42, 10).c_str());   
					if (aux2 != 0 || ramo.bsh[i] != 0)
						ramo.bsh[i] = 1 / (aux2 + ramo.bsh[i]); 
					else
						ramo.bsh[i] = 0;
				}
			}
			else
				break; 
		}
		CDF.close();
	}
	else {
		std::cout << "ERRO: arquivo .cdf não pode ser aberto!" << std::endl;
		return 1;
	}

	return 0;
}

void lerArquivoEAlocarMemoria(sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo) {
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
		global::isMatpower = false;
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
		global::isMatpower = false;
	}
	else if (ext == "m") {
		
		matPowerDataType mpData = lerMatPowerEAlocarMemoria(global::arq_entrada, sistema, barra, ramo, iterativo);

		InitCsrPhi(sistema, ramo);
		global::isMatpower = true;
	}
	else {
		printf("nda == %s", ext.c_str());
	}
}
