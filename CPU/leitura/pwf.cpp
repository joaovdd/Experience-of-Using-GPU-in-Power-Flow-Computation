
#include "pwf.h"

bool readPWF(sistemaType& sistema, barraType& barra, ramoType& ramo, DCTEtype& DCTE, DBARtype& DBAR, DLINtype& DLIN, DGBTtype& DGBT, DCSCtype& DCSC) {
	sistema.nB = 0;
	sistema.nPQ = 0;
	sistema.nPV = 0;
	sistema.barraVO = 0;

	sistema.baseMVA = DCTE.BASE;

	for (size_t i = 1; i <= DBAR.numero.size(); i++) {
		barra.id[IDX1F(i)] = DBAR.numero[IDX1F(i)];

		switch (DBAR.tipo[IDX1F(i)]) {
		case '3':
			std::cout << "AVISO: Barra com limites de tensão inserida" << std::endl; 
		case '0': 
			sistema.barrasPQ[sistema.nPQ] = i;
			sistema.nPQ++;
			break;
		case '1': 
			sistema.barrasPV[sistema.nPV] = i;
			sistema.nPV++;

			barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];

			sistema.limQsup[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMaxima[IDX1F(i)] / sistema.baseMVA;  
			sistema.limQinf[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMinima[IDX1F(i)] / sistema.baseMVA;  

			break;
		case '2':
			if (sistema.barraVO == 0) {
				sistema.barraVO = i;

				barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];

				barra.theta[IDX1F(i)] = DBAR.angulo[IDX1F(i)] * 3.14159265358979323846 / 180.;  
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
		
		barra.Pload[IDX1F(i)] = DBAR.cargaAtiva[IDX1F(i)] / sistema.baseMVA;  
		barra.Qload[IDX1F(i)] = DBAR.cargaReativa[IDX1F(i)] / sistema.baseMVA; 
		barra.Pg[IDX1F(i)] = DBAR.geracaoAtiva[IDX1F(i)] / sistema.baseMVA;  
		barra.Qg[IDX1F(i)] = DBAR.geracaoReativa[IDX1F(i)] / sistema.baseMVA;  
		barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[IDX1F(i)]];  
		
		barra.bsh[IDX1F(i)] = 0; 
		barra.gsh[IDX1F(i)] = 0; 

		barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  
		barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  
	}

	sistema.nB = sistema.nPQ + sistema.nPV + 1;

	sistema.nL = 0;
	int nRamosDuplicatas = 0;
	for (size_t i = 0; i < DLIN.daBarra.size(); i++) 
	{
		int auxde = id2i(DLIN.daBarra[i], sistema, barra); 
		int auxpara = id2i(DLIN.paraBarra[i], sistema, barra); 
		bool flgRamoNovo = 1;
		int j; 

		for (j = 0; j < sistema.nL; j++) {
			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
				
				flgRamoNovo = 0;
				nRamosDuplicatas++;
				break;
			}
		}
		if (flgRamoNovo) { 
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; 
			ramo.para[IDX1F(sistema.nL)] = auxpara; 

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(DLIN.resistencia[i] / 100,  
				DLIN.reatancia[i] / 100); 

			float_type auxA = DLIN.tap[i]; 
			ramo.phi[IDX1F(sistema.nL)] = DLIN.defasagem[i] * 3.14159265358979323846 / 180.; 

			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); 

			ramo.bsh[IDX1F(sistema.nL)] = DLIN.susceptancia[i] / sistema.baseMVA; 
		}
		else { 
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
			}

			float_type auxA = DLIN.tap[i]; 
			float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
			complex_type aux;
			aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); 

			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %zu)\n", i);
				}
			}
			
			aux = _mkComplex(DLIN.resistencia[i] / 100,  
				DLIN.reatancia[i] / 100); 

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); 

			float_type aux2 = DLIN.susceptancia[i];   
			if (aux2 != 0 || ramo.bsh[j] != 0)
				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
			else
				ramo.bsh[j] = 0;
		}
	}

	for (size_t i = 0; i < DCSC.daBarra.size(); i++) 
	{
		int auxde = id2i(DCSC.daBarra[i], sistema, barra); 
		int auxpara = id2i(DCSC.paraBarra[i], sistema, barra); 
		bool flgRamoNovo = 1;
		int j; 

		for (j = 0; j < sistema.nL; j++) {
			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
				
				flgRamoNovo = 0;
				nRamosDuplicatas++;
				break;
			}
		}
		if (flgRamoNovo) { 
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; 
			ramo.para[IDX1F(sistema.nL)] = auxpara; 

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(0.0,  
				                                             DCSC.valorEspecificado[i] / 100); 

			ramo.phi[IDX1F(sistema.nL)] = 0.0; 

			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(1.0, 0.0); 

			ramo.bsh[IDX1F(sistema.nL)] = 0.0; 
		}
		else { 
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
			}

			complex_type aux;
			aux = _mkComplex(1.0, 0.0); 

			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %zu)\n", i);
				}
			}
			
			aux = _mkComplex(0.0,  
				                       DCSC.valorEspecificado[i] / 100); 

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); 

			float_type aux2 = 0.0;   
			if (aux2 != 0 || ramo.bsh[j] != 0)
				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
			else
				ramo.bsh[j] = 0;
		}
	}

	return 0;
}

std::string ifBlanckSetDefault(std::string str, std::string defaultValue) {
	auto aux = std::string(str);
	aux.erase(remove_if(aux.begin(), aux.end(), ::isspace), aux.end());
	if (aux == "" ) {
		return defaultValue;
	}
	return (str);
}

char ifBlanckSetDefault(char str, char defaultValue) {
	if (str == '\0' || str == ' ') {
		return defaultValue;
	}
	return (str);
}

float parseFloatPtoDecImp(std::string str, int casaUnitariaDoPontoDecimalImplicito) {
	auto aux = str.find('.');
	if (aux != std::string::npos) {
		return atof(str.c_str());
	}
	else {
		return atoi(str.c_str()) / (pow(10, casaUnitariaDoPontoDecimalImplicito));
	}
}

void print(std::string objeto, std::string& str_pwf) {
	auto pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); 

	if (objeto == "TITU") {
		
		while (true) {
			if (str_pwf[pnt] == '(') {
				
				pnt = str_pwf.find("\n", pnt) + 1;
				continue;
			}
			break;
		}

		std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);
	}
	else {
		
		while (true) {
			if (str_pwf[pnt] == '(') {
				
				pnt = str_pwf.find("\n", pnt) + 1;
				continue;
			}

			std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

			pnt = str_pwf.find("\n", pnt) + 1; 

			if (pnt == str_pwf.find("99999", pnt)) {
				
				break;
			}
		}
	}

}

void parseDBARline(DBARtype& DBAR, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio , fim - inicio + 1);

	DBAR.numero.push_back(atoi(linha.substr(0, 5).c_str()));
	DBAR.operacao.push_back(ifBlanckSetDefault(linha[5], 'A'));
	DBAR.estado.push_back(ifBlanckSetDefault(linha[6], 'L'));
	DBAR.tipo.push_back(ifBlanckSetDefault(linha[7], '0'));

	DBAR.grupoDeBaseDeTensao.push_back(ifBlanckSetDefault(linha.substr(8, 2), "0")); 
	DBAR.nome.push_back(linha.substr(10, 11));

	DBAR.grupoDeLimiteDeTensao.push_back(ifBlanckSetDefault(linha.substr(2, 2), "0")); 

	DBAR.tensao.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(24, 4), "1000"), 3)); 
	DBAR.angulo.push_back(atof(ifBlanckSetDefault(linha.substr(28, 4), "0.0").c_str())); 
	DBAR.geracaoAtiva.push_back(atof(ifBlanckSetDefault(linha.substr(32, 5), "0.0").c_str())); 
	DBAR.geracaoReativa.push_back(atof(ifBlanckSetDefault(linha.substr(37, 5), "0.0").c_str())); 
	DBAR.geracaoReativaMinima.push_back(atof(ifBlanckSetDefault(linha.substr(42, 5), "-9999.0").c_str())); 
	DBAR.geracaoReativaMaxima.push_back(atof(ifBlanckSetDefault(linha.substr(47, 5), "99999.0").c_str()));
	DBAR.barraControlada.push_back(atoi(linha.substr(52, 6).c_str())); 
	DBAR.cargaAtiva.push_back(atof(ifBlanckSetDefault(linha.substr(58, 5), "0.0").c_str())); 
	DBAR.cargaReativa.push_back(atof(ifBlanckSetDefault(linha.substr(63, 5), "0.0").c_str())); 
	DBAR.capacitorReator.push_back(atof(ifBlanckSetDefault(linha.substr(68, 5), "0.0").c_str())); 
	DBAR.area.push_back(atoi(ifBlanckSetDefault(linha.substr(73, 3), "1").c_str()));

	DBAR.tensaoParaDefinicaoDeCarga.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(76, 4), "1000"), 3)); 
	DBAR.modoDeVisualizacao.push_back(ifBlanckSetDefault(linha[80], '0'));

}

void parseDLINline(DLINtype& DLIN, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio , fim - inicio + 1);

	DLIN.daBarra.push_back(atoi(linha.substr(0, 5).c_str()));
	DLIN.aberturaDaBarra.push_back(ifBlanckSetDefault(linha[5], 'L'));
	DLIN.operacao.push_back(ifBlanckSetDefault(linha[7], 'A'));
	DLIN.aberturaParaABarra.push_back(ifBlanckSetDefault(linha[9], 'L'));
	DLIN.paraBarra.push_back(atoi(linha.substr(10, 5).c_str()));
	DLIN.circuito.push_back(atoi(linha.substr(15, 2).c_str())); 
	DLIN.estado.push_back(ifBlanckSetDefault(linha[17], 'L'));
	DLIN.proprietario.push_back(ifBlanckSetDefault(linha[18], 'F'));

	DLIN.resistencia.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(20, 6), "0.0"), 2)); 

	DLIN.reatancia.push_back(parseFloatPtoDecImp(linha.substr(26, 6), 2)); 

	DLIN.susceptancia.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(32, 6), "0.0"), 3)); 

	DLIN.tap.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(38, 5), "1.0"), 3)); 
	DLIN.tapMinimo.push_back(atof(linha.substr(43, 5).c_str())); 
	DLIN.tapMaximo.push_back(atof(linha.substr(48, 5).c_str())); 

	DLIN.defasagem.push_back(parseFloatPtoDecImp(linha.substr(53, 5), 2)); 
	DLIN.barraControlada.push_back(atoi(linha.substr(58, 6).c_str())); 
	DLIN.capacidadeNormal.push_back(atof(linha.substr(64, 4).c_str())); 
	DLIN.capacidadeEmEmergencia.push_back(atof(linha.substr(68, 4).c_str())); 
	DLIN.numeroDeTaps.push_back(atoi(linha.substr(72, 2).c_str())); 
	DLIN.capacidadeDeEquipamento.push_back(atof(linha.substr(74, 4).c_str())); 

}

void parseDCTEline(DCTEtype& DCTE, std::string& str_pwf, unsigned long long& inicio) { 
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio , fim - inicio + 1);

}

void parseDCSCline(DCSCtype& DCSC, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio , fim - inicio + 1);

	DCSC.daBarra.push_back(atoi(linha.substr(0, 5).c_str()));
	DCSC.operacao.push_back(ifBlanckSetDefault(linha[6], 'A'));
	DCSC.paraBarra.push_back(atoi(linha.substr(9, 5).c_str()));
	DCSC.circuito.push_back(atoi(linha.substr(14, 2).c_str())); 
	DCSC.estado.push_back(ifBlanckSetDefault(linha[16], 'L'));

	DCSC.bypass.push_back(ifBlanckSetDefault(linha[18], 'D'));

	DCSC.modoDeControle.push_back(ifBlanckSetDefault(linha[43], 'X'));
	DCSC.valorEspecificado.push_back(atof(ifBlanckSetDefault(linha.substr(45, 6), "0.0").c_str()));

}

void parseDGBTline(DGBTtype& DGBT, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio, fim - inicio + 1);

	auto aux = linha.substr(0, 2); 
	aux.erase(remove_if(aux.begin(), aux.end(), ::isspace), aux.end());

	DGBT.mapa[aux] = atof(linha.substr(3, 7).c_str()); 
}

void parse(DBARtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DBAR";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); 

	while (true) {
		if (str_pwf[pnt] == '(') {
			
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		parseDBARline(estrutura, str_pwf, pnt);

		pnt = str_pwf.find("\n", pnt) + 1; 

		if (pnt == str_pwf.find("99999", pnt)) {
			
			break;
		}
	}
}

void parse(DLINtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DLIN";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); 

	while (true) {
		if (str_pwf[pnt] == '(') {
			
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		parseDLINline(estrutura, str_pwf, pnt);

		pnt = str_pwf.find("\n", pnt) + 1; 

		if (pnt == str_pwf.find("99999", pnt)) {
			
			break;
		}
	}
}

void parse(DCTEtype& estrutura, std::string& str_pwf) {
	auto inicio = str_pwf.find("DCTE\n") + std::string("DCTE\n").size();

	auto fim = str_pwf.find("99999", inicio);

	std::string auxDCTEstr = str_pwf.substr(inicio, fim - inicio);

	auto aux = auxDCTEstr.find("BASE ");
	if (aux != std::string::npos) {
		aux += 5;
		estrutura.BASE = atof(auxDCTEstr.substr(aux, 6).c_str());
	}
}

void parse(DGBTtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DGBT";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); 

	while (true) {
		if (str_pwf[pnt] == '(') {
			
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		parseDGBTline(estrutura, str_pwf, pnt);

		pnt = str_pwf.find("\n", pnt) + 1; 

		if (pnt == str_pwf.find("99999", pnt)) {
			
			break;
		}
	}
}

void parse(DCSCtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DCSC";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); 

	while (true) {
		if (str_pwf[pnt] == '(') {
			
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		parseDCSCline(estrutura, str_pwf, pnt);

		pnt = str_pwf.find("\n", pnt) + 1; 

		if (pnt == str_pwf.find("99999", pnt)) {
			
			break;
		}
	}
}

void getPWF(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo) {
	std::fstream PWF;
	PWF.open(pwfFile.c_str(), std::ios::in);

	std::string str_pwf((std::istreambuf_iterator<char>(PWF)), (std::istreambuf_iterator<char>()));

	DBARtype DBAR;

	parse(DBAR, str_pwf);

	DLINtype DLIN;

	parse(DLIN, str_pwf);

	DCTEtype DCTE;

	parse(DCTE, str_pwf);

	DGBTtype DGBT;

	parse(DGBT, str_pwf);

	DCSCtype DCSC;

	parse(DCSC, str_pwf);

	readPWF(sistema, barra, ramo, DCTE, DBAR, DLIN, DGBT, DCSC);
}

void lerPWFEAlocarMemoria(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	std::fstream PWF;
	PWF.open(pwfFile.c_str(), std::ios::in);

	std::string str_pwf((std::istreambuf_iterator<char>(PWF)), (std::istreambuf_iterator<char>()));

	PWF.close();

	DBARtype DBAR;

	parse(DBAR, str_pwf);

	DLINtype DLIN;

	parse(DLIN, str_pwf);

	DCTEtype DCTE;

	parse(DCTE, str_pwf);

	DGBTtype DGBT;

	parse(DGBT, str_pwf);

	DCSCtype DCSC;

	parse(DCSC, str_pwf);

	sistema.nB  = DBAR.numero.size();
	sistema.nL  = DLIN.daBarra.size() + DCSC.daBarra.size(); 
	sistema.nPQ = std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '0') + std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '3');
	sistema.nPV = std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '1');

	initBranch(sistema, ramo);
	initSistema(sistema);
	initBus(sistema, barra);
	initIter(sistema, iterativo);

	if (readPWF(sistema, barra, ramo, DCTE, DBAR, DLIN, DGBT, DCSC)) { 
		printf("Deu ruim...");
	}
}