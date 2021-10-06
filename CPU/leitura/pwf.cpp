// Baseado no Manual do Programa de Análise de Redes (Anarede) V10.02.03
//
// (C) 2020 João Victor Daher Daibes


// Manual: file:///C:/Cepel/Anarede/V100203/Manuais/Manual-Anarede.pdf
// DBAR -> p82 (2.13.4)
// DLIN -> p153(2.54.4)

// DAGR -> p73 (2.9.4)				Leitura dos dados de agregadores genéricos.
// DOPC -> p167(2.62.4)				Leitura e modificação dos dados das Opções de Controle de Execução padrão.
// DCTE -> p111(2.29.1)				Leitura e modificação dos dados de constantes utilizadas no programa.
// DCSC -> p109(2.28.4)				Leitura dos dados de CSC (Compensador Série Controlável).
// DBSH // FBAN -> p87  (2.15.4)	Leitura  dos  dados  de  bancos  de  capacitores  e/ou  reatores  individualizados  conectados  a  barras  CA  ou  linhas  de  transmissão. 
// DSHL -> p176(2.69.4)				Leitura dos dados de dispositivos shunt de circuito CA.
// DCER -> p100(2.22.4)				Leitura dos dados de compensador estático de reativos. 
// DCTR -> p123(2.31.4)				Leitura  dos  dados  complementares  de  transformadores.
// DGBT -> p135(2.40.4)				Leitura dos dados de grupos de base de tensão de barras CA.
// DGLT -> p145
// DARE -> p77
// DTPF -> p180
// DELO -> p125
// DCBA -> p95
// DCLI -> p103
// DCNV ->

#include "pwf.h"

bool readPWF(sistemaType& sistema, barraType& barra, ramoType& ramo, DCTEtype& DCTE, DBARtype& DBAR, DLINtype& DLIN, DGBTtype& DGBT, DCSCtype& DCSC) {

	sistema.nB = 0;
	sistema.nPQ = 0;
	sistema.nPV = 0;
	sistema.barraVO = 0;

	sistema.baseMVA = DCTE.BASE;

	for (size_t i = 1; i <= DBAR.numero.size(); i++) {

		barra.id[IDX1F(i)] = DBAR.numero[IDX1F(i)];

		//barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[i]];  // Base kV // repetida...

		switch (DBAR.tipo[IDX1F(i)]) {
		case '3':
			std::cout << "AVISO: Barra com limites de tensão inserida" << std::endl; //   barra  de  carga  com  limite  de  tensão  (PQ  -  Injeções  de  potências  ativa  e  reativa  fixas  enquanto  a  magnitude  de  tensão  permanecer  entre os valores limites)
		case '0': // tipo PQ [Unregulated (load, PQ)]
			sistema.barrasPQ[sistema.nPQ] = i;
			sistema.nPQ++;
			break;
		case '1': // tipo PV [Hold voltage within VAR limits (gen, PV)]
			sistema.barrasPV[sistema.nPV] = i;
			sistema.nPV++;

			barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];

			//Limites
			sistema.limQsup[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMaxima[IDX1F(i)] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
			sistema.limQinf[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMinima[IDX1F(i)] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]

			break;
		case '2':// tipo VO [Hold voltage and angle (swing, V-Theta)]
			if (sistema.barraVO == 0) {
				sistema.barraVO = i;

				barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];

				barra.theta[IDX1F(i)] = DBAR.angulo[IDX1F(i)] * 3.14159265358979323846 / 180.;  // Final angle, degrees
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
		//std::string auxx = line.substr(40, 9).c_str();
		barra.Pload[IDX1F(i)] = DBAR.cargaAtiva[IDX1F(i)] / sistema.baseMVA;  // Load MW 
		barra.Qload[IDX1F(i)] = DBAR.cargaReativa[IDX1F(i)] / sistema.baseMVA; // Load MVAR
		barra.Pg[IDX1F(i)] = DBAR.geracaoAtiva[IDX1F(i)] / sistema.baseMVA;  // Generation MW // [59, 66]
		barra.Qg[IDX1F(i)] = DBAR.geracaoReativa[IDX1F(i)] / sistema.baseMVA;  // Generation MVAR // [67, 74]
		barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[IDX1F(i)]];  // Base KV
		// onde estao esses valores sh??
		barra.bsh[IDX1F(i)] = 0; // 
		barra.gsh[IDX1F(i)] = 0; // Generation MVAR // [67, 74]

		barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
		barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos
	}

	// número de barras
	sistema.nB = sistema.nPQ + sistema.nPV + 1;

	sistema.nL = 0;
	int nRamosDuplicatas = 0;
	for (size_t i = 0; i < DLIN.daBarra.size(); i++) // i itera ramos
	{
		int auxde = id2i(DLIN.daBarra[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
		int auxpara = id2i(DLIN.paraBarra[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
		bool flgRamoNovo = 1;
		int j; // j aponta número do ramo repetido

		for (j = 0; j < sistema.nL; j++) {
			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
				// ramo atual é o mesmo que o i-ésimo ramo
				// soma elementos em paralelo

				flgRamoNovo = 0;
				nRamosDuplicatas++;
				break;
			}
		}
		if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
			ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(DLIN.resistencia[i] / 100,  // ramo resistance R, per unit [19, 28] ,
				DLIN.reatancia[i] / 100); // ramo reactance  X, per unit [29, 39]

			float_type auxA = DLIN.tap[i]; // Transformer final turns ratio (F) // em p.u.?
			ramo.phi[IDX1F(sistema.nL)] = DLIN.defasagem[i] * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)

			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi

			ramo.bsh[IDX1F(sistema.nL)] = DLIN.susceptancia[i] / sistema.baseMVA; // Line charging B, per unit (F) * (total line charging, +B)
		}
		else { // ramo i é duplicata do nL lido
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
			}

			// conferência do tap
			float_type auxA = DLIN.tap[i]; // tap
			float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
			complex_type aux;
			aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); // a * cis phi

			// se valor lido é diferente do já armazenado: atenção!
			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %zu)\n", i);
				}
			}
			// impedância série
			// aux = _mkComplex(.0, .0);

			aux = _mkComplex(DLIN.resistencia[i] / 100,  // ramo resistance R, per unit [19, 28] ,
				DLIN.reatancia[i] / 100); // ramo reactance  X, per unit [29, 39]

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)

			// bsh
			float_type aux2 = DLIN.susceptancia[i];   // Line charging B, per unit (F) * (total line charging, +B)
			if (aux2 != 0 || ramo.bsh[j] != 0)
				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
			else
				ramo.bsh[j] = 0;
		}
	}

	for (size_t i = 0; i < DCSC.daBarra.size(); i++) // i itera Compensadores Série Controláveis
	{
		int auxde = id2i(DCSC.daBarra[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
		int auxpara = id2i(DCSC.paraBarra[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
		bool flgRamoNovo = 1;
		int j; // j aponta número do ramo repetido

		for (j = 0; j < sistema.nL; j++) {
			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
				// ramo atual é o mesmo que o i-ésimo ramo
				// soma elementos em paralelo

				flgRamoNovo = 0;
				nRamosDuplicatas++;
				break;
			}
		}
		if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
			ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
				                                             DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]

			//float_type auxA = DLIN.tap[i]; // Transformer final turns ratio (F) // em p.u.?
			ramo.phi[IDX1F(sistema.nL)] = 0.0; // Transformer (phase shifter) final angle (F)

			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(1.0, 0.0); // a * cis phi

			ramo.bsh[IDX1F(sistema.nL)] = 0.0; // Line charging B, per unit (F) * (total line charging, +B)
		}
		else { // ramo i é duplicata do nL lido
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
			}

			// conferência do tap
			//float_type auxA = DLIN.tap[i]; // tap
			//float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
			complex_type aux;
			aux = _mkComplex(1.0, 0.0); // a * cis phi

			// se valor lido é diferente do já armazenado: atenção!
			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %zu)\n", i);
				}
			}
			// impedância série
			// aux = _mkComplex(.0, .0);

			aux = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
				                       DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)

			// bsh
			float_type aux2 = 0.0;   // Line charging B, per unit (F) * (total line charging, +B)
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
	if (aux == "" /*|| aux == " "*/) {
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

// lê valor real com possibilidade de ponto decimal implícito parsePtoDecImp
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
	auto pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); // aponta para o primeiro caractere da próxima linha

	if (objeto == "TITU") {
		// single line, not "99999" terminated
		while (true) {
			if (str_pwf[pnt] == '(') {
				// linha comentada
				pnt = str_pwf.find("\n", pnt) + 1;
				continue;
			}
			break;
		}

		// parse line

		std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);
	}
	else {
		// "99999" terminated
		while (true) {
			if (str_pwf[pnt] == '(') {
				// linha comentada
				pnt = str_pwf.find("\n", pnt) + 1;
				continue;
			}

			// parse line

			std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

			pnt = str_pwf.find("\n", pnt) + 1; // próxima linha

			if (pnt == str_pwf.find("99999", pnt)) {
				// atingiu o fim do bloco
				break;
			}
		}
	}

}

void parseDBARline(DBARtype& DBAR, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio/* + 1*/);
	auto linha = str_pwf.substr(inicio /*+ 1*/, fim - inicio + 1);

	DBAR.numero.push_back(atoi(linha.substr(0, 5).c_str()));
	DBAR.operacao.push_back(ifBlanckSetDefault(linha[5], 'A'));
	DBAR.estado.push_back(ifBlanckSetDefault(linha[6], 'L'));
	DBAR.tipo.push_back(ifBlanckSetDefault(linha[7], '0'));
	// DBAR.grupoDeBaseDeTensao.push_back(std::pair <char, char>(linha[8], linha[9])); //** DGBT
	DBAR.grupoDeBaseDeTensao.push_back(ifBlanckSetDefault(linha.substr(8, 2), "0")); //** DGBT
	DBAR.nome.push_back(linha.substr(10, 11));
	// DBAR.grupoDeLimiteDeTensao.push_back(std::pair <char, char>(linha[22], linha[23])); //**
	DBAR.grupoDeLimiteDeTensao.push_back(ifBlanckSetDefault(linha.substr(2, 2), "0")); //**
	// DBAR.tensao.push_back(atoi(ifBlanckSetDefault(linha.substr(24, 4), "1000").c_str()) / 1000.f); // valor em p.u. // ponto decimal implícito sempre considerado
	DBAR.tensao.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(24, 4), "1000"), 3)); // valor em p.u. // ponto decimal implícito sempre considerado
	DBAR.angulo.push_back(atof(ifBlanckSetDefault(linha.substr(28, 4), "0.0").c_str())); // valor em graus
	DBAR.geracaoAtiva.push_back(atof(ifBlanckSetDefault(linha.substr(32, 5), "0.0").c_str())); // valor em MW
	DBAR.geracaoReativa.push_back(atof(ifBlanckSetDefault(linha.substr(37, 5), "0.0").c_str())); // MVAr
	DBAR.geracaoReativaMinima.push_back(atof(ifBlanckSetDefault(linha.substr(42, 5), "-9999.0").c_str())); // ifBlanckSetDefault não considera caso em que apenas um está preenchido
	DBAR.geracaoReativaMaxima.push_back(atof(ifBlanckSetDefault(linha.substr(47, 5), "99999.0").c_str()));
	DBAR.barraControlada.push_back(atoi(linha.substr(52, 6).c_str())); // ifBlanckSetDefault não considerado
	DBAR.cargaAtiva.push_back(atof(ifBlanckSetDefault(linha.substr(58, 5), "0.0").c_str())); // MW
	DBAR.cargaReativa.push_back(atof(ifBlanckSetDefault(linha.substr(63, 5), "0.0").c_str())); // MVAr
	DBAR.capacitorReator.push_back(atof(ifBlanckSetDefault(linha.substr(68, 5), "0.0").c_str())); // MVAr
	DBAR.area.push_back(atoi(ifBlanckSetDefault(linha.substr(73, 3), "1").c_str()));
	// DBAR.tensaoParaDefinicaoDeCarga.push_back(atoi(ifBlanckSetDefault(linha.substr(76, 4), "1000").c_str())/1000.f); // p.u. // ponto decimal implícito sempre considerado
	DBAR.tensaoParaDefinicaoDeCarga.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(76, 4), "1000"), 3)); // p.u. // ponto decimal implícito sempre considerado
	DBAR.modoDeVisualizacao.push_back(ifBlanckSetDefault(linha[80], '0'));
	//DBAR.agregadores.push_back(atoi(linha.substr(81, 3).c_str())); //?
}

void parseDLINline(DLINtype& DLIN, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio/* + 1*/);
	auto linha = str_pwf.substr(inicio /*+ 1*/, fim - inicio + 1);

	//std::cout << linha;

	DLIN.daBarra.push_back(atoi(linha.substr(0, 5).c_str()));
	DLIN.aberturaDaBarra.push_back(ifBlanckSetDefault(linha[5], 'L'));
	DLIN.operacao.push_back(ifBlanckSetDefault(linha[7], 'A'));
	DLIN.aberturaParaABarra.push_back(ifBlanckSetDefault(linha[9], 'L'));
	DLIN.paraBarra.push_back(atoi(linha.substr(10, 5).c_str()));
	DLIN.circuito.push_back(atoi(linha.substr(15, 2).c_str())); // ifBlanckSetDefault não implementado
	DLIN.estado.push_back(ifBlanckSetDefault(linha[17], 'L'));
	DLIN.proprietario.push_back(ifBlanckSetDefault(linha[18], 'F'));
	// DLIN.resistencia.push_back(atof(ifBlanckSetDefault(linha.substr(20, 6), "0.0").c_str())); // % // ponto decimal implícito não implementado
	DLIN.resistencia.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(20, 6), "0.0"), 2)); // % 
	// DLIN.reatancia.push_back(atof(linha.substr(26, 6).c_str())); // % // ponto decimal implícito não implementado
	DLIN.reatancia.push_back(parseFloatPtoDecImp(linha.substr(26, 6), 2)); // %
	// DLIN.susceptancia.push_back(atof(ifBlanckSetDefault(linha.substr(32, 6), "0.0").c_str())); // % // ponto decimal implícito não implementado
	DLIN.susceptancia.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(32, 6), "0.0"), 3)); // %
	// DLIN.tap.push_back(atof(linha.substr(38, 5).c_str())); // p.u.// ponto decimal implícito não implementado
	DLIN.tap.push_back(parseFloatPtoDecImp(ifBlanckSetDefault(linha.substr(38, 5), "1.0"), 3)); // p.u.
	DLIN.tapMinimo.push_back(atof(linha.substr(43, 5).c_str())); // ponto decimal implícito não implementado
	DLIN.tapMaximo.push_back(atof(linha.substr(48, 5).c_str())); // ponto decimal implícito não implementado
	// DLIN.defasagem.push_back(atof(linha.substr(53, 5).c_str())); // graus
	DLIN.defasagem.push_back(parseFloatPtoDecImp(linha.substr(53, 5), 2)); // graus
	DLIN.barraControlada.push_back(atoi(linha.substr(58, 6).c_str())); // número // valor padrao nao implementado
	DLIN.capacidadeNormal.push_back(atof(linha.substr(64, 4).c_str())); // MVA // valor padrao nao implementado
	DLIN.capacidadeEmEmergencia.push_back(atof(linha.substr(68, 4).c_str())); // MVA // valor padrao nao implementado
	DLIN.numeroDeTaps.push_back(atoi(linha.substr(72, 2).c_str())); // número // valor padrao nao implementado
	DLIN.capacidadeDeEquipamento.push_back(atof(linha.substr(74, 4).c_str())); // Capacidade de carregamento do equipamento com menor capacidade de carregamento conectado ao circuito.  // valor padrao nao implementado
	//DLIN.agregadores.push_back(atof(linha.substr(, ).c_str()));
}

void parseDCTEline(DCTEtype& DCTE, std::string& str_pwf, unsigned long long& inicio) { // inacabado
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio /*+ 1*/, fim - inicio + 1);

	// TODO: 
}

// incompleto: checa apenas modo de controle sendo 'X' e lê Valor de Reatância Especificado.
void parseDCSCline(DCSCtype& DCSC, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio/* + 1*/);
	auto linha = str_pwf.substr(inicio /*+ 1*/, fim - inicio + 1);

	DCSC.daBarra.push_back(atoi(linha.substr(0, 5).c_str()));
	DCSC.operacao.push_back(ifBlanckSetDefault(linha[6], 'A'));
	DCSC.paraBarra.push_back(atoi(linha.substr(9, 5).c_str()));
	DCSC.circuito.push_back(atoi(linha.substr(14, 2).c_str())); // ifBlanckSetDefault não implementado'
	DCSC.estado.push_back(ifBlanckSetDefault(linha[16], 'L'));

	// ...

	DCSC.bypass.push_back(ifBlanckSetDefault(linha[18], 'D'));

	// ... 

	DCSC.modoDeControle.push_back(ifBlanckSetDefault(linha[43], 'X'));
	DCSC.valorEspecificado.push_back(atof(ifBlanckSetDefault(linha.substr(45, 6), "0.0").c_str()));

	// ...
}


void parseDGBTline(DGBTtype& DGBT, std::string& str_pwf, unsigned long long& inicio) {
	unsigned long long fim = str_pwf.find("\n", inicio);
	auto linha = str_pwf.substr(inicio, fim - inicio + 1);

	auto aux = linha.substr(0, 2); // retira espacos da chave do mapa
	aux.erase(remove_if(aux.begin(), aux.end(), ::isspace), aux.end());

	DGBT.mapa[aux] = atof(linha.substr(3, 7).c_str()); // kV
}

void parse(DBARtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DBAR";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); // aponta para o primeiro caractere da próxima linha

	while (true) {
		if (str_pwf[pnt] == '(') {
			// linha comentada
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		// parse line
		parseDBARline(estrutura, str_pwf, pnt);

		//std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

		pnt = str_pwf.find("\n", pnt) + 1; // próxima linha

		if (pnt == str_pwf.find("99999", pnt)) {
			// atingiu o fim do bloco
			break;
		}
	}
}

void parse(DLINtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DLIN";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); // aponta para o primeiro caractere da próxima linha

	while (true) {
		if (str_pwf[pnt] == '(') {
			// linha comentada
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		// parse line
		parseDLINline(estrutura, str_pwf, pnt);

		//std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

		pnt = str_pwf.find("\n", pnt) + 1; // próxima linha

		if (pnt == str_pwf.find("99999", pnt)) {
			// atingiu o fim do bloco
			break;
		}
	}
}

// incompleto
void parse(DCTEtype& estrutura, std::string& str_pwf) {
	// TODO: implementar estrutura DCTE com acesso por posições

	//std::array<const char*, 65> constants/*const char* constants[]*/ = 
	//{	"TEPA",	"TEPR",	"TLPR",	"TLVC",	"TLTC",	"TETP",	"TBPA",	"TSFR",	"TUDC",
	//	"TADC",	"BASE",	"DASE",	"ZMAX",	"ACIT",	"LPIT",	"LFLP",	"LFIT",	"DCIT",
	//	"VSIT",	"LCRT",	"LPRT",	"LFCV",	"TPST",	"QLST",	"EXST",	"TLPP",	"TSBZ",
	//	"TSBA",	"PGER",	"VDVN",	"VDVM",	"ASTP",	"VSTP",	"CSTP",	"VFLD",	"HIST",
	//	"ZMIN",	"PDIT",	"ICIT",	"FDIV",	"DMAX",	"ICMN",	"VART",	"TSTP",	"TSDC",
	//	"ASDC",	"ICMV",	"APAS",	"CPAR",	"VAVT",	"VAVF",	"VMVF",	"VPVT",	"VPVF",
	//	"VPMF",	"VSVF",	"VINF",	"VSUP",	"TLSI",	"NDIR",	"STIR",	"STTR",	"TRPT",
	//	"BFPO",	"LFPO"};

	auto inicio = str_pwf.find("DCTE\n") + std::string("DCTE\n").size();

	auto fim = str_pwf.find("99999", inicio);

	std::string auxDCTEstr = str_pwf.substr(inicio, fim - inicio);

	//for (const char* str : constants) {
	//	auto aux = auxDCTEstr.find(str);

	//	if (aux != std::string::npos) {

	//	} 

	//}

	// BASE
	auto aux = auxDCTEstr.find("BASE ");
	if (aux != std::string::npos) {
		aux += 5;
		estrutura.BASE = atof(auxDCTEstr.substr(aux, 6).c_str());
	}
}

void parse(DGBTtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DGBT";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); // aponta para o primeiro caractere da próxima linha

	while (true) {
		if (str_pwf[pnt] == '(') {
			// linha comentada
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		// parse line
		parseDGBTline(estrutura, str_pwf, pnt);

		//std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

		pnt = str_pwf.find("\n", pnt) + 1; // próxima linha

		if (pnt == str_pwf.find("99999", pnt)) {
			// atingiu o fim do bloco
			break;
		}
	}
}

void parse(DCSCtype& estrutura, std::string& str_pwf) {
	std::string objeto = "DCSC";
	unsigned long long pnt = str_pwf.find(objeto + std::string("\n")) + (objeto + std::string("\n")).size(); // aponta para o primeiro caractere da próxima linha

	while (true) {
		if (str_pwf[pnt] == '(') {
			// linha comentada
			pnt = str_pwf.find("\n", pnt) + 1;
			continue;
		}

		// parse line
		parseDCSCline(estrutura, str_pwf, pnt);

		//std::cout << str_pwf.substr(pnt, str_pwf.find("\n", pnt) - pnt + 1);

		pnt = str_pwf.find("\n", pnt) + 1; // próxima linha

		if (pnt == str_pwf.find("99999", pnt)) {
			// atingiu o fim do bloco
			break;
		}
	}
}

void getPWF(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo) {
	// entrada
	// std::string pwfFile = "c:/Users/PICHAU/Desktop/PD2029-MEDIA_NORTE SECO_2022.PWF"; sistema sistema; barra barra; ramo ramo;

	std::fstream PWF;
	PWF.open(pwfFile.c_str(), std::ios::in);

	std::string str_pwf((std::istreambuf_iterator<char>(PWF)), (std::istreambuf_iterator<char>()));

	//print("TITU", str_pwf);

	//print("DBAR", str_pwf);

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
	// entrada
	// std::string pwfFile = "c:/Users/PICHAU/Desktop/PD2029-MEDIA_NORTE SECO_2022.PWF"; sistema sistema; barra barra; ramo ramo;

	std::fstream PWF;
	PWF.open(pwfFile.c_str(), std::ios::in);

	std::string str_pwf((std::istreambuf_iterator<char>(PWF)), (std::istreambuf_iterator<char>()));

	PWF.close();

	//print("TITU", str_pwf);

	//print("DBAR", str_pwf);

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
	sistema.nL  = DLIN.daBarra.size() + DCSC.daBarra.size(); // n�mero de linhas + n�mero de compensadores s�rie
	sistema.nPQ = std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '0') + std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '3');
	sistema.nPV = std::count(DBAR.tipo.begin(), DBAR.tipo.end(), '1');

	initBranch(sistema, ramo);
	initSistema(sistema);
	initBus(sistema, barra);
	initIter(sistema, iterativo);

	
	if (readPWF(sistema, barra, ramo, DCTE, DBAR, DLIN, DGBT, DCSC)) { // guarda elementos lidos nas estruturas do Flumen
		printf("Deu ruim...");
	}
}