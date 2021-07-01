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

		std::getline(CDF, line); std::getline(CDF, line); // pula linha de cabeçalho
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				sistema.nB++;	// conta barras	
				switch (atoi(line.substr(25, 2).c_str())) { // [24, 25]
				case 0: // tipo PQ [Unregulated (load, PQ)]
					sistema.nPQ++;
					break;
				case 1: // tipo PQ [Hold MVAR generation within voltage limits, (PQ)]
					sistema.nPQ++;
					break;
				case 2: // tipo PV [Hold voltage within VAR limits (gen, PV)]
					sistema.nPV++;
					break;
				case 3:// tipo VO [Hold voltage and angle (swing, V-Theta)]
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
				break; // se -999 a se��o Bus Data acabou!
		}
		std::getline(CDF, line); // pula duas linhas
		sistema.nL = 0;
		std::vector<unsigned int> de, para;

		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				unsigned int auxde = atoi(line.substr(0, 5).c_str()); // ramo da linha i // [0, 3]
				unsigned int auxpara = atoi(line.substr(6, 5).c_str()); // at� a j // [4, 8]
				bool flgRamoNovo = 1;

				for (unsigned int i = 0; i < sistema.nL; i++) {
					if ((de[i] == auxde) && (para[i] == auxpara)) {
						// ramo atual é o mesmo que o i-ésimo ramo
						// soma elementos em paralelo

						flgRamoNovo = 0;
						break;
					}
				}
				if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
					de.push_back(auxde);
					para.push_back(auxpara);
					sistema.nL++; // conta ramos
				}
			}
			else
				break; // se -999 a se��o ramo Data acabou!
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
	//std::vector <std::vector <branch>> vecBranch;

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

		//std::cout << "linha: " << atoi(line.substr(31, 6).c_str()) << std::endl;
		sistema.baseMVA = atoi(line.substr(31, 6).c_str());

		unsigned int i = 0; // i é o número da barra que está sendo lida

		std::getline(CDF, line); // pula duas linhas
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				i++; // nova barra encontrada

				barra.id[IDX1F(i)] = atoi(line.substr(0, 5).c_str());

				barra.Vbase[IDX1F(i)] = atof(line.substr(77, 7).c_str());  // Base KV // [76, 82]

				switch (atoi(line.substr(25, 2).c_str())) { // [24, 25]
				case 0: // tipo PQ [Unregulated (load, PQ)]
					sistema.barrasPQ[sistema.nPQ] = i;
					sistema.nPQ++;
					break;
				case 1: // tipo PQ [Hold MVAR generation within voltage limits, (PQ)]
					sistema.barrasPQ[sistema.nPQ] = i;
					sistema.nPQ++;
					break;
				case 2: // tipo PV [Hold voltage within VAR limits (gen, PV)]
					sistema.barrasPV[sistema.nPV] = i;
					sistema.nPV++;

					barra.V[IDX1F(i)] = atof(line.substr(28, 6).c_str());
					sistema.VfixadaPV[IDX1F(sistema.nPV)] = barra.V[IDX1F(i)]; // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat

					//Limites
					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(91, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(99, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]

					break;
				case 3:// tipo VO [Hold voltage and angle (swing, V-Theta)]
					if (sistema.barraVO == 0) {
						sistema.barraVO = i;

						barra.V[IDX1F(i)] = atof(line.substr(28, 6).c_str());

						barra.theta[IDX1F(i)] = atof(line.substr(34, 7).c_str()) * 3.14159265358979323846 / 180.;  // Final angle, degrees
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
				barra.Pload[IDX1F(i)] = atof(line.substr(41, 9).c_str()) / sistema.baseMVA;  // Load MW // [40, 48]
				barra.Qload[IDX1F(i)] = atof(line.substr(50, 10).c_str()) / sistema.baseMVA; // Load MVAR // [49, 58]
				barra.Pg[IDX1F(i)] = atof(line.substr(60, 8).c_str()) / sistema.baseMVA;  // Generation MW // [59, 66]
				barra.Qg[IDX1F(i)] = atof(line.substr(68, 8).c_str()) / sistema.baseMVA;  // Generation MVAR // [67, 74]
				barra.Vbase[IDX1F(i)] = atof(line.substr(77, 7).c_str());  // Base KV // [76, 82]
				barra.bsh[IDX1F(i)] = atof(line.substr(115, 8).c_str()); // 
				barra.gsh[IDX1F(i)] = atof(line.substr(107, 8).c_str()); // Generation MVAR // [67, 74]

				barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
				barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos

			}
			else
				break; // se -999 a se��o Bus Data acabou!
		}

		// número de barras
		sistema.nB = sistema.nPQ + sistema.nPV + 1;

		std::getline(CDF, line); // pula duas linhas
		sistema.nL = 0;
		unsigned int nRamosDuplicatas = 0;
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				// Parse branch	
				unsigned int auxde = id2i(atoi(line.substr(0, 5).c_str()), sistema, barra); // ramo da linha i 
				unsigned int auxpara = id2i(atoi(line.substr(6, 5).c_str()), sistema, barra); // at� a j 
				bool flgRamoNovo = 1;
				unsigned int i;

				for (i = 0; i < sistema.nL; i++) {
					if ((ramo.de[i] == auxde) && (ramo.para[i] == auxpara)) {
						// ramo atual é o mesmo que o i-ésimo ramo
						// soma elementos em paralelo

						flgRamoNovo = 0;
						nRamosDuplicatas++;
						break;
					}
				}
				if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
					sistema.nL++;
					ramo.de[IDX1F(sistema.nL)] = auxde;
					ramo.para[IDX1F(sistema.nL)] = auxpara;

					ramo.z[IDX1F(sistema.nL)] = _mkComplex(atof(line.substr(21, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
						atof(line.substr(31, 11).c_str())); // ramo reactance  X, per unit [29, 39]

					float_type auxA = atof(line.substr(78, 6).c_str()); // Transformer final turns ratio (F)
					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(85, 7).c_str()) * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)
					if (auxA != 0.) {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
					}
					else {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
					}

					ramo.bsh[IDX1F(sistema.nL)] = atof(line.substr(42, 10).c_str()); // Line charging B, per unit (F) * (total line charging, +B)
				}
				else { // ramo i é duplicata do nL lido
					if (!global::laconic_mode) {
						printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", ramo.de[IDX1F(sistema.nL)], ramo.para[IDX1F(sistema.nL)]);
					}

					// conferência do tap
					float_type auxA = atof(line.substr(78, 6).c_str()); // tap
					float_type auxPhi = atof(line.substr(85, 7).c_str()) * 3.14159265358979323846 / 180.;
					complex_type aux;
					if ((auxA != 0.)) {
						aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); // a * cis phi
					}
					else {
						aux = _mkComplex(cos(auxPhi), sin(auxPhi));
					}
					// se valor lido é diferente do já armazenado: atenção!
					if ((aux.x != ramo.tap[i].x) || (aux.y != ramo.tap[i].y))
						printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);

					// impedância série
					aux = _mkComplex(.0, .0);

					aux = _mkComplex(atof(line.substr(21, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
						atof(line.substr(31, 11).c_str())); // ramo reactance  X, per unit [29, 39]

					ramo.z[i] = _cuDiv(_cuMul(aux, ramo.z[i]), _cuAdd(ramo.z[i], aux)); // A//B = (A*B)/(A+B)

					// bsh
					float_type aux2 = atof(line.substr(42, 10).c_str());   // Line charging B, per unit (F) * (total line charging, +B)
					if (aux2 != 0 || ramo.bsh[i] != 0)
						ramo.bsh[i] = 1 / (aux2 + ramo.bsh[i]); //errado:= ramo.bsh[i]*aux2/(ramo.bsh[i]+aux2); // A//B = (A*B)/(A+B)
					else
						ramo.bsh[i] = 0;
				}
			}
			else
				break; // se -999 a se��o ramo Data acabou!
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
	else {
		printf("nda == %s", ext.c_str());
	}
}
