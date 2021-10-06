#include "m.h"

#include <iostream>
#include <string>
#include <regex>

#include <fstream>
#include <streambuf>

#define INF_VALUE "999999"
#define NAN_VALUE "0"

bool isToken(char chr, const std::string& tokens) {
	for (auto& i : tokens) {
		if (chr == i)
			return true;
	}
	return false;
}

int getNextToken(const std::string& str, const std::string& tokens, const int& startPos) {
	size_t pos = startPos;

	while (pos < str.size()) {
		if (isToken(str[pos], tokens))
			return pos;
		pos++;
	}

	return str.size();
}

std::vector<std::string> tokenize(const std::string& str, const std::string& tokens) {
	std::vector<std::string> ans;

	size_t pos = 0;
	int posFim;
	while (pos < str.size()) {
		posFim = getNextToken(str, tokens, pos);
		if (pos != std::string::npos) {
			int sz = posFim - pos;
			if (sz)
				ans.emplace_back(str.substr(pos, sz));
			pos = posFim + 1;
		}
		else {
			break;
		}
	}
	// imprime
	// for (auto& i : ans) {
	//     if (i.find("\n") != std::string::npos)
	//         std::cout << "\\n" << ", ";
	//     else
	//         std::cout << i << ", ";
	// }
	return ans;
}

std::vector<std::vector<std::string>> readMatrix(const std::string& name, const std::string& fileStr) {
	//std::string name = "bus";
	std::string regStr = std::string("mpc\\.") + name + std::string("[ \\t]*=[ \\t]*\\[([\\n\\t\\w\\d\\.\\-\\+]+)\\]");
	//std::cout << regStr;
	std::regex r(regStr);

	std::smatch matches;

	std::vector<std::vector<std::string>> ans;
	if (std::regex_search(fileStr, matches, r)) {
		std::cout << "Match found\n";

		auto strMatrix = matches[1].str();
		// apaga quebras de linha extras
		strMatrix = std::regex_replace(strMatrix, std::regex("^\n|\n$"), "");

		auto strMatrixLines = tokenize(strMatrix, "\n");
		//dd::print(strMatrix);
		for (size_t i = 0; i < strMatrixLines.size(); i++) {
			ans.emplace_back(tokenize(strMatrixLines[i], "\t"));
		}
	}
	else {
		std::cout << "Match not found\n";
	}
	// dd::print(ans);
	return ans;
}

std::string readValue(const std::string& name, const std::string& fileStr) {
	//std::string name = "bus";
	std::string regStr = std::string("mpc\\.") + name + std::string("[ \\t]*=[ \\t]*'?([\\w\\d\\.\\-\\+]+)'?[ \\t]*\\n");
	//std::cout << regStr;
	std::regex r(regStr);

	std::smatch matches;

	std::string ans;
	if (std::regex_search(fileStr, matches, r)) {
		std::cout << "Match found\n";

		ans = matches[1].str();
		// apaga quebras de linha extras
	/*    strMatrix = std::regex_replace(strMatrix, std::regex("^\n|\n$"), "");

		auto strMatrixLines = tokenize(strMatrix, "\n");
		dd::print(strMatrix);
		for (int i = 0; i < strMatrixLines.size(); i++) {
			ans.emplace_back(tokenize(strMatrixLines[i], "\t"));
		}*/
	}
	else {
		std::cout << "Match not found\n";
	}
	// dd::print(ans);
	return ans;
}

std::vector<double> parseVector(std::vector<std::string> vec) {
	std::vector<double> ans;
	for (auto& i : vec) {
		i = std::regex_replace(i, std::regex("inf", std::regex_constants::icase), INF_VALUE);
		i = std::regex_replace(i, std::regex("nan", std::regex_constants::icase), NAN_VALUE);
		ans.emplace_back(std::stod(i));
	}
	return ans;
}

matPowerDataType lerArquivoMatPOWER(std::string filePath) {
	// Lê arquivo
	std::ifstream ifs(filePath);
	std::string str((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
	ifs.close();

	// apaga comentários
	str = std::regex_replace(str, std::regex("%.*\n"), "\n");

	// apaga ;
	str = std::regex_replace(str, std::regex("; *\n"), "\n");

	// apaga excesso de \n
	str = std::regex_replace(str, std::regex("\n\n+"), "\n");

	auto busData = readMatrix("bus", str);
	// dd::print("busData");
	// dd::print(busData);

	auto genData = readMatrix("gen", str);
	// dd::print("genData");
	// dd::print(genData);

	auto branchData = readMatrix("branch", str);
	// dd::print("branchData");
	// dd::print(branchData);

	auto gencostData = readMatrix("gencost", str);
	// dd::print("gencostData");
	// dd::print(gencostData);


	auto baseMVAData = readValue("baseMVA", str);
	// dd::print("baseMVAData");
	// dd::print(baseMVAData);

	auto versionData = readValue("version", str);
	// dd::print("versionData");
	// dd::print(versionData);

	matPowerDataType mpData;

	for (auto& i : busData) {
		auto vec = parseVector(i);
		mpData.bus.BUS_I.emplace_back((int)vec[0]);
		mpData.bus.BUS_TYPE.emplace_back((int)vec[1]);
		mpData.bus.PD.emplace_back(vec[2]);
		mpData.bus.QD.emplace_back(vec[3]);
		mpData.bus.GS.emplace_back(vec[4]);
		mpData.bus.BS.emplace_back(vec[5]);
		mpData.bus.BUS_AREA.emplace_back((int)vec[6]);
		mpData.bus.VM.emplace_back(vec[7]);
		mpData.bus.VA.emplace_back(vec[8]);
		mpData.bus.BASE_KV.emplace_back(vec[9]);
		mpData.bus.ZONE.emplace_back((int)vec[10]);
		mpData.bus.VMAX.emplace_back(vec[11]);
		mpData.bus.VMIN.emplace_back(vec[12]);
		//mpData.bus.LAM_P.emplace_back(vec[13]);
		//mpData.bus.LAM_Q.emplace_back(vec[14]);
		//mpData.bus.MU_VMAX.emplace_back(vec[15]);
		//mpData.bus.MU_VMIN.emplace_back(vec[16]);
	}

	for (auto& i : genData) {
		auto vec = parseVector(i);
		mpData.gen.GEN_BUS.emplace_back((int)vec[0]);
		mpData.gen.PG.emplace_back(vec[1]);
		mpData.gen.QG.emplace_back(vec[2]);
		mpData.gen.QMAX.emplace_back(vec[3]);
		mpData.gen.QMIN.emplace_back(vec[4]);
		mpData.gen.VG.emplace_back(vec[5]);
		mpData.gen.MBASE.emplace_back(vec[6]);
		mpData.gen.GEN_STATUS.emplace_back(vec[7]);
		mpData.gen.PMAX.emplace_back(vec[8]);
		mpData.gen.PMIN.emplace_back(vec[9]);
		mpData.gen.PC1.emplace_back(vec[10]);
		mpData.gen.PC2.emplace_back(vec[11]);
		mpData.gen.QC1MIN.emplace_back(vec[12]);
		mpData.gen.QC1MAX.emplace_back(vec[13]);
		mpData.gen.QC2MIN.emplace_back(vec[14]);
		mpData.gen.QC2MAX.emplace_back(vec[15]);
		mpData.gen.RAMP_AGC.emplace_back(vec[16]);
		mpData.gen.RAMP_10.emplace_back(vec[17]);
		mpData.gen.RAMP_30.emplace_back(vec[18]);
		mpData.gen.RAMP_Q.emplace_back(vec[19]);
		mpData.gen.APF.emplace_back(vec[20]);
		//mpData.gen.MU_PMAX.emplace_back(vec[21]);
		//mpData.gen.MU_PMIN.emplace_back(vec[22]);
		//mpData.gen.MU_QMAX.emplace_back(vec[23]);
		//mpData.gen.MU_QMIN.emplace_back(vec[24]);
	}

	for (auto& i : branchData) {
		auto vec = parseVector(i);
		mpData.branch.F_BUS.emplace_back((int)vec[0]);
		mpData.branch.T_BUS.emplace_back((int)vec[1]);
		mpData.branch.BR_R.emplace_back(vec[2]);
		mpData.branch.BR_X.emplace_back(vec[3]);
		mpData.branch.BR_B.emplace_back(vec[4]);
		mpData.branch.RATE_A.emplace_back(vec[5]);
		mpData.branch.RATE_B.emplace_back(vec[6]);
		mpData.branch.RATE_C.emplace_back(vec[7]);
		mpData.branch.TAP.emplace_back(vec[8]);
		mpData.branch.SHIFT.emplace_back(vec[9]);
		mpData.branch.BR_STATUS.emplace_back(vec[10]);
		mpData.branch.ANGMIN.emplace_back(vec[11]);
		mpData.branch.ANGMAX.emplace_back(vec[12]);
		//mpData.branch.PF.emplace_back(vec[13]);
		//mpData.branch.QF.emplace_back(vec[14]);
		//mpData.branch.PT.emplace_back(vec[15]);
		//mpData.branch.QT.emplace_back(vec[16]);
		//mpData.branch.MU_SF.emplace_back(vec[17]);
		//mpData.branch.MU_ST.emplace_back(vec[18]);
		//mpData.branch.MU_ANGMIN.emplace_back(vec[19]);
		//mpData.branch.MU_ANGMAX.emplace_back(vec[20]);
	}

	// incluir dados de genCost...

	mpData.version = std::stod(versionData);
	mpData.baseMVA = std::stod(baseMVAData);

	return mpData;
}

void cleanLine(std::string& line) {
	// apaga comentários
	line = std::regex_replace(line, std::regex("%.*"), "");
	// apaga ;
	line = std::regex_replace(line, std::regex("[ \\t]*;?[ \\t]*$"), "");
}

std::vector<std::vector<std::string>> readMatrix(std::ifstream& file, std::string& line) {
	std::string regEnd = std::string("\\]"); // ([\\t\\w\\d\\.\\ - \\ + ] + )
	std::regex r(regEnd);

	std::vector<std::vector<std::string>> ans;

	bool isLastLine;

	while (true) {
		isLastLine = std::regex_search(line, r);

		if (isLastLine) {
			line = std::regex_replace(line, r, "");
		}

		cleanLine(line);

		bool linhaVazia = line == "";
		if (!linhaVazia) {
			ans.emplace_back(tokenize(line, "\t"));
		}

		if (isLastLine) {
			// dd::print(ans);
			return ans;
		}

		if (!std::getline(file, line)) {
			break;
		}
	}
	return ans;
}

matPowerDataType lerArquivoMatPOWEReficiente(const std::string& filePath) {
	std::vector<std::vector<std::string>> busData;
	std::vector<std::vector<std::string>> genData;
	std::vector<std::vector<std::string>> branchData;
	std::vector<std::vector<std::string>> gencostData;

	std::string baseMVAData;
	std::string versionData;
	
	std::ifstream file(filePath);
	if (file.is_open()) {
		std::string line;
		while (std::getline(file, line)) {
			cleanLine(line);

			bool linhaVazia = line == ""; // std::regex_match(line, std::regex(""));
			if (linhaVazia) continue;

			std::vector<std::string> names{ "bus", "gen", "branch" };
			std::vector<std::vector<std::vector<std::string>>*> vars{ &busData, &genData, &branchData };

			for (size_t i = 0; i < names.size(); i++) {
				std::string regStr = std::string("mpc\\.") + names[i] + std::string("[ \\t]*=[ \\t]*\\[");
				std::regex r(regStr);
				if (std::regex_search(line, r)) {
					line = std::regex_replace(line, r, "");
					*(vars[i]) = readMatrix(file, line);
					break;
				}
			}

			std::vector<std::string> scalarNames{ "baseMVA", "version" };
			std::vector<std::string*> scalarVars{ &baseMVAData, &versionData };
			
			for (size_t i = 0; i < scalarNames.size(); i++) {
				std::string regStr = std::string("mpc\\.") + scalarNames[i] + std::string("[ \\t]*=[ \\t]*");
				std::regex r(regStr);
				if (std::regex_search(line, r)) {
					line = std::regex_replace(line, r, "");
					line = std::regex_replace(line, std::regex("[\"\']"), "");
					*(scalarVars[i]) = line;
					break;
				}
			}
		}
		file.close();
	}
	else {
		std::cout << "Erro na leitura do arquivo. Ele existe?\n";
		matPowerDataType mpData;
		return mpData;
	}

	matPowerDataType mpData;

	for (auto& i : busData) {
		auto vec = parseVector(i);
		mpData.bus.BUS_I.emplace_back((int)vec[0]);
		mpData.bus.BUS_TYPE.emplace_back((int)vec[1]);
		mpData.bus.PD.emplace_back(vec[2]);
		mpData.bus.QD.emplace_back(vec[3]);
		mpData.bus.GS.emplace_back(vec[4]);
		mpData.bus.BS.emplace_back(vec[5]);
		mpData.bus.BUS_AREA.emplace_back((int)vec[6]);
		mpData.bus.VM.emplace_back(vec[7]);
		mpData.bus.VA.emplace_back(vec[8]);
		mpData.bus.BASE_KV.emplace_back(vec[9]);
		mpData.bus.ZONE.emplace_back((int)vec[10]);
		mpData.bus.VMAX.emplace_back(vec[11]);
		mpData.bus.VMIN.emplace_back(vec[12]);
		//mpData.bus.LAM_P.emplace_back(vec[13]);
		//mpData.bus.LAM_Q.emplace_back(vec[14]);
		//mpData.bus.MU_VMAX.emplace_back(vec[15]);
		//mpData.bus.MU_VMIN.emplace_back(vec[16]);
	}

	for (auto& i : genData) {
		auto vec = parseVector(i);
		mpData.gen.GEN_BUS.emplace_back((int)vec[0]);
		mpData.gen.PG.emplace_back(vec[1]);
		mpData.gen.QG.emplace_back(vec[2]);
		mpData.gen.QMAX.emplace_back(vec[3]);
		mpData.gen.QMIN.emplace_back(vec[4]);
		mpData.gen.VG.emplace_back(vec[5]);
		mpData.gen.MBASE.emplace_back(vec[6]);
		mpData.gen.GEN_STATUS.emplace_back(vec[7]);
		mpData.gen.PMAX.emplace_back(vec[8]);
		mpData.gen.PMIN.emplace_back(vec[9]);
		mpData.gen.PC1.emplace_back(vec[10]);
		mpData.gen.PC2.emplace_back(vec[11]);
		mpData.gen.QC1MIN.emplace_back(vec[12]);
		mpData.gen.QC1MAX.emplace_back(vec[13]);
		mpData.gen.QC2MIN.emplace_back(vec[14]);
		mpData.gen.QC2MAX.emplace_back(vec[15]);
		mpData.gen.RAMP_AGC.emplace_back(vec[16]);
		mpData.gen.RAMP_10.emplace_back(vec[17]);
		mpData.gen.RAMP_30.emplace_back(vec[18]);
		mpData.gen.RAMP_Q.emplace_back(vec[19]);
		mpData.gen.APF.emplace_back(vec[20]);
		//mpData.gen.MU_PMAX.emplace_back(vec[21]);
		//mpData.gen.MU_PMIN.emplace_back(vec[22]);
		//mpData.gen.MU_QMAX.emplace_back(vec[23]);
		//mpData.gen.MU_QMIN.emplace_back(vec[24]);
	}

	for (auto& i : branchData) {
		auto vec = parseVector(i);
		mpData.branch.F_BUS.emplace_back((int)vec[0]);
		mpData.branch.T_BUS.emplace_back((int)vec[1]);
		mpData.branch.BR_R.emplace_back(vec[2]);
		mpData.branch.BR_X.emplace_back(vec[3]);
		mpData.branch.BR_B.emplace_back(vec[4]);
		mpData.branch.RATE_A.emplace_back(vec[5]);
		mpData.branch.RATE_B.emplace_back(vec[6]);
		mpData.branch.RATE_C.emplace_back(vec[7]);
		mpData.branch.TAP.emplace_back(vec[8]);
		mpData.branch.SHIFT.emplace_back(vec[9]);
		mpData.branch.BR_STATUS.emplace_back(vec[10]);
		mpData.branch.ANGMIN.emplace_back(vec[11]);
		mpData.branch.ANGMAX.emplace_back(vec[12]);
		//mpData.branch.PF.emplace_back(vec[13]);
		//mpData.branch.QF.emplace_back(vec[14]);
		//mpData.branch.PT.emplace_back(vec[15]);
		//mpData.branch.QT.emplace_back(vec[16]);
		//mpData.branch.MU_SF.emplace_back(vec[17]);
		//mpData.branch.MU_ST.emplace_back(vec[18]);
		//mpData.branch.MU_ANGMIN.emplace_back(vec[19]);
		//mpData.branch.MU_ANGMAX.emplace_back(vec[20]);
	}

	// incluir dados de genCost...

	mpData.version = std::stod(versionData);
	mpData.baseMVA = std::stod(baseMVAData);

	return mpData;
}

bool storeMatPOWER(sistemaType& sistema, barraType& barra, ramoType& ramo, matPowerDataType& mpData) {
	sistema.nB = 0;
	sistema.nPQ = 0;
	sistema.nPV = 0;
	sistema.barraVO = 0;

	sistema.baseMVA = mpData.baseMVA;

	for (size_t i = 1; i <= mpData.bus.BUS_I.size(); i++) {
		
		long long int iGen = -1;

		barra.id[IDX1F(i)] = mpData.bus.BUS_I[IDX1F(i)];

		//barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[i]];  // Base kV // repetida...

		switch (mpData.bus.BUS_TYPE[IDX1F(i)]) {
		case 1: // tipo PQ [Unregulated (load, PQ)]
			sistema.barrasPQ[sistema.nPQ] = i;
			sistema.nPQ++;
			break;
		case 2: // tipo PV [Hold voltage within VAR limits (gen, PV)]
			sistema.barrasPV[sistema.nPV] = i;
			sistema.nPV++;

			barra.V[IDX1F(i)] = mpData.bus.VM[IDX1F(i)];

			// get generator number
			{
				std::vector<int>::iterator it = std::find(mpData.gen.GEN_BUS.begin(), mpData.gen.GEN_BUS.end(), barra.id[IDX1F(i)]);
				if (it == mpData.gen.GEN_BUS.end()) {
					std::cout << "barra " << barra.id[IDX1F(i)] << ", de tipo PV, não possui de geração!\n";
				}
				else {
					iGen = std::distance(mpData.gen.GEN_BUS.begin(), it); // estudar uso de nPV...

					//Limites
					sistema.limQsup[IDX1F(sistema.nPV)] = mpData.gen.QMAX[iGen] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
					sistema.limQinf[IDX1F(sistema.nPV)] = mpData.gen.QMIN[iGen] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]
				}
			}

			break;
		case 3: // tipo VO [Hold voltage and angle (swing, V-Theta)]
			if (sistema.barraVO == 0) {
				sistema.barraVO = i;

				barra.V[IDX1F(i)] = mpData.bus.VM[IDX1F(i)];

				barra.theta[IDX1F(i)] = mpData.bus.VA[IDX1F(i)] * 3.14159265358979323846 / 180.;  // Final angle, degrees

				{
					std::vector<int>::iterator it = std::find(mpData.gen.GEN_BUS.begin(), mpData.gen.GEN_BUS.end(), barra.id[IDX1F(i)]);
					if (it == mpData.gen.GEN_BUS.end()) {
						std::cout << "barra " << barra.id[IDX1F(i)] << ", de tipo VO, não possui de geração!\n";
					}
					else {
						iGen = std::distance(mpData.gen.GEN_BUS.begin(), it);
					}
				}
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
		barra.Pload[IDX1F(i)] = mpData.bus.PD[IDX1F(i)] / sistema.baseMVA;  // Load MW 
		barra.Qload[IDX1F(i)] = mpData.bus.QD[IDX1F(i)] / sistema.baseMVA; // Load MVAR
		barra.Pg[IDX1F(i)] = iGen == -1 ? 0 : mpData.gen.PG[iGen] / sistema.baseMVA;  // Generation MW // [59, 66]
		barra.Qg[IDX1F(i)] = iGen == -1 ? 0 : mpData.gen.QG[iGen] / sistema.baseMVA;  // Generation MVAR // [67, 74]
		barra.Vbase[IDX1F(i)] = mpData.bus.BASE_KV[IDX1F(i)];  // Base KV
		// onde estao esses valores sh??
		barra.bsh[IDX1F(i)] = mpData.bus.BS[IDX1F(i)] / sistema.baseMVA; // 
		barra.gsh[IDX1F(i)] = mpData.bus.GS[IDX1F(i)] / sistema.baseMVA; // Generation MVAR // [67, 74]

		barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
		barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos
	}

	// número de barras
	sistema.nB = sistema.nPQ + sistema.nPV + 1;

	sistema.nL = 0;
	// int nRamosDuplicatas = 0;
	for (size_t i = 0; i < mpData.branch.F_BUS.size(); i++) // i itera ramos
	{
		int auxde = id2i(mpData.branch.F_BUS[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
		int auxpara = id2i(mpData.branch.T_BUS[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]

		//if (i == 50)
		//{
		//	std::cin.get();
		//}

		bool flgRamoNovo = 1;
		int j=0; // j aponta número do ramo repetido

		//for (j = 0; j < sistema.nL; j++) {
		//	if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
		//		// ramo atual é o mesmo que o i-ésimo ramo
		//		// soma elementos em paralelo

		//		flgRamoNovo = 0;
		//		nRamosDuplicatas++;
		//		break;
		//	}
		//}
		if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
			ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(mpData.branch.BR_R[i],  // ramo resistance R, per unit [19, 28] ,
				mpData.branch.BR_X[i]); // ramo reactance  X, per unit [29, 39]

			float_type auxA = mpData.branch.TAP[i]; // Transformer final turns ratio (F) // em p.u.?
			ramo.phi[IDX1F(sistema.nL)] = mpData.branch.SHIFT[i] * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)
			
			if (auxA != 0.) {
				ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
			}
			else {
				ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
			}

			ramo.bsh[IDX1F(sistema.nL)] = mpData.branch.BR_B[i]; // Line charging B, per unit (F) * (total line charging, +B)
		}
		else { // ramo i é duplicata do nL lido
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[j])], barra.id[IDX1F(ramo.para[j])]);
			}

			// conferência do tap
			float_type auxA = mpData.branch.TAP[i]; // tap
			float_type auxPhi = mpData.branch.SHIFT[i] * 3.14159265358979323846 / 180.;
			complex_type aux;
			if (auxA != 0.) {
				aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); // a * cis phi
			}
			else {
				aux = _mkComplex(cos(auxPhi), sin(auxPhi));
			}

			// se valor lido é diferente do já armazenado: atenção!
			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %zu)\n", i + 1);
				}
			}
			// impedância série
			// aux = _mkComplex(.0, .0);

			aux = _mkComplex(mpData.branch.BR_R[i],  // ramo resistance R, per unit [19, 28] ,
				mpData.branch.BR_X[i]); // ramo reactance  X, per unit [29, 39]

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)

			// bsh
			float_type aux2 = mpData.branch.BR_B[i];   // Line charging B, per unit (F) * (total line charging, +B)
			if (aux2 != 0 || ramo.bsh[j] != 0)
				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
			else
				ramo.bsh[j] = 0;
		}
	}

	//for (int i = 0; i < DCSC.daBarra.size(); i++) // i itera Compensadores Série Controláveis
	//{
	//	int auxde = id2i(DCSC.daBarra[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
	//	int auxpara = id2i(DCSC.paraBarra[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
	//	bool flgRamoNovo = 1;
	//	int j; // j aponta número do ramo repetido

	//	for (j = 0; j < sistema.nL; j++) {
	//		if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
	//			// ramo atual é o mesmo que o i-ésimo ramo
	//			// soma elementos em paralelo

	//			flgRamoNovo = 0;
	//			nRamosDuplicatas++;
	//			break;
	//		}
	//	}
	//	if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
	//		sistema.nL++;
	//		ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
	//		ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);

	//		ramo.z[IDX1F(sistema.nL)] = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
	//			DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]

	////float_type auxA = DLIN.tap[i]; // Transformer final turns ratio (F) // em p.u.?
	//		ramo.phi[IDX1F(sistema.nL)] = 0.0; // Transformer (phase shifter) final angle (F)

	//		ramo.tap[IDX1F(sistema.nL)] = _mkComplex(1.0, 0.0); // a * cis phi

	//		ramo.bsh[IDX1F(sistema.nL)] = 0.0; // Line charging B, per unit (F) * (total line charging, +B)
	//	}
	//	else { // ramo i é duplicata do nL lido
	//		if (!global::laconic_mode) {
	//			printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
	//		}

	//		// conferência do tap
	//		//float_type auxA = DLIN.tap[i]; // tap
	//		//float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
	//		complex_type aux;
	//		aux = _mkComplex(1.0, 0.0); // a * cis phi

	//		// se valor lido é diferente do já armazenado: atenção!
	//		if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
	//			if (!global::laconic_mode) {
	//				printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);
	//			}
	//		}
	//		// impedância série
	//		// aux = _mkComplex(.0, .0);

	//		aux = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
	//			DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]

	//		ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)

	//		// bsh
	//		float_type aux2 = 0.0;   // Line charging B, per unit (F) * (total line charging, +B)
	//		if (aux2 != 0 || ramo.bsh[j] != 0)
	//			ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
	//		else
	//			ramo.bsh[j] = 0;
	//	}
	//}

	return 0;
}

matPowerDataType lerMatPowerEAlocarMemoria(std::string mFile, sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	// entrada
	matPowerDataType mpData = lerArquivoMatPOWEReficiente(mFile);

	sistema.nB = mpData.bus.BUS_I.size();
	sistema.nL = mpData.branch.F_BUS.size();
	sistema.nPQ = std::count(mpData.bus.BUS_TYPE.begin(), mpData.bus.BUS_TYPE.end(), 1);
	sistema.nPV = std::count(mpData.bus.BUS_TYPE.begin(), mpData.bus.BUS_TYPE.end(), 2);

	initBranch(sistema, ramo);
	initSistema(sistema);
	initBus(sistema, barra);
	initIter(sistema, iterativo);


	if (storeMatPOWER(sistema, barra, ramo, mpData)) { // guarda elementos lidos nas estruturas do Flumen
		printf("Erro ao abrir o arquivo...");
	}
	return mpData;
}

//bool lerTamanhosMatPOWER(std::string cdfFile, sistemaType &sistema){
//
//	sistema.nB = 0;
//	sistema.nPQ = 0;
//	sistema.nPV = 0;
//	sistema.barraVO = 0;
//
//	sistema.baseMVA = DCTE.BASE;
//
//	for (int i = 1; i <= DBAR.numero.size(); i++) {
//
//		barra.id[IDX1F(i)] = DBAR.numero[IDX1F(i)];
//
//		//barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[i]];  // Base kV // repetida...
//
//		switch (DBAR.tipo[IDX1F(i)]) {
//		case '3':
//			std::cout << "AVISO: Barra com limites de tensão inserida" << std::endl; //   barra  de  carga  com  limite  de  tensão  (PQ  -  Injeções  de  potências  ativa  e  reativa  fixas  enquanto  a  magnitude  de  tensão  permanecer  entre os valores limites)
//		case '0': // tipo PQ [Unregulated (load, PQ)]
//			sistema.barrasPQ[sistema.nPQ] = i;
//			sistema.nPQ++;
//			break;
//		case '1': // tipo PV [Hold voltage within VAR limits (gen, PV)]
//			sistema.barrasPV[sistema.nPV] = i;
//			sistema.nPV++;
//
//			barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];
//
//			//Limites
//			sistema.limQsup[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMaxima[IDX1F(i)] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
//			sistema.limQinf[IDX1F(sistema.nPV)] = DBAR.geracaoReativaMinima[IDX1F(i)] / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]
//
//			break;
//		case '2':// tipo VO [Hold voltage and angle (swing, V-Theta)]
//			if (sistema.barraVO == 0) {
//				sistema.barraVO = i;
//
//				barra.V[IDX1F(i)] = DBAR.tensao[IDX1F(i)];
//
//				barra.theta[IDX1F(i)] = DBAR.angulo[IDX1F(i)] * 3.14159265358979323846 / 180.;  // Final angle, degrees
//			}
//			else {
//				std::cout << "ERRO: mais de uma barrra swing definida!" << std::endl;
//				return 1;
//			}
//			break;
//		default:
//			std::cout << "ERRO: valor de tipo de barra inválido!" << std::endl;
//			return 1;
//		}
//		//std::string auxx = line.substr(40, 9).c_str();
//		barra.Pload[IDX1F(i)] = DBAR.cargaAtiva[IDX1F(i)] / sistema.baseMVA;  // Load MW 
//		barra.Qload[IDX1F(i)] = DBAR.cargaReativa[IDX1F(i)] / sistema.baseMVA; // Load MVAR
//		barra.Pg[IDX1F(i)] = DBAR.geracaoAtiva[IDX1F(i)] / sistema.baseMVA;  // Generation MW // [59, 66]
//		barra.Qg[IDX1F(i)] = DBAR.geracaoReativa[IDX1F(i)] / sistema.baseMVA;  // Generation MVAR // [67, 74]
//		barra.Vbase[IDX1F(i)] = DGBT.mapa[DBAR.grupoDeBaseDeTensao[IDX1F(i)]];  // Base KV
//		// onde estao esses valores sh??
//		barra.bsh[IDX1F(i)] = 0; // 
//		barra.gsh[IDX1F(i)] = 0; // Generation MVAR // [67, 74]
//
//		barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
//		barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos
//	}
//
//	// número de barras
//	sistema.nB = sistema.nPQ + sistema.nPV + 1;
//
//	sistema.nL = 0;
//	int nRamosDuplicatas = 0;
//	for (int i = 0; i < DLIN.daBarra.size(); i++) // i itera ramos
//	{
//		int auxde = id2i(DLIN.daBarra[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
//		int auxpara = id2i(DLIN.paraBarra[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
//		bool flgRamoNovo = 1;
//		int j; // j aponta número do ramo repetido
//
//		for (j = 0; j < sistema.nL; j++) {
//			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
//				// ramo atual é o mesmo que o i-ésimo ramo
//				// soma elementos em paralelo
//
//				flgRamoNovo = 0;
//				nRamosDuplicatas++;
//				break;
//			}
//		}
//		if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
//			sistema.nL++;
//			ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
//			ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);
//
//			ramo.z[IDX1F(sistema.nL)] = _mkComplex(DLIN.resistencia[i] / 100,  // ramo resistance R, per unit [19, 28] ,
//				DLIN.reatancia[i] / 100); // ramo reactance  X, per unit [29, 39]
//
//			float_type auxA = DLIN.tap[i]; // Transformer final turns ratio (F) // em p.u.?
//			ramo.phi[IDX1F(sistema.nL)] = DLIN.defasagem[i] * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)
//
//			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
//
//			ramo.bsh[IDX1F(sistema.nL)] = DLIN.susceptancia[i] / sistema.baseMVA; // Line charging B, per unit (F) * (total line charging, +B)
//		}
//		else { // ramo i é duplicata do nL lido
//			if (!global::laconic_mode) {
//				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
//			}
//
//			// conferência do tap
//			float_type auxA = DLIN.tap[i]; // tap
//			float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
//			complex_type aux;
//			aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); // a * cis phi
//
//			// se valor lido é diferente do já armazenado: atenção!
//			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
//				if (!global::laconic_mode) {
//					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);
//				}
//			}
//			// impedância série
//			// aux = _mkComplex(.0, .0);
//
//			aux = _mkComplex(DLIN.resistencia[i] / 100,  // ramo resistance R, per unit [19, 28] ,
//				DLIN.reatancia[i] / 100); // ramo reactance  X, per unit [29, 39]
//
//			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)
//
//			// bsh
//			float_type aux2 = DLIN.susceptancia[i];   // Line charging B, per unit (F) * (total line charging, +B)
//			if (aux2 != 0 || ramo.bsh[j] != 0)
//				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
//			else
//				ramo.bsh[j] = 0;
//		}
//	}
//
//	for (int i = 0; i < DCSC.daBarra.size(); i++) // i itera Compensadores Série Controláveis
//	{
//		int auxde = id2i(DCSC.daBarra[i], sistema, barra); // atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
//		int auxpara = id2i(DCSC.paraBarra[i], sistema, barra); // atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
//		bool flgRamoNovo = 1;
//		int j; // j aponta número do ramo repetido
//
//		for (j = 0; j < sistema.nL; j++) {
//			if ((ramo.de[j] == auxde) && (ramo.para[j] == auxpara)) {
//				// ramo atual é o mesmo que o i-ésimo ramo
//				// soma elementos em paralelo
//
//				flgRamoNovo = 0;
//				nRamosDuplicatas++;
//				break;
//			}
//		}
//		if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
//			sistema.nL++;
//			ramo.de[IDX1F(sistema.nL)] = auxde; // id2i(auxde, sistema, barra);
//			ramo.para[IDX1F(sistema.nL)] = auxpara; // id2i(auxpara, sistema, barra);
//
//			ramo.z[IDX1F(sistema.nL)] = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
//				DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]
//
////float_type auxA = DLIN.tap[i]; // Transformer final turns ratio (F) // em p.u.?
//			ramo.phi[IDX1F(sistema.nL)] = 0.0; // Transformer (phase shifter) final angle (F)
//
//			ramo.tap[IDX1F(sistema.nL)] = _mkComplex(1.0, 0.0); // a * cis phi
//
//			ramo.bsh[IDX1F(sistema.nL)] = 0.0; // Line charging B, per unit (F) * (total line charging, +B)
//		}
//		else { // ramo i é duplicata do nL lido
//			if (!global::laconic_mode) {
//				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[IDX1F(sistema.nL)])], barra.id[IDX1F(ramo.para[IDX1F(sistema.nL)])]);
//			}
//
//			// conferência do tap
//			//float_type auxA = DLIN.tap[i]; // tap
//			//float_type auxPhi = DLIN.defasagem[i] * 3.14159265358979323846 / 180.;
//			complex_type aux;
//			aux = _mkComplex(1.0, 0.0); // a * cis phi
//
//			// se valor lido é diferente do já armazenado: atenção!
//			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
//				if (!global::laconic_mode) {
//					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);
//				}
//			}
//			// impedância série
//			// aux = _mkComplex(.0, .0);
//
//			aux = _mkComplex(0.0,  // ramo resistance R, per unit [19, 28] ,
//				DCSC.valorEspecificado[i] / 100); // ramo reactance  X, per unit [29, 39]
//
//			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); // A//B = (A*B)/(A+B)
//
//			// bsh
//			float_type aux2 = 0.0;   // Line charging B, per unit (F) * (total line charging, +B)
//			if (aux2 != 0 || ramo.bsh[j] != 0)
//				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
//			else
//				ramo.bsh[j] = 0;
//		}
//	}
//
//	return 0;
//}
//
//bool readMatPOWER(std::string cdfFile, sistemaType& sistema, barraType& barra, ramoType& ramo) {
//
//	std::string line;
//	std::fstream MatPOWER;
//	MatPOWER.open(cdfFile.c_str(), std::ios::in);
//	if (MatPOWER.is_open())
//	{
//		sistema.nB = 0;
//		sistema.nPQ = 0;
//		sistema.nPV = 0;
//		sistema.barraVO = 0;
//		std::getline(MatPOWER, line);
//
//		//std::cout << "linha: " << atoi(line.substr(31, 6).c_str()) << std::endl;
//		sistema.baseMVA = atoi(line.substr(31, 6).c_str());
//
//		int i = 0; // i é o número da barra que está sendo lida
//
//		std::getline(MatPOWER, line); // pula duas linhas
//		while (std::getline(MatPOWER, line)) // atualiza line a cada itera��o
//		{
//			if (line.find("-999") == std::string::npos) {
//				i++; // nova barra encontrada
//
//				barra.id[IDX1F(i)] = atoi(line.substr(0, 4).c_str()); // [0, 3]
//
//				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]
//
//				switch (atoi(line.substr(24, 2).c_str())) { // [24, 25]
//				case 0: // tipo PQ [Unregulated (load, PQ)]
//					sistema.barrasPQ[sistema.nPQ] = i;
//					sistema.nPQ++;
//					barra.pntPQ[IDX1F(i)] = sistema.nPQ;
//					break;
//				case 1: // tipo PQ [Hold MVAR generation within voltage limits, (PQ)]
//					sistema.barrasPQ[sistema.nPQ] = i;
//					sistema.nPQ++;
//					barra.pntPQ[IDX1F(i)] = sistema.nPQ;
//					break;
//				case 2: // tipo PV [Hold voltage within VAR limits (gen, PV)]
//					sistema.barrasPV[sistema.nPV] = i;
//					sistema.nPV++;
//
//					barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
//					sistema.VfixadaPV[IDX1F(sistema.nPV)] = barra.V[IDX1F(i)]; // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat
//
//					//Limites
//					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(90, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
//					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(98, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]
//
//					break;
//				case 3:// tipo VO [Hold voltage and angle (swing, V-Theta)]
//					if (sistema.barraVO == 0) {
//						sistema.barraVO = i;
//
//						barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
//
//						barra.theta[IDX1F(i)] = atof(line.substr(33, 7).c_str()) * 3.14159265358979323846 / 180.;  // Final angle, degrees
//					}
//					else {
//						std::cout << "ERRO: mais de uma barrra swing definida!" << std::endl;
//						return 1;
//					}
//					break;
//				default:
//					std::cout << "ERRO: valor de tipo de barra inválido!" << std::endl;
//					return 1;
//				}
//				//std::string auxx = line.substr(40, 9).c_str();
//				barra.Pload[IDX1F(i)] = atof(line.substr(40, 9).c_str()) / sistema.baseMVA;  // Load MW // [40, 48]
//				barra.Qload[IDX1F(i)] = atof(line.substr(49, 10).c_str()) / sistema.baseMVA; // Load MVAR // [49, 58]
//				barra.Pg[IDX1F(i)] = atof(line.substr(59, 8).c_str()) / sistema.baseMVA;  // Generation MW // [59, 66]
//				barra.Qg[IDX1F(i)] = atof(line.substr(67, 8).c_str()) / sistema.baseMVA;  // Generation MVAR // [67, 74]
//				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]
//				barra.bsh[IDX1F(i)] = atof(line.substr(114, 8).c_str()); // 
//				barra.gsh[IDX1F(i)] = atof(line.substr(106, 8).c_str()); // Generation MVAR // [67, 74]
//
//				barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
//				barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos
//
//			}
//			else
//				break; // se -999 a se��o Bus Data acabou!
//		}
//
//		// número de barras
//		sistema.nB = sistema.nPQ + sistema.nPV + 1;
//
//		std::getline(MatPOWER, line); // pula duas linhas
//		sistema.nL = 0;
//		int nRamosDuplicatas = 0;
//		while (std::getline(MatPOWER, line)) // atualiza line a cada itera��o
//		{
//			if (line.find("-999") == std::string::npos) {
//				// Parse branch	
//				int auxde = id2i(atoi(line.substr(0, 4).c_str()), sistema, barra); // ramo da linha i // [0, 3]
//				int auxpara = id2i(atoi(line.substr(4, 5).c_str()), sistema, barra); // at� a j // [4, 8]
//				bool flgRamoNovo = 1;
//				int i;
//
//				for (i = 0; i < sistema.nL; i++) {
//					if ((ramo.de[i] == auxde) && (ramo.para[i] == auxpara)) {
//						// ramo atual é o mesmo que o i-ésimo ramo
//						// soma elementos em paralelo
//
//						flgRamoNovo = 0;
//						nRamosDuplicatas++;
//						break;
//					}
//				}
//				if (flgRamoNovo) { // novo ramo deve ser acrescentado à lista
//					sistema.nL++;
//					ramo.de[IDX1F(sistema.nL)] = auxde;
//					ramo.para[IDX1F(sistema.nL)] = auxpara;
//
//					ramo.z[IDX1F(sistema.nL)] = _mkComplex(atof(line.substr(19, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
//						atof(line.substr(29, 11).c_str())); // ramo reactance  X, per unit [29, 39]
//
//					float_type auxA = atof(line.substr(76, 6).c_str()); // Transformer final turns ratio (F)
//					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)
//					if (auxA != 0.) {
//						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
//					}
//					else {
//						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
//					}
//
//					ramo.bsh[IDX1F(sistema.nL)] = atof(line.substr(40, 10).c_str()); // Line charging B, per unit (F) * (total line charging, +B)
//				}
//				else { // ramo i é duplicata do nL lido
//					if (!global::laconic_mode) {
//						printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", ramo.de[IDX1F(sistema.nL)], ramo.para[IDX1F(sistema.nL)]);
//					}
//
//					// conferência do tap
//					float_type auxA = atof(line.substr(76, 6).c_str()); // tap
//					float_type auxPhi = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.;
//					complex_type aux;
//					if ((auxA != 0.)) {
//						aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); // a * cis phi
//					}
//					else {
//						aux = _mkComplex(cos(auxPhi), sin(auxPhi));
//					}
//					// se valor lido é diferente do já armazenado: atenção!
//					if ((aux.x != ramo.tap[i].x) || (aux.y != ramo.tap[i].y))
//						printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);
//
//					// impedância série
//					aux = _mkComplex(.0, .0);
//
//					aux = _mkComplex(atof(line.substr(19, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
//						atof(line.substr(29, 11).c_str())); // ramo reactance  X, per unit [29, 39]
//
//					ramo.z[i] = _cuDiv(_cuMul(aux, ramo.z[i]), _cuAdd(ramo.z[i], aux)); // A//B = (A*B)/(A+B)
//
//					// bsh
//					float_type aux2 = atof(line.substr(40, 10).c_str());   // Line charging B, per unit (F) * (total line charging, +B)
//					if (aux2 != 0 || ramo.bsh[i] != 0)
//						ramo.bsh[i] = 1 / (aux2 + ramo.bsh[i]); //errado:= ramo.bsh[i]*aux2/(ramo.bsh[i]+aux2); // A//B = (A*B)/(A+B)
//					else
//						ramo.bsh[i] = 0;
//				}
//			}
//			else
//				break; // se -999 a se��o ramo Data acabou!
//		}
//		MatPOWER.close();
//	}
//	else {
//		std::cout << "ERRO: arquivo .m não pode ser aberto!" << std::endl;
//		return 1;
//	}
//
//	return 0;
//}