#pragma once
#include <iostream>
#include <string>
#include <regex>

#include <fstream>
#include <streambuf>

// #include "../idx.h"
// #include "../estruturas.h"
// #include "../global.h"
// #include "id2i.h"

// #include "config.h"

struct mpBusType
{
    std::vector<int> BUS_I;      //   bus number (positive integer)
    std::vector<int> BUS_TYPE;   //   bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
    std::vector<double> PD;      //   real power demand (MW)
    std::vector<double> QD;      //   reactive power demand (MVAr)
    std::vector<double> GS;      //   shunt conductance (MW demanded at V = 1.0 p.u.)
    std::vector<double> BS;      //   shunt susceptance (MVAr injected at V = 1.0 p.u.)
    std::vector<int> BUS_AREA;   //   area number (positive integer)
    std::vector<double> VM;      //   voltage magnitude (p.u.)
    std::vector<double> VA;      //   voltage angle (degrees)
    std::vector<double> BASE_KV; //   base voltage (kV)
    std::vector<int> ZONE;       //   loss zone (positive integer)
    std::vector<double> VMAX;    // † maximum voltage magnitude (p.u.)
    std::vector<double> VMIN;    // † minimum voltage magnitude (p.u.)
    std::vector<double> LAM_P;   // † Lagrange multiplier on real power mismatch (u/MW)
    std::vector<double> LAM_Q;   // † Lagrange multiplier on reactive power mismatch (u/MVAr)
    std::vector<double> MU_VMAX; // † Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    std::vector<double> MU_VMIN; // † Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
};

struct mpGenType
{
    std::vector<int> GEN_BUS;       //   1 bus number
    std::vector<double> PG;         //   2 real power output (MW)
    std::vector<double> QG;         //   3 reactive power output (MVAr)
    std::vector<double> QMAX;       //   4 maximum reactive power output (MVAr)
    std::vector<double> QMIN;       //   5 minimum reactive power output (MVAr)
    std::vector<double> VG;         // ‡ 6 voltage magnitude setpoint (p.u.)
    std::vector<double> MBASE;      //   7 total MVA base of machine, defaults to baseMVA
    std::vector<double> GEN_STATUS; //   8 machine status, > 0 = machine in-service, <= 0 = machine out-of-service
    std::vector<double> PMAX;       //   9 maximum real power output (MW)
    std::vector<double> PMIN;       //   10 minimum real power output (MW)
    std::vector<double> PC1;        // * 11 lower real power output of PQ capability curve (MW)
    std::vector<double> PC2;        // * 12 upper real power output of PQ capability curve (MW)
    std::vector<double> QC1MIN;     // * 13 minimum reactive power output at PC1 (MVAr)
    std::vector<double> QC1MAX;     // * 14 maximum reactive power output at PC1 (MVAr)
    std::vector<double> QC2MIN;     // * 15 minimum reactive power output at PC2 (MVAr)
    std::vector<double> QC2MAX;     // * 16 maximum reactive power output at PC2 (MVAr)
    std::vector<double> RAMP_AGC;   // * 17 ramp rate for load following/AGC (MW/min)
    std::vector<double> RAMP_10;    // * 18 ramp rate for 10 minute reserves (MW)
    std::vector<double> RAMP_30;    // * 19 ramp rate for 30 minute reserves (MW)
    std::vector<double> RAMP_Q;     // * 20 ramp rate for reactive power (2 sec timescale) (MVAr/min)
    std::vector<double> APF;        // * 21 area participation factor
    std::vector<double> MU_PMAX;    // † 22 Kuhn-Tucker multiplier on upper Pg limit (u/MW)
    std::vector<double> MU_PMIN;    // † 23 Kuhn-Tucker multiplier on lower Pg limit (u/MW)
    std::vector<double> MU_QMAX;    // † 24 Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
    std::vector<double> MU_QMIN;    // † 25 Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)
};

struct mpBranchType
{
    std::vector<double> F_BUS;     //   1 "from" bus number
    std::vector<double> T_BUS;     //   2 "to" bus number
    std::vector<double> BR_R;      //   3 resistance (p.u.)
    std::vector<double> BR_X;      //   4 reactance (p.u.)
    std::vector<double> BR_B;      //   5 total line charging susceptance (p.u.)
    std::vector<double> RATE_A;    // * 6 MVA rating A (long term rating), set to 0 for unlimited
    std::vector<double> RATE_B;    // * 7 MVA rating B (short term rating), set to 0 for unlimited
    std::vector<double> RATE_C;    // * 8 MVA rating C (emergency rating), set to 0 for unlimited
    std::vector<double> TAP;       //   9 transformer o nominal turns ratio, if non-zero (taps at \from" bus, impedance at \to" bus, i.e. if r = x = b = 0, tap = jVf j jVtj ; tap = 0 used to indicate transmission line rather than transformer, i.e. mathematically equivalent to transformer with tap = 1)
    std::vector<double> SHIFT;     //   10 transformer phase shift angle (degrees), positive ) delay
    std::vector<double> BR_STATUS; //   11 initial branch status, 1 = in-service, 0 = out-of-service
    std::vector<double> ANGMIN;    // † 12 minimum angle dierence, f 􀀀 t (degrees)
    std::vector<double> ANGMAX;    // † 13 maximum angle dierence, f 􀀀 t (degrees)
    std::vector<double> PF;        // ‡ 14 real power injected at \from" bus end (MW)
    std::vector<double> QF;        // ‡ 15 reactive power injected at \from" bus end (MVAr)
    std::vector<double> PT;        // ‡ 16 real power injected at \to" bus end (MW)
    std::vector<double> QT;        // ‡ 17 reactive power injected at \to" bus end (MVAr)
    std::vector<double> MU_SF;     // § 18 Kuhn-Tucker multiplier on MVA limit at \from" bus (u/MVA)
    std::vector<double> MU_ST;     // § 19 Kuhn-Tucker multiplier on MVA limit at \to" bus (u/MVA)
    std::vector<double> MU_ANGMIN; // § 20 Kuhn-Tucker multiplier lower angle dierence limit (u/degree)
    std::vector<double> MU_ANGMAX; // § 21 Kuhn-Tucker multiplier upper angle dierence limit (u/degree)
};

struct matPowerDataType
{
    mpBusType bus;
    mpGenType gen;
    mpBranchType branch;
    double baseMVA;
    double version;
};

matPowerDataType lerArquivoMatPOWER(std::string filePath);

matPowerDataType lerMatPowerEAlocarMemoria(std::string mFile, sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo);

bool lerTamanhosMatPOWER(matPowerDataType data, sistema &sistema);

// bool readMatPOWER(matPowerDataType data, sistema& sistema, barra& barra, ramo& ramo);

bool storeMatPOWER(sistema& sistema, barra& barra, ramo& ramo, matPowerDataType& mpData);

bool isToken(char chr, const std::string& tokens) {
	for (auto& i : tokens) {
		if (chr == i)
			return true;
	}
	return false;
}

int getNextToken(const std::string& str, const std::string& tokens, const int& startPos) {
	int pos = startPos;

	while (pos < str.size()) {
		if (isToken(str[pos], tokens))
			return pos;
		pos++;
	}

	return str.size();
}

std::vector<std::string> tokenize(const std::string& str, const std::string& tokens) {
	std::vector<std::string> ans;

	int pos = 0;
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
		for (int i = 0; i < strMatrixLines.size(); i++) {
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
		// std::cout << i << '\n';
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

	std::cout << versionData << '\n';
	std::cout << baseMVAData << '\n';

	mpData.version = std::stod(versionData);
	mpData.baseMVA = std::stod(baseMVAData);

	return mpData;
}

bool storeMatPOWER(sistema& sistema, barra& barra, ramo& ramo, matPowerDataType& mpData) {
	sistema.nB = 0;
	sistema.nPQ = 0;
	sistema.nPV = 0;
	sistema.barraVO = 0;

	sistema.baseMVA = mpData.baseMVA;

	for (int i = 1; i <= mpData.bus.BUS_I.size(); i++) {
		
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
	for (int i = 0; i < mpData.branch.F_BUS.size(); i++) // i itera ramos
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
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i + 1);
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

	return 0;
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

			for (int i = 0; i < names.size(); i++) {
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
			
			for (int i = 0; i < scalarNames.size(); i++) {
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

matPowerDataType lerMatPowerEAlocarMemoria(std::string mFile, sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo) {
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
		printf("Deu ruim...");
	}
	return mpData;
}