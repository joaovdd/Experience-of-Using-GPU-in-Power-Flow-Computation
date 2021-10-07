#pragma once
#include <iostream>
#include <string>
#include <regex>

#include <fstream>
#include <streambuf>

struct mpBusType
{
    std::vector<int> BUS_I;      
    std::vector<int> BUS_TYPE;   
    std::vector<double> PD;      
    std::vector<double> QD;      
    std::vector<double> GS;      
    std::vector<double> BS;      
    std::vector<int> BUS_AREA;   
    std::vector<double> VM;      
    std::vector<double> VA;      
    std::vector<double> BASE_KV; 
    std::vector<int> ZONE;       
    std::vector<double> VMAX;    
    std::vector<double> VMIN;    
    std::vector<double> LAM_P;   
    std::vector<double> LAM_Q;   
    std::vector<double> MU_VMAX; 
    std::vector<double> MU_VMIN; 
};

struct mpGenType
{
    std::vector<int> GEN_BUS;       
    std::vector<double> PG;         
    std::vector<double> QG;         
    std::vector<double> QMAX;       
    std::vector<double> QMIN;       
    std::vector<double> VG;         
    std::vector<double> MBASE;      
    std::vector<double> GEN_STATUS; 
    std::vector<double> PMAX;       
    std::vector<double> PMIN;       
    std::vector<double> PC1;        
    std::vector<double> PC2;        
    std::vector<double> QC1MIN;     
    std::vector<double> QC1MAX;     
    std::vector<double> QC2MIN;     
    std::vector<double> QC2MAX;     
    std::vector<double> RAMP_AGC;   
    std::vector<double> RAMP_10;    
    std::vector<double> RAMP_30;    
    std::vector<double> RAMP_Q;     
    std::vector<double> APF;        
    std::vector<double> MU_PMAX;    
    std::vector<double> MU_PMIN;    
    std::vector<double> MU_QMAX;    
    std::vector<double> MU_QMIN;    
};

struct mpBranchType
{
    std::vector<double> F_BUS;     
    std::vector<double> T_BUS;     
    std::vector<double> BR_R;      
    std::vector<double> BR_X;      
    std::vector<double> BR_B;      
    std::vector<double> RATE_A;    
    std::vector<double> RATE_B;    
    std::vector<double> RATE_C;    
    std::vector<double> TAP;       
    std::vector<double> SHIFT;     
    std::vector<double> BR_STATUS; 
    std::vector<double> ANGMIN;    
    std::vector<double> ANGMAX;    
    std::vector<double> PF;        
    std::vector<double> QF;        
    std::vector<double> PT;        
    std::vector<double> QT;        
    std::vector<double> MU_SF;     
    std::vector<double> MU_ST;     
    std::vector<double> MU_ANGMIN; 
    std::vector<double> MU_ANGMAX; 
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

	return ans;
}

std::vector<std::vector<std::string>> readMatrix(const std::string& name, const std::string& fileStr) {
	std::string regStr = std::string("mpc\\.") + name + std::string("[ \\t]*=[ \\t]*\\[([\\n\\t\\w\\d\\.\\-\\+]+)\\]");

	std::regex r(regStr);

	std::smatch matches;

	std::vector<std::vector<std::string>> ans;
	if (std::regex_search(fileStr, matches, r)) {
		std::cout << "Match found\n";

		auto strMatrix = matches[1].str();
		
		strMatrix = std::regex_replace(strMatrix, std::regex("^\n|\n$"), "");

		auto strMatrixLines = tokenize(strMatrix, "\n");
		
		for (int i = 0; i < strMatrixLines.size(); i++) {
			ans.emplace_back(tokenize(strMatrixLines[i], "\t"));
		}
	}
	else {
		std::cout << "Match not found\n";
	}

	return ans;
}

std::string readValue(const std::string& name, const std::string& fileStr) {
	std::string regStr = std::string("mpc\\.") + name + std::string("[ \\t]*=[ \\t]*'?([\\w\\d\\.\\-\\+]+)'?[ \\t]*\\n");

	std::regex r(regStr);

	std::smatch matches;

	std::string ans;
	if (std::regex_search(fileStr, matches, r)) {
		std::cout << "Match found\n";

		ans = matches[1].str();
		
	}
	else {
		std::cout << "Match not found\n";
	}

	return ans;
}

std::vector<double> parseVector(std::vector<std::string> vec) {
	std::vector<double> ans;
	for (auto& i : vec) {
		
		ans.emplace_back(std::stod(i));
	}
	return ans;
}

matPowerDataType lerArquivoMatPOWER(std::string filePath) {
	std::ifstream ifs(filePath);
	std::string str((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
	ifs.close();

	str = std::regex_replace(str, std::regex("%.*\n"), "\n");

	str = std::regex_replace(str, std::regex("; *\n"), "\n");

	str = std::regex_replace(str, std::regex("\n\n+"), "\n");

	auto busData = readMatrix("bus", str);

	auto genData = readMatrix("gen", str);

	auto branchData = readMatrix("branch", str);

	auto gencostData = readMatrix("gencost", str);

	auto baseMVAData = readValue("baseMVA", str);

	auto versionData = readValue("version", str);

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
		
	}

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

		switch (mpData.bus.BUS_TYPE[IDX1F(i)]) {
		case 1: 
			sistema.barrasPQ[sistema.nPQ] = i;
			sistema.nPQ++;
			break;
		case 2: 
			sistema.barrasPV[sistema.nPV] = i;
			sistema.nPV++;

			barra.V[IDX1F(i)] = mpData.bus.VM[IDX1F(i)];

			{
				std::vector<int>::iterator it = std::find(mpData.gen.GEN_BUS.begin(), mpData.gen.GEN_BUS.end(), barra.id[IDX1F(i)]);
				if (it == mpData.gen.GEN_BUS.end()) {
					std::cout << "barra " << barra.id[IDX1F(i)] << ", de tipo PV, não possui de geração!\n";
				}
				else {
					iGen = std::distance(mpData.gen.GEN_BUS.begin(), it); 

					sistema.limQsup[IDX1F(sistema.nPV)] = mpData.gen.QMAX[iGen] / sistema.baseMVA;  
					sistema.limQinf[IDX1F(sistema.nPV)] = mpData.gen.QMIN[iGen] / sistema.baseMVA;  
				}
			}

			break;
		case 3: 
			if (sistema.barraVO == 0) {
				sistema.barraVO = i;

				barra.V[IDX1F(i)] = mpData.bus.VM[IDX1F(i)];

				barra.theta[IDX1F(i)] = mpData.bus.VA[IDX1F(i)] * 3.14159265358979323846 / 180.;  

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
		
		barra.Pload[IDX1F(i)] = mpData.bus.PD[IDX1F(i)] / sistema.baseMVA;  
		barra.Qload[IDX1F(i)] = mpData.bus.QD[IDX1F(i)] / sistema.baseMVA; 
		barra.Pg[IDX1F(i)] = iGen == -1 ? 0 : mpData.gen.PG[iGen] / sistema.baseMVA;  
		barra.Qg[IDX1F(i)] = iGen == -1 ? 0 : mpData.gen.QG[iGen] / sistema.baseMVA;  
		barra.Vbase[IDX1F(i)] = mpData.bus.BASE_KV[IDX1F(i)];  
		
		barra.bsh[IDX1F(i)] = mpData.bus.BS[IDX1F(i)] / sistema.baseMVA; 
		barra.gsh[IDX1F(i)] = mpData.bus.GS[IDX1F(i)] / sistema.baseMVA; 

		barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  
		barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  
	}

	sistema.nB = sistema.nPQ + sistema.nPV + 1;

	sistema.nL = 0;

	for (int i = 0; i < mpData.branch.F_BUS.size(); i++) 
	{
		int auxde = id2i(mpData.branch.F_BUS[i], sistema, barra); 
		int auxpara = id2i(mpData.branch.T_BUS[i], sistema, barra); 

		bool flgRamoNovo = 1;
		int j=0; 

		if (flgRamoNovo) { 
			sistema.nL++;
			ramo.de[IDX1F(sistema.nL)] = auxde; 
			ramo.para[IDX1F(sistema.nL)] = auxpara; 

			ramo.z[IDX1F(sistema.nL)] = _mkComplex(mpData.branch.BR_R[i],  
				mpData.branch.BR_X[i]); 

			float_type auxA = mpData.branch.TAP[i]; 
			ramo.phi[IDX1F(sistema.nL)] = mpData.branch.SHIFT[i] * 3.14159265358979323846 / 180.; 
			
			if (auxA != 0.) {
				ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); 
			}
			else {
				ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
			}

			ramo.bsh[IDX1F(sistema.nL)] = mpData.branch.BR_B[i]; 
		}
		else { 
			if (!global::laconic_mode) {
				printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", barra.id[IDX1F(ramo.de[j])], barra.id[IDX1F(ramo.para[j])]);
			}

			float_type auxA = mpData.branch.TAP[i]; 
			float_type auxPhi = mpData.branch.SHIFT[i] * 3.14159265358979323846 / 180.;
			complex_type aux;
			if (auxA != 0.) {
				aux = _mkComplex(auxA * cos(auxPhi), auxA * sin(auxPhi)); 
			}
			else {
				aux = _mkComplex(cos(auxPhi), sin(auxPhi));
			}

			if ((aux.x != ramo.tap[j].x) || (aux.y != ramo.tap[j].y)) {
				if (!global::laconic_mode) {
					printf("\nAtencao! Os parametros lidos para os transformadores de potencia das linhas de transmissão paralelas são diferentes! (linha %d)\n", i + 1);
				}
			}
			
			aux = _mkComplex(mpData.branch.BR_R[i],  
				mpData.branch.BR_X[i]); 

			ramo.z[j] = _cuDiv(_cuMul(aux, ramo.z[j]), _cuAdd(ramo.z[j], aux)); 

			float_type aux2 = mpData.branch.BR_B[i];   
			if (aux2 != 0 || ramo.bsh[j] != 0)
				ramo.bsh[j] = 1 / (aux2 + ramo.bsh[j]);
			else
				ramo.bsh[j] = 0;
		}
	}

	return 0;
}

void cleanLine(std::string& line) {
	line = std::regex_replace(line, std::regex("%.*"), "");

	line = std::regex_replace(line, std::regex("[ \\t]*;?[ \\t]*$"), "");
}

std::vector<std::vector<std::string>> readMatrix(std::ifstream& file, std::string& line) {
	std::string regEnd = std::string("\\]"); 
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

			bool linhaVazia = line == ""; 
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
		
	}

	mpData.version = std::stod(versionData);
	mpData.baseMVA = std::stod(baseMVAData);

	return mpData;
}

matPowerDataType lerMatPowerEAlocarMemoria(std::string mFile, sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo) {
	matPowerDataType mpData = lerArquivoMatPOWEReficiente(mFile);

	sistema.nB = mpData.bus.BUS_I.size();
	sistema.nL = mpData.branch.F_BUS.size();
	sistema.nPQ = std::count(mpData.bus.BUS_TYPE.begin(), mpData.bus.BUS_TYPE.end(), 1);
	sistema.nPV = std::count(mpData.bus.BUS_TYPE.begin(), mpData.bus.BUS_TYPE.end(), 2);

	initBranch(sistema, ramo);
	initSistema(sistema);
	initBus(sistema, barra);
	initIter(sistema, iterativo);

	if (storeMatPOWER(sistema, barra, ramo, mpData)) { 
		printf("Deu ruim...");
	}
	return mpData;
}