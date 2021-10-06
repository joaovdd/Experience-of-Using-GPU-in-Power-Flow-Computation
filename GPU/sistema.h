#pragma once
#include <cuComplex.h>

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

#include <vector>

#include "sistInfo.h"

#include <complex>
//#include "Eigen/Dense"
#include "Eigen/Sparse"

enum subMatJ {
	//nda = 0,
	H = 1,
	L = 2,
	M = 3,
	N = 4
};

struct sistema {
	int* barrasPV;//[NPV]; // indexador simples (iterador)
	int nPV;// = 0;
	int* barrasPQ;//[NPQ]; // indexador simples (iterador)
	int nPQ;// = 0;
	int barraVO;// = 0;
	int nB;
	int nL;
	float_type baseMVA;

	complex_type* Y;//[NB*NB]; // elementos da matriz admitância (via IDX2F)

	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>* spY;
	complex_type* spYval;
	std::complex<float_type>* spYvalE; // uso do eigen
	int* csrRowPtrY;
	int* csrColIndY;
	int* cooRowIndY;
	int nnzY;

	//limite de injeção de reativos
	float_type* limQinf; // [NPV] // limite inferior de injeção de reativos de cada barra PV
	float_type* limQsup; // [NPV] // limite superior de injeção de reativos de cada barra PV

	float_type* VfixadaPV; // [NPV] // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat
};

struct h_sparse {
	std::vector<float_type> spJval;
	/* int* */ std::vector<subMatJ> spJsubmatType;
	/* int* */ std::vector<int> cooColIndSubMatJ;
	/* int* */ std::vector<int> cooRowIndSubMatJ;
	/* int* */ std::vector<int> cooColIndJ;
	/* int* */ std::vector<int> cooRowIndJ;
	std::vector<int> csrRowPtrJ;
	int nnzJ;

	std::vector<int> Hpos;
	std::vector<int> Lpos;
	std::vector<int> Mpos;
	std::vector<int> Npos;
};

struct d_sparse {
	float_type* spJval;
	subMatJ* spJsubmatType;
	int* cooColIndSubMatJ;
	int* cooRowIndSubMatJ;
	int* cooColIndJ;
	int* cooRowIndJ;
	int* csrRowPtrJ;
	int nnzJ;

	int* Hpos;
	int* Lpos;
	int* Mpos;
	int* Npos;
	int nnzH;
	int nnzL;
	int nnzM;
	int nnzN;
};

// indexador é o numero da barra (via IDX1F)
struct barra {
	int* id; // bus identifier number

	float_type* V;//[NB];
	float_type* theta;//[NB];

	float_type* Pliq;//[NB];
	float_type* Qliq;//[NB];
	float_type* Pload;//[NB];
	float_type* Qload;//[NB];
	float_type* Pg;//[NB];
	float_type* Qg;//[NB];

	float_type* Vbase;//[NB];

	float_type* gsh;//[NB];
	float_type* bsh;//[NB];

	float_type* phi;
};

// indexador simples (iterador): No do ramo
struct ramo {
	complex_type* z;//[NL];
	float_type* bsh;//[NL];// = 0;
	complex_type* tap;//[NL];// = 0;

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>* eigen_phi;
	int* d_csrColIndPhi; // gpu only
	int* d_csrRowPtrPhi; // gpu only
	float_type* phiVal;               // gpu only
	int nnzPhi;          // gpu only

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	int* de;//[NL];
	int* para;//[NL];
};

struct iterativo {
	float_type* Pcalc;//[NB];
	float_type* Qcalc;//[NB];

	//float_type* J;//[(NB-1+NPQ)*(NB-1+NPQ)]; -> Jlim
	int iteracao, noMax;

	//float_type* deltaP;
	//float_type* deltaQ;
	//float_type* g; -> gLim

	// limite de injeção de reativos
	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	//limite de injeção de reativos
	bool* limQ; // FUTURO?:dois bits: o 1o é 0 se o limite atingido é o superior, 1 se o inferior. O segundo é
				// 0 se o limite de injeção de reativos não foi atingido, 1 caso contrário
				// nPV entradas
	int* barrasPVlim; //[NPV]; // vetor contendo as barras de tipo PV após a análise de limite de reativos
	int nPVlim; // número de barras de tipo PV após a análise de limite de reativos
	int* barrasPQlim; //[NPV]; // vetor contendo as barras de tipo PV após a análise de limite de reativos
	int nPQlim; // número de barras de tipo PV após a análise de limite de reativos

	float_type* QliqLim;       // atualização do limite do valor de Qliq para as novas barras de tipo PQ
	bool flgMudancaLimInjReat;
};

// popula estrutura csr com os dados de defasagem angular dos ramos do sistema e inicializa vetores
void InitCsrPhi(sistema& sistema, ramo& ramo) {

	std::vector<Eigen::Triplet<float_type>> valores;
	for (int i = 0; i < sistema.nL; i++) {
		valores.push_back(Eigen::Triplet<float_type>(IDX1F(ramo.de[i]), IDX1F(ramo.para[i]), ramo.phi[i]));
	}

	// ramo.eigen_phi = new Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>(sistema.nL, sistema.nL);

	ramo.eigen_phi = new Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>;

	ramo.eigen_phi->resize(sistema.nL, sistema.nL);

	ramo.eigen_phi->setFromTriplets(valores.begin(), valores.end());

	ramo.nnzPhi = ramo.eigen_phi->nonZeros();

	//ramo.csrColIndPhi = (int*)malloc(ramo.phi->nonZeros() * sizeof(int));
	//for (int j = 0; j < ramo.eigen_phi->nonZeros(); j++) {
	//	ramo.csrColIndPhi[j] = ramo.eigen_phi->innerIndexPtr()[j] /*+ 1*/;
	//}

	//ramo.csrRowPtrPhi = (int*)malloc((ramo.eigen_phi->outerSize() + 1) * sizeof(int));
	//for (int j = 0; j <= ramo.eigen_phi->outerSize(); j++) {
	//	ramo.csrRowPtrPhi[j] = ramo.eigen_phi->outerIndexPtr()[j] /*+ 1*/;
	//}
}

void initBranch(sistema& sistema, ramo& ramo) {
	ramo.z = (complex_type*)malloc(sistema.nL * sizeof(complex_type));
	ramo.bsh = (float_type*)malloc(sistema.nL * sizeof(float_type));
	ramo.tap = (complex_type*)malloc(sistema.nL * sizeof(complex_type));
	ramo.phi = (float_type*)malloc(sistema.nL * sizeof(float_type));

	ramo.de = (int*)malloc(sistema.nL * sizeof(float_type));
	ramo.para = (int*)malloc(sistema.nL * sizeof(float_type));

	ramo.Pdp = (float_type*)malloc(sistema.nL * sizeof(float_type));
	ramo.Qdp = (float_type*)malloc(sistema.nL * sizeof(float_type));
	ramo.Ppd = (float_type*)malloc(sistema.nL * sizeof(float_type));
	ramo.Qpd = (float_type*)malloc(sistema.nL * sizeof(float_type));


	for (int i = 0; i < sistema.nL; i++) {
		ramo.z[i] = _mkComplex(0., 0.);
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.bsh[i] = 0.;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.tap[i] = _mkComplex(1., 0.);
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.de[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.para[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Pdp[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Qdp[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Ppd[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Qpd[i] = 0;
	}
}

void initSistema(sistema& sistema) {
	//sistema.nPV = 0;
	//sistema.nPQ = 0;
	sistema.barraVO = 0;
	//sistema.nL = 0;
	//sistema.nB = 0;
	//sistema.baseMVA = 0;

	if (global::metodo == metodo::denso) {
		sistema.Y = (complex_type*)malloc(sistema.nB * sistema.nB * sizeof(complex_type));
		for (int i = 0; i < sistema.nB * sistema.nB; i++) {
			sistema.Y[i] = _mkComplex(0., 0.);
		}
	}
	else {
		sistema.Y = nullptr;
	}

	sistema.barrasPV = (int*)malloc(sistema.nPV * sizeof(int));
	sistema.barrasPQ = (int*)malloc(sistema.nPQ * sizeof(int));


	//limite de injeção de reativos
	sistema.limQinf = (float_type*)malloc(sistema.nPV * sizeof(float_type));
	sistema.limQsup = (float_type*)malloc(sistema.nPV * sizeof(float_type));

	sistema.VfixadaPV = (float_type*)malloc(sistema.nPV * sizeof(float_type));

	sistema.spY = nullptr;
	sistema.csrRowPtrY = nullptr;
	sistema.csrColIndY = nullptr;
}

void initBus(sistema& sistema, barra& barra) {
	barra.id = (int*)malloc(sistema.nB * sizeof(int));

	barra.V = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.theta = (float_type*)malloc(sistema.nB * sizeof(float_type));

	barra.Pliq = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.Qliq = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.Pload = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.Qload = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.Pg = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.Qg = (float_type*)malloc(sistema.nB * sizeof(float_type));

	barra.Vbase = (float_type*)malloc(sistema.nB * sizeof(float_type));

	barra.gsh = (float_type*)malloc(sistema.nB * sizeof(float_type));
	barra.bsh = (float_type*)malloc(sistema.nB * sizeof(float_type));

	barra.phi = (float_type*)malloc(sistema.nB * sistema.nB * sizeof(float_type));

	for (int i = 0; i < sistema.nB; i++) {
		barra.V[i] = global::v_inicial;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.theta[i] = global::theta_inicial;
	}

	for (int i = 0; i < sistema.nB; i++) {
		barra.Pliq[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.Qliq[i] = 0.;
	}

	for (int i = 0; i < sistema.nB; i++) {
		barra.Pload[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.Qload[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.Pg[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.Qg[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.Vbase[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.gsh[i] = 0.;
	}
	for (int i = 0; i < sistema.nB; i++) {
		barra.bsh[i] = 0.;
	}
	for (int i = 0; i < sistema.nB * sistema.nB; i++) {
		barra.phi[i] = 0.;
	}
}

void initIter(sistema& sistema, iterativo& iterativo) {
	iterativo.Pcalc = (float_type*)malloc(sistema.nB * sizeof(float_type));
	iterativo.Qcalc = (float_type*)malloc(sistema.nB * sizeof(float_type));

	//iterativo.g = NULL;
	//iterativo.g = (float_type*)malloc((sistema.nPV + sistema.nPQ + sistema.nPQ) * sizeof(float_type));

	//iterativo.J = (float_type*)malloc((sistema.nB - 1 + sistema.nPQ) * (sistema.nB - 1 + sistema.nPQ) * sizeof(float_type));

	iterativo.iteracao = 0;
	iterativo.noMax = global::no_max_iter;
	for (int i = 0; i < sistema.nB - 1; i++) {
		iterativo.Pcalc[i] = 0.;
	}
	for (int i = 0; i < sistema.nPQ; i++) {
		iterativo.Qcalc[i] = 0.;
	}
	//for (int i = 0; i < (sistema.nB - 1 + sistema.nPQ) * (sistema.nB - 1 + sistema.nPQ); i++) {
	//	iterativo.J[i] = 0.;
	//}


	// limite injeção de reativos
	iterativo.gLim = (float_type*)malloc((sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * sizeof(float_type));
	if (global::metodo == metodo::denso) {
		iterativo.Jlim = (float_type*)malloc((sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * sizeof(float_type));
		for (int i = 0; i < (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ) * (sistema.nPV + sistema.nPV + sistema.nPQ + sistema.nPQ); i++) {
			iterativo.Jlim[i] = 0.;
		}
	}
	else {
		iterativo.Jlim = nullptr;
	}

	//limite de injeção de reativos
	iterativo.limQ = (bool*)malloc(sistema.nPV * sizeof(bool));
	iterativo.barrasPVlim = (int*)malloc(sistema.nPV * sizeof(int));
	iterativo.barrasPQlim = (int*)malloc((sistema.nPV + sistema.nPQ) * sizeof(int));
	iterativo.QliqLim = (float_type*)malloc(sistema.nB * sizeof(float_type));
}

void finSistema(sistema& sistema) {
	if (sistema.Y != nullptr) { free(sistema.Y); }

	if (sistema.barrasPV != nullptr) { free(sistema.barrasPV); }
	if (sistema.barrasPQ != nullptr) { free(sistema.barrasPQ); }

	if (sistema.limQinf != nullptr) { free(sistema.limQinf); }
	if (sistema.limQsup != nullptr) { free(sistema.limQsup); }

	if (sistema.VfixadaPV == nullptr) { free(sistema.VfixadaPV); }

	if (sistema.spY != nullptr) { delete(sistema.spY); }
	// if (sistema.csrRowPtrY != nullptr) { free(sistema.csrRowPtrY); }
	// if (sistema.csrColIndY != nullptr) { free(sistema.csrColIndY); }
}

void finBranch(ramo& ramo) {
	free(ramo.z);
	free(ramo.bsh);
	free(ramo.tap);
	free(ramo.phi);
	free(ramo.Pdp);
	free(ramo.Qdp);
	free(ramo.Ppd);
	free(ramo.Qpd);

	free(ramo.de);
	free(ramo.para);

	if (ramo.eigen_phi != nullptr) { delete(ramo.eigen_phi); }

}

void finBus(barra& barra) {
	if (barra.id == nullptr) { free(barra.id); }
	free(barra.V);
	free(barra.theta);

	free(barra.Pliq);
	free(barra.Qliq);
	free(barra.Pload);
	free(barra.Qload);
	free(barra.Pg);
	free(barra.Qg);

	free(barra.Vbase);

	free(barra.gsh);
	free(barra.bsh);
}

void finIter(iterativo& iterativo) {
	free(iterativo.Pcalc);
	free(iterativo.Qcalc);

	// limite injeção de reativos
	free(iterativo.gLim);
	if (iterativo.Jlim != nullptr) { free(iterativo.Jlim); }
	free(iterativo.limQ);
	free(iterativo.barrasPVlim);
	free(iterativo.barrasPQlim);
	free(iterativo.QliqLim);
}

bool lerTamanhos(std::string cdfFile, sistema& sistema) {
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
				switch (atoi(line.substr(24, 2).c_str())) { // [24, 25]
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
		std::vector<int> de, para;
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {

				int auxde = atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
				int auxpara = atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
				bool flgRamoNovo = 1;

				for (int i = 0; i < sistema.nL; i++) {
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

//bool readCDF(std::string cdfFile, sistema& sistema, barra& barra, ramo& ramo) {
//	//std::vector <std::vector <branch>> vecBranch;
//
//	std::string line;
//	std::fstream CDF;
//	CDF.open(cdfFile.c_str(), std::ios::in);
//	if (CDF.is_open())
//	{
//		sistema.nB = 0;
//		sistema.nPQ = 0;
//		sistema.nPV = 0;
//		sistema.barraVO = 0;
//		std::getline(CDF, line);
//
//		//std::cout << "linha: " << atoi(line.substr(31, 6).c_str()) << std::endl;
//		sistema.baseMVA = atoi(line.substr(31, 6).c_str());
//
//		std::getline(CDF, line); // pula duas linhas
//		while (std::getline(CDF, line)) // atualiza line a cada itera��o
//		{
//			if (line.find("-999") == std::string::npos) {
//				//std::cout << "Val: " << line.substr(0, 4) << std::endl;
//				int i = atoi(line.substr(0, 4).c_str()); // [0, 3]
//
//				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]
//
//				switch (atoi(line.substr(24, 2).c_str())) { // [24, 25]
//				case 0: // tipo PQ [Unregulated (load, PQ)]
//					sistema.barrasPQ[sistema.nPQ] = i;
//					sistema.nPQ++;
//					break;
//				case 1: // tipo PQ [Hold MVAR generation within voltage limits, (PQ)]
//					sistema.barrasPQ[sistema.nPQ] = i;
//					sistema.nPQ++;
//					break;
//				case 2: // tipo PV [Hold voltage within VAR limits (gen, PV)]
//					sistema.barrasPV[sistema.nPV] = i;
//					sistema.nPV++;
//
//					barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
//
//					//Limites
//					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(90, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
//					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(98, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]
//
//					//if (barra.Vbase[IDX1F(i)] == 0.)
//					//	barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());  // Vbase é zero quando não está sendo usada...
//					//else
//					//	barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str())*barra.Vbase[IDX1F(i)];  // Final voltage, p.u. (F)
//					break;
//				case 3:// tipo VO [Hold voltage and angle (swing, V-Theta)]
//					if (sistema.barraVO == 0) {
//						sistema.barraVO = i;
//
//						barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
//
//						//if (barra.Vbase[IDX1F(i)] == 0.)
//						//	barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());  // Vbase é zero quando não está sendo usada...
//						//else
//						//	barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str())*barra.Vbase[IDX1F(i)];  // Final voltage, p.u. (F)
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
//				barra.Pload[IDX1F(i)] = atof(line.substr(40, 9).c_str()) / sistema.baseMVA;  // Load MW // [40, 48]
//				barra.Qload[IDX1F(i)] = atof(line.substr(49, 10).c_str()) / sistema.baseMVA; // Load MVAR // [49, 58]
//				barra.Pg[IDX1F(i)] = atof(line.substr(59, 8).c_str()) / sistema.baseMVA;  // Generation MW // [59, 66]
//				barra.Qg[IDX1F(i)] = atof(line.substr(67, 8).c_str()) / sistema.baseMVA;  // Generation MVAR // [67, 74]
//				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]
//				barra.bsh[IDX1F(i)] = atof(line.substr(114, 8).c_str()); // 
//				barra.gsh[IDX1F(i)] = atof(line.substr(106, 8).c_str()); // Generation MVAR // [67, 74]
//
//
//				barra.Pliq[IDX1F(i)] = barra.Pg[IDX1F(i)] - barra.Pload[IDX1F(i)];  // MW líquidos
//				barra.Qliq[IDX1F(i)] = barra.Qg[IDX1F(i)] - barra.Qload[IDX1F(i)];  // MVAR líquidos
//
//			}
//			else
//				break; // se -999 a se��o Bus Data acabou!
//		}
//
//		//barra.Pliq[IDX1F(sistema.barraVO)] = 0;
//		//barra.Qliq[IDX1F(sistema.barraVO)] = 0;
//
//
//		std::getline(CDF, line); // pula duas linhas
//		sistema.nL = 0;
//		int nRamosDuplicatas = 0;
//		while (std::getline(CDF, line)) // atualiza line a cada itera��o
//		{
//			if (line.find("-999") == std::string::npos) {
//				// Parse branch
//
//				//ramo.de[IDX1F(sistema.nL)]   = atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
//				//ramo.para[IDX1F(sistema.nL)] = atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
//
//				int auxde = atoi(line.substr(0, 4).c_str()); // ramo da linha i // [0, 3]
//				int auxpara = atoi(line.substr(4, 5).c_str()); // at� a j // [4, 8]
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
//					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(83, 7).c_str()); // Transformer (phase shifter) final angle (F)
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
//					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(83, 7).c_str());
//					complex_type aux;
//					if ((auxA != 0.)) {
//						aux = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
//					}
//					else {
//						/*ramo.tap[IDX1F(sistema.nL)]*/ aux = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
//					}
//					// se valor lido é diferente do já armazenado: atenção!
//					if ((aux.x != ramo.tap[i].x) || (aux.y != ramo.tap[i].y))
//						printf("Atenção! Os parâmetros lidos para os transformadores de potência das linhas de transmissão paralelas são diferentes! (linha %d)\n", i);
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
//		CDF.close();
//	}
//	else {
//		std::cout << "ERRO: arquivo .cdf não pode ser aberto!" << std::endl;
//		return 1;
//	}
//	// número de barras
//	sistema.nB = sistema.nPQ + sistema.nPV + 1;
//
//	//barra.Pg[IDX1F(sistema.barraVO)] = 0;
//	//barra.Qg[IDX1F(sistema.barraVO)] = 0;
//	//barra.Pload[IDX1F(sistema.barraVO)] = 0;
//	//barra.Qload[IDX1F(sistema.barraVO)] = 0;
//
//	//barra.Pliq[IDX1F(sistema.barraVO)] = 0;
//	//barra.Qliq[IDX1F(sistema.barraVO)] = 0;
//
//	//for (int i = 1; i < sistema.nL; i++) {
//	//	ramo.tap[IDX1F(i)] = _mkComplex(1., 0.);
//	//}
//
//	return 0;
//}

int id2i(int id, sistema& sistema, barra& barra) {
	// TODO: pesquisar por id em barra.id[] eretornar posição i
	for (int i = 0; i <= sistema.nB; i++) {
		if (barra.id[i] == id) {
			return i + 1;
		}
	}
	throw - 1;
}

bool readCDF(std::string cdfFile, sistema& sistema, barra& barra, ramo& ramo) {
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

		int i = 0; // i é o número da barra que está sendo lida

		std::getline(CDF, line); // pula duas linhas
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				i++; // nova barra encontrada

				barra.id[IDX1F(i)] = atoi(line.substr(0, 4).c_str()); // [0, 3]

				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]

				switch (atoi(line.substr(24, 2).c_str())) { // [24, 25]
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

					barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
					sistema.VfixadaPV[IDX1F(sistema.nPV)] = barra.V[IDX1F(i)]; // preserva VfixadaPV caso a barra se torne PQ por violacao de seus LimInjReat

					//Limites
					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(90, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [91, 98]
					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(98, 7).c_str()) / sistema.baseMVA;  // Maximum MVAR or voltage limit (F) // [99, 106]

					break;
				case 3:// tipo VO [Hold voltage and angle (swing, V-Theta)]
					if (sistema.barraVO == 0) {
						sistema.barraVO = i;

						barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());

						barra.theta[IDX1F(i)] = atof(line.substr(33, 7).c_str()) * 3.14159265358979323846 / 180.;  // Final angle, degrees
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
				barra.Pload[IDX1F(i)] = atof(line.substr(40, 9).c_str()) / sistema.baseMVA;  // Load MW // [40, 48]
				barra.Qload[IDX1F(i)] = atof(line.substr(49, 10).c_str()) / sistema.baseMVA; // Load MVAR // [49, 58]
				barra.Pg[IDX1F(i)] = atof(line.substr(59, 8).c_str()) / sistema.baseMVA;  // Generation MW // [59, 66]
				barra.Qg[IDX1F(i)] = atof(line.substr(67, 8).c_str()) / sistema.baseMVA;  // Generation MVAR // [67, 74]
				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  // Base KV // [76, 82]
				barra.bsh[IDX1F(i)] = atof(line.substr(114, 8).c_str()); // 
				barra.gsh[IDX1F(i)] = atof(line.substr(106, 8).c_str()); // Generation MVAR // [67, 74]

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
		int nRamosDuplicatas = 0;
		while (std::getline(CDF, line)) // atualiza line a cada itera��o
		{
			if (line.find("-999") == std::string::npos) {
				// Parse branch	
				int auxde = id2i(atoi(line.substr(0, 4).c_str()), sistema, barra); // ramo da linha i // [0, 3]
				int auxpara = id2i(atoi(line.substr(4, 5).c_str()), sistema, barra); // at� a j // [4, 8]
				bool flgRamoNovo = 1;
				int i;

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

					ramo.z[IDX1F(sistema.nL)] = _mkComplex(atof(line.substr(19, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
						atof(line.substr(29, 11).c_str())); // ramo reactance  X, per unit [29, 39]

					float_type auxA = atof(line.substr(76, 6).c_str()); // Transformer final turns ratio (F)
					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.; // Transformer (phase shifter) final angle (F)
					if (auxA != 0.) {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); // a * cis phi
					}
					else {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
					}

					ramo.bsh[IDX1F(sistema.nL)] = atof(line.substr(40, 10).c_str()); // Line charging B, per unit (F) * (total line charging, +B)
				}
				else { // ramo i é duplicata do nL lido
					if (!global::laconic_mode) {
						printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", ramo.de[IDX1F(sistema.nL)], ramo.para[IDX1F(sistema.nL)]);
					}

					// conferência do tap
					float_type auxA = atof(line.substr(76, 6).c_str()); // tap
					float_type auxPhi = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.;
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

					aux = _mkComplex(atof(line.substr(19, 10).c_str()),  // ramo resistance R, per unit [19, 28] ,
						atof(line.substr(29, 11).c_str())); // ramo reactance  X, per unit [29, 39]

					ramo.z[i] = _cuDiv(_cuMul(aux, ramo.z[i]), _cuAdd(ramo.z[i], aux)); // A//B = (A*B)/(A+B)

					// bsh
					float_type aux2 = atof(line.substr(40, 10).c_str());   // Line charging B, per unit (F) * (total line charging, +B)
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

// OK!
void calcYbus2(sistema& sistema, barra& barra, ramo& ramo) {
	//inicializar Y

	//percorre ramos
	for (int i = 1; i <= sistema.nL; i++) {
		complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
		aux = _cuMul(_mkComplex(-1., 0.), aux); // -t*
		complex_type aux1 = _cuDiv(aux, ramo.z[IDX1F(i)]);	// -t*/z
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux1; // linha única entre de[i] e para[i]

		aux = _cuMul(ramo.z[IDX1F(i)], aux);// -t* * z
		aux1 = _cuDiv(_mkComplex(1., 0.), aux); // -1/(t* * z)
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux1;

		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
		dAux = dAux * dAux; // |t_ft|^2

		aux = _mkComplex(dAux, 0.); // |t_ft|^2
		aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_ft|^2*y_ft
		aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

		aux = _mkComplex(1 / dAux, 0.); // |t_ft|^-2 = |t_tf|^2
		//aux = _mkComplex(0., 0.);
		aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_tf|^2*(y_ft) = |t_tf|^2*(y_tf)
		aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); // Y_tt + |t_tf|^2*y_tf
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + |t_tf|^2*y_tf + j*bsh_ft/2 = Y_tt + |t_tf|^2*y_tf + j*bsh_tf/2
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
	}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
	}
}

void calcYbus(sistema& sistema, barra& barra, ramo& ramo) {
	// percorre ramos
		for (int i = 1; i <= sistema.nL; i++) {
			complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
			aux = _cuDiv(_mkComplex(-1., 0.), aux); // -1/t*

			complex_type aux2 = _cuDiv(_mkComplex(-1., 0.), ramo.tap[IDX1F(i)]); // -1/t
			aux2 = _cuDiv(aux2, ramo.z[IDX1F(i)]); // -1/(t * z)

			sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = _cuDiv(aux, ramo.z[IDX1F(i)]); // linha única entre de[i] e para[i]
			sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux2;

			float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
			dAux = 1 / (dAux * dAux); // |t_ft|^2

			aux = _mkComplex(dAux, 0.); // 1/|t_ft|^2
			aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft/|t_ft|^2
			aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2/y_ft
			aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
			sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

			aux = _cuDiv(_mkComplex(1., 0.), ramo.z[IDX1F(i)]);
			aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); // Y_tt + y_tf
			aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + y_tf + j*bsh_ft/2 = Y_tt + y_tf + j*bsh_tf/2
			sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
		}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
	}
	
	//não matpower

	////percorre ramos
	//for (int i = 1; i <= sistema.nL; i++) {
	//	complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
	//	aux = _cuMul(_mkComplex(-1., 0.), aux); // -t*
	//	complex_type aux2 = _cuMul(ramo.z[IDX1F(i)], aux); // -t* * z
	//	aux2 = _cuDiv(_mkComplex(1., 0.), aux2); // -1/(t* * z)
	//	sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = _cuDiv(aux, ramo.z[IDX1F(i)]); // linha única entre de[i] e para[i]
	//	sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux2;

	//	float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
	//	dAux = (dAux * dAux); // |t_ft|^2

	//	aux = _mkComplex(dAux, 0.); // |t_ft|^2
	//	aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft*|t_ft|^2
	//	aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
	//	aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
	//	sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

	//	aux = _mkComplex(1 / dAux, 0.); // |t_ft|^-2 = |t_tf|^2
	//	aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_tf|^2*(y_ft) = |t_tf|^2*(y_tf)
	//	aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); // Y_tt + |t_tf|^2*y_tf
	//	aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + |t_tf|^2*y_tf + j*bsh_ft/2 = Y_tt + |t_tf|^2*y_tf + j*bsh_tf/2
	//	sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
	//}

	////percorre barras
	//for (int i = 1; i <= sistema.nB; i++) {
	//	sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
	//}





	////inicializar Y

	////percorre ramos
	//for (int i = 1; i <= sistema.nL; i++) {
	//	complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
	//	aux = _cuMul(_mkComplex(-1., 0.), aux); // -t*
	//	aux = _cuMul(ramo.z[IDX1F(i)], aux); // -t* * z
	//	aux = _cuDiv(_mkComplex(1., 0.), aux); // -1/(t* * z)
	//	sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux; // linha única entre de[i] e para[i]
	//	sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

	//	float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
	//	dAux = 1 / (dAux * dAux); // 1/|t_ft|^2

	//	aux = _mkComplex(dAux, 0.); // 1/|t_ft|^2
	//	aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft/|t_ft|^2
	//	aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
	//	aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
	//	sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

	//	aux = _mkComplex(1, 0.); // |t_ft|^-2 = |t_tf|^2
	//	aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_tf|^2*(y_ft) = |t_tf|^2*(y_tf)
	//	aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); // Y_tt + |t_tf|^2*y_tf
	//	aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + |t_tf|^2*y_tf + j*bsh_ft/2 = Y_tt + |t_tf|^2*y_tf + j*bsh_tf/2
	//	sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
	//}

	////percorre barras
	//for (int i = 1; i <= sistema.nB; i++) {
	//	sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
	//}
}

//void calcYbusSp(sistema& sistema, barra& barra, ramo& ramo) {
//	//inicializar Y
//
//	using namespace std::complex_literals;
//
//	Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);
//
//	//percorre ramos
//	for (int i = 1; i <= sistema.nL; i++) {
//
//		std::complex<float_type> aux = -1. / ((ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) * std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i));
//		dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux;
//		dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux;
//
//		// complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
//		// aux = _cuMul(_mkComplex(-1., 0.), aux); // -t*
//		// aux = _cuMul(ramo.z[IDX1F(i)], aux); // -t* * z
//		// aux = _cuDiv(_mkComplex(1., 0.), aux); // -1/(t* * z)
//		// // sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux; // linha única entre de[i] e para[i]
//		// // sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux.x + aux.y*1i;
//		// dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux.x + aux.y*1i;
//
//		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
//		dAux = 1 / (dAux * dAux); // 1/|t_ft|^2
//
//		dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
//
//		// aux = _mkComplex(dAux, 0.); // 1/|t_ft|^2
//		// aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft/|t_ft|^2
//		// aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
//		// aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
//		// // sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux.x + aux.y*1i;
//
//		dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = 1. / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
//
//		// aux = _mkComplex(1 , 0.); // |t_ft|^-2 = |t_tf|^2
//		// aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_tf|^2*(y_ft) = |t_tf|^2*(y_tf)
//		// aux = _cuAdd(/*sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)]*/dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])), aux); // Y_tt + |t_tf|^2*y_tf
//		// aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + |t_tf|^2*y_tf + j*bsh_ft/2 = Y_tt + |t_tf|^2*y_tf + j*bsh_tf/2
//		// // sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux.x + aux.y*1i;
//	}
//
//	//percorre barras
//	for (int i = 1; i <= sistema.nB; i++) {
//		// sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
//		dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)] * 1i;
//	}
//
//	///* Eigen::MatrixXcd */ Eigen::SparseMatrix<std::complex<float_type>> spY = dnY.sparseView(); // converte para esparsa;
//
//	//Eigen::SparseMatrix<std::complex<float_type>, RowMajor> spY = dnY.sparseView(); // converte para esparsa;
//
//	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>;
//	*sistema.spY = dnY.sparseView(); // converte para esparsa;
//
//
//	//spY.SparseMatrix::makeCompressed();
//
//	// sistema.spYval = new Eigen::SparseMatrix<float_type>;
//	// *(sistema.spYval) = spY;
//
//	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
//	sistema.nnzY = sistema.spY->nonZeros();
//	sistema.spYvalE = sistema.spY->valuePtr();
//	//sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
//	//sistema.csrColIndY = sistema.spY->outerIndexPtr();
//
//	// conversão de sistema.spYvalE de std::complex<float_type> para complex_type
//	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
//	for (int i = 0; i < sistema.nnzY; i++) {
//		sistema.spYval[i].x = real(sistema.spYvalE[i]);
//		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
//	}
//
//	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
//	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
//		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] + 1;
//		// std::cout << sistema.csrColIndY[j] << ' ';
//	}
//	// std::cout << std::endl;
//
//	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
//	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
//		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] + 1;
//		// std::cout << sistema.csrRowPtrY[j] << ' ';
//	}
//
//	sistema.cooRowIndY = (int*)malloc(dnY.nonZeros() * sizeof(int));
//	{
//		int acc = 0;
//		// csr to coo
//		for (int j = 0; j < sistema.spY->outerSize(); j++) {
//			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
//				sistema.cooRowIndY[acc] = j + 1;
//				// std::cout << j /*+ 1*/ << ' ';
//			}
//		}
//		//std::cout << std::endl;
//	}
//}
//
//void calcYbusSp_0based(sistema& sistema, barra& barra, ramo& ramo) {
//	//inicializar Y
//
//	using namespace std::complex_literals;
//
//	Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);
//
//	//percorre ramos
//	for (int i = 1; i <= sistema.nL; i++) {
//
//		std::complex<float_type> aux = -1. / ((ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) * std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i));
//		dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux;
//		dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux;
//
//		// complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); // t*
//		// aux = _cuMul(_mkComplex(-1., 0.), aux); // -t*
//		// aux = _cuMul(ramo.z[IDX1F(i)], aux); // -t* * z
//		// aux = _cuDiv(_mkComplex(1., 0.), aux); // -1/(t* * z)
//		// // sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux; // linha única entre de[i] e para[i]
//		// // sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux.x + aux.y*1i;
//		// dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux.x + aux.y*1i;
//
//		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
//		dAux = 1 / (dAux * dAux); // 1/|t_ft|^2
//
//		dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
//
//		// aux = _mkComplex(dAux, 0.); // 1/|t_ft|^2
//		// aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft/|t_ft|^2
//		// aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
//		// aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
//		// // sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux.x + aux.y*1i;
//
//		dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = 1. / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
//
//		// aux = _mkComplex(1 , 0.); // |t_ft|^-2 = |t_tf|^2
//		// aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // |t_tf|^2*(y_ft) = |t_tf|^2*(y_tf)
//		// aux = _cuAdd(/*sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)]*/dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])), aux); // Y_tt + |t_tf|^2*y_tf
//		// aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_tt + |t_tf|^2*y_tf + j*bsh_ft/2 = Y_tt + |t_tf|^2*y_tf + j*bsh_tf/2
//		// // sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux.x + aux.y*1i;
//	}
//
//	//percorre barras
//	for (int i = 1; i <= sistema.nB; i++) {
//		// sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
//		dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)] * 1i;
//	}
//
//	///* Eigen::MatrixXcd */ Eigen::SparseMatrix<std::complex<float_type>> spY = dnY.sparseView(); // converte para esparsa;
//
//	//Eigen::SparseMatrix<std::complex<float_type>, RowMajor> spY = dnY.sparseView(); // converte para esparsa;
//
//	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>;
//	*sistema.spY = dnY.sparseView(); // converte para esparsa;
//
//
//	//spY.SparseMatrix::makeCompressed();
//
//	// sistema.spYval = new Eigen::SparseMatrix<float_type>;
//	// *(sistema.spYval) = spY;
//
//	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
//	sistema.nnzY = sistema.spY->nonZeros();
//	sistema.spYvalE = sistema.spY->valuePtr();
//	//sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
//	//sistema.csrColIndY = sistema.spY->outerIndexPtr();
//
//	// conversão de sistema.spYvalE de std::complex<float_type> para complex_type
//	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
//	for (int i = 0; i < sistema.nnzY; i++) {
//		sistema.spYval[i].x = real(sistema.spYvalE[i]);
//		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
//	}
//
//	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
//	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
//		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j]/* + 1*/;
//		// std::cout << sistema.csrColIndY[j] << ' ';
//	}
//	// std::cout << std::endl;
//
//	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
//	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
//		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j]/* + 1*/;
//		// std::cout << sistema.csrRowPtrY[j] << ' ';
//	}
//
//	sistema.cooRowIndY = (int*)malloc(dnY.nonZeros() * sizeof(int));
//	{
//		int acc = 0;
//		// csr to coo
//		for (int j = 0; j < sistema.spY->outerSize(); j++) {
//			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
//				sistema.cooRowIndY[acc] = j/* + 1*/;
//				// std::cout << j /*+ 1*/ << ' ';
//			}
//		}
//		//std::cout << std::endl;
//	}
//}

void calcYbusSp_eficinte(sistema& sistema, barra& barra, ramo& ramo) {
	// not like matpower

	//using namespace std::complex_literals;

	//std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	////Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);

	////percorre ramos
	//for (int i = 1; i <= sistema.nL; i++) {

	//	std::complex<float_type> aux = -1. / ((ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i));
	//	// std::complex<float_type> aux = -1. / ((ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) * std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i));
	//	//dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = aux;
	//	//dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux;
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux * std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i)));
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), aux * 1. / std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i)));


	//	float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
	//	dAux = (dAux * dAux); // |t_ft|^2

	//	//dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + (ramo.bsh[IDX1F(i)] * 1i) / 2.));

	//	dAux = 1 / dAux; // 1/|t_ft|^2 = |t_tf|^2 

	//	//dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = 1. / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + (ramo.bsh[IDX1F(i)] * 1i) / 2.));
	//}

	////percorre barras
	//for (int i = 1; i <= sistema.nB; i++) {
	//	//dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)] * 1i;
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), barra.bsh[IDX1F(i)] * 1i));
	//}
	
	
	
	
	
	////inicializar Y

	//using namespace std::complex_literals;

	//std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	////percorre ramos
	//for (int i = 1; i <= sistema.nL; i++) {

	//	std::complex<float_type> aux = -1. / ((ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) * std::conj(ramo.tap[IDX1F(i)].x + ramo.tap[IDX1F(i)].y * 1i));
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux));
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), aux));


	//	float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); // |t_ft|
	//	dAux = 1 / (dAux * dAux); // 1/|t_ft|^2

	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + (ramo.bsh[IDX1F(i)] * 1i) / 2.));

	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), 1. / (ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y * 1i) + (ramo.bsh[IDX1F(i)] * 1i) / 2.));
	//}

	////percorre barras
	//for (int i = 1; i <= sistema.nB; i++) {
	//	valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), barra.bsh[IDX1F(i)] * 1i));
	//}

	//using namespace std::complex_literals;

	std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	//Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);

	//percorre ramos
	for (int i = 1; i <= sistema.nL; i++) {

		std::complex<float_type> aux = std::complex <float_type>(1., 0.) / (std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), -aux / std::conj(std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y))));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), -aux / std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y)));


		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]);
		dAux = 1. / (dAux * dAux);

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux * aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));
	}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		//dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)] * 1i;
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])));
	}

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);
	sistema.spY->setFromTriplets(valores.begin(), valores.end());

	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();

	// conversão de sistema.spYvalE de std::complex<float_type> para complex_type
	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	for (int i = 0; i < sistema.nnzY; i++) {
		sistema.spYval[i].x = real(sistema.spYvalE[i]);
		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	}

	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] /*+ 1*/;
	}

	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] /*+ 1*/;
	}

	sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	{
		int acc = 0;
		// csr to coo
		for (int j = 0; j < sistema.spY->outerSize(); j++) {
			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
				sistema.cooRowIndY[acc] = j /*+ 1*/;
			}
		}
	}
}

void calcYbusSp_Matpower(sistema& sistema, barra& barra, ramo& ramo) {
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplCt, triplCf;
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplYff, triplYft, triplYtf, triplYtt, triplYsh;

	triplCt.reserve(sistema.nL);
	triplCf.reserve(sistema.nL);
	triplYff.reserve(sistema.nL);
	triplYft.reserve(sistema.nL);
	triplYtf.reserve(sistema.nL);
	triplYtt.reserve(sistema.nL);
	triplYsh.reserve(sistema.nB);

	//percorre ramos
	//using namespace std::complex_literals;
	for (int i = 1; i <= sistema.nL; i++) {
		// construção das matrizes de conexao
		triplCf.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.de[IDX1F(i)]), 1));
		triplCt.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.para[IDX1F(i)]), 1));

		// construcao das matrizes diagonais
		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		auto daux = _cuAbs(ramo.tap[IDX1F(i)]); // cuCmul(_cuCon(ramo.tap[IDX1F(i)]), ramo.tap[IDX1F(i)]);
		aux = _cuMul(aux, _mkComplex(1/(daux*daux), 0));
		triplYff.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), _cuCon(ramo.tap[IDX1F(i)])), ramo.z[IDX1F(i)]);
		triplYft.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), ramo.tap[IDX1F(i)]), ramo.z[IDX1F(i)]);
		triplYtf.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		triplYtt.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));
	}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		auto aux = _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]);
		triplYsh.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));
	}

	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Cf(sistema.nL, sistema.nB), Ct(sistema.nL, sistema.nB);
	Cf.setFromTriplets(triplCf.begin(), triplCf.end());
	Ct.setFromTriplets(triplCt.begin(), triplCt.end());
	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yff(sistema.nL, sistema.nL), Ytt(sistema.nL, sistema.nL),
		Yft(sistema.nL, sistema.nL), Ytf(sistema.nL, sistema.nL), Ysh(sistema.nB, sistema.nB);
	Yff.setFromTriplets(triplYff.begin(), triplYff.end());
	Ytt.setFromTriplets(triplYtt.begin(), triplYtt.end());
	Yft.setFromTriplets(triplYft.begin(), triplYft.end());
	Ytf.setFromTriplets(triplYtf.begin(), triplYtf.end());
	Ysh.setFromTriplets(triplYsh.begin(), triplYsh.end());

	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yf = Yff * Cf + Yft * Ct;
	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yt = Ytf * Cf + Ytt * Ct;

	auto aux = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Cf.transpose());
	auto aux1 = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Ct.transpose());
	auto aux2 = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>((aux * Yf) + (aux1 * Yt));
	
	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> ans = aux2 + Ysh;

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

	*(sistema.spY) = ans;

	// Cópia da outra função

	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();
	//sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
	//sistema.csrColIndY = sistema.spY->outerIndexPtr();

    // conversão de sistema.spYvalE de std::complex<float_type> para complex_type
	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	for (int i = 0; i < sistema.nnzY; i++) {
		sistema.spYval[i].x = real(sistema.spYvalE[i]);
		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	}

	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j];
	}

	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j];
	}

	sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	{
		int acc = 0;
		// csr to coo
		for (int j = 0; j < sistema.spY->outerSize(); j++) {
			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
				sistema.cooRowIndY[acc] = j;
			}
		}
	}
}

#define MKL_INT int
#define MKL_Complex16 std::complex<float_type>

#include "mkl_types.h"
#include "mkl_spblas.h"

void printMKL_z_csr(const sparse_matrix_t *source) {
    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rows_start, *rows_end, *col_indx;
    MKL_Complex16 *values;

    sparse_status_t status = mkl_sparse_z_export_csr(*source, &indexing, &rows, &cols, &rows_start,
        &rows_end, &col_indx, &values);

    if(status == SPARSE_STATUS_SUCCESS) {
        printf("nRows: %d\nnCols: %d\nnnz: %d\n", rows, cols, rows_end[rows - 1] - indexing);
        printf("Row Pointer Vector:\n");
        for(int i = 0; i < rows - 1; i++) {
            printf("%d, ", rows_start[i]);
        }
        printf("%d\n", rows_start[rows - 1]);
		printf("Column Index Vector:\n");
        for(int i = 0; i < rows_end[rows - 1] - indexing - 1; i++){
            printf("%d, ", col_indx[i]);
        }
        printf("%d\n", col_indx[rows_end[rows - 1] - indexing - 1]);
		printf("Values Vector:\n");
        for(int i = 0; i < rows_end[rows - 1] - indexing - 1; i++){
            printf(" %f + j * %f\n", std::real(values[i]), std::imag(values[i]));
        }
        printf(" %f + j * %f\n", std::real(values[rows_end[rows - 1] - indexing - 1]), std::imag(values[rows_end[rows - 1] - indexing - 1]));
    }
    else {
        printf("ERRO\n");
    }
}

void calcYbusSp_Matpower_(sistema& sistema, barra& barra, ramo& ramo) {
	// std::vector<Eigen::Triplet<std::complex<float_type>>> triplCt, triplCf;
	// std::vector<Eigen::Triplet<std::complex<float_type>>> triplYff, triplYft, triplYtf, triplYtt, triplYsh;

    // struct matrix_descr descr;
	// descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    // descr.mode = SPARSE_FILL_MODE_UPPER;
    // descr.diag = SPARSE_DIAG_NON_UNIT;

	// sparse_matrix_t cooCt, cooCf, cooYff, cooYft, cooYtf, cooYtt, cooYsh;
	std::vector<MKL_INT>                  cooCt_rows, cooCf_rows, cooYff_rows, cooYft_rows, cooYtf_rows, cooYtt_rows, cooYsh_rows,
	                                      cooCt_cols, cooCf_cols, cooYff_cols, cooYft_cols, cooYtf_cols, cooYtt_cols, cooYsh_cols;
	std::vector<std::complex<float_type>> cooCt_vals, cooCf_vals, cooYff_vals, cooYft_vals, cooYtf_vals, cooYtt_vals, cooYsh_vals;
	// mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, n, n, ia, ia + 1, ja, a);

	// mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, descr, matrix_vectors[0], SPARSE_OPERATION_NON_TRANSPOSE, general_type, matrix_vectors[1], SPARSE_STAGE_NNZ_COUNT, &result_vectors[0]);

	//percorre ramos
	//using namespace std::complex_literals;
	for (int i = 1; i <= sistema.nL; i++) {
		// construção das matrizes de conexao
		cooCf_rows.push_back(IDX1F(i));
		cooCf_cols.push_back(IDX1F(ramo.de[IDX1F(i)]));
		cooCf_vals.push_back(1);

		cooCt_rows.push_back(IDX1F(i));
		cooCt_cols.push_back(IDX1F(ramo.para[IDX1F(i)]));
		cooCt_vals.push_back(1);

		// construcao das matrizes diagonais
		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));

		// aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		cooYtt_rows.push_back(IDX1F(i));
		cooYtt_cols.push_back(IDX1F(i));
		cooYtt_vals.emplace_back(std::complex<float_type>(_cuReal(aux), _cuImag(aux)));	

		auto daux = _cuAbs(ramo.tap[IDX1F(i)]);
		aux = _cuMul(aux, _mkComplex(1/(daux*daux), 0));
		cooYff_rows.push_back(IDX1F(i));
		cooYff_cols.push_back(IDX1F(i));
		cooYff_vals.emplace_back(std::complex<float_type>(_cuReal(aux), _cuImag(aux)));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), _cuCon(ramo.tap[IDX1F(i)])), ramo.z[IDX1F(i)]);
		cooYft_rows.push_back(IDX1F(i));
		cooYft_cols.push_back(IDX1F(i));
		cooYft_vals.emplace_back(std::complex<float_type>(_cuReal(aux), _cuImag(aux)));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), ramo.tap[IDX1F(i)]), ramo.z[IDX1F(i)]);
		cooYtf_rows.push_back(IDX1F(i));
		cooYtf_cols.push_back(IDX1F(i));
		cooYtf_vals.emplace_back(std::complex<float_type>(_cuReal(aux), _cuImag(aux)));		
	}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		// auto aux = _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]);
		// triplYsh.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
		// 	std::complex<float_type>(_cuReal(aux), _cuImag(aux))));
		cooYsh_rows.push_back(IDX1F(i));
		cooYsh_cols.push_back(IDX1F(i));
		cooYsh_vals.emplace_back(std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]));
	}

	// Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Cf(sistema.nL, sistema.nB), Ct(sistema.nL, sistema.nB);
	// Cf.setFromTriplets(triplCf.begin(), triplCf.end());
	// Ct.setFromTriplets(triplCt.begin(), triplCt.end());
	// Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yff(sistema.nL, sistema.nL), Ytt(sistema.nL, sistema.nL), Yft(sistema.nL, sistema.nL), Ytf(sistema.nL, sistema.nL), Ysh(sistema.nB, sistema.nB);
	// Yff.setFromTriplets(triplYff.begin(), triplYff.end());
	// Ytt.setFromTriplets(triplYtt.begin(), triplYtt.end());
	// Yft.setFromTriplets(triplYft.begin(), triplYft.end());
	// Ytf.setFromTriplets(triplYtf.begin(), triplYtf.end());
	// Ysh.setFromTriplets(triplYsh.begin(), triplYsh.end());

    // sparse_status_t mkl_sparse_d_create_csr( sparse_matrix_t           *A,
    //                                          const sparse_index_base_t indexing, /* indexing: C-style or Fortran-style */
    //                                          const MKL_INT             rows,
    //                                          const MKL_INT             cols,
    //                                                MKL_INT             *rows_start,
    //                                                MKL_INT             *rows_end,
    //                                                MKL_INT             *col_indx,
    //                                                double              *values );

	// mkl_sparse_z_create_coo(A, indexing, rows, cols, nnz, row_indx, col_indx, values )

	sparse_status_t status;
	sparse_matrix_t cooCt, cooCf, cooYff, cooYft, cooYtf, cooYtt, cooYsh, cooYf, cooYt, aux1, aux2, aux3;
	status = mkl_sparse_z_create_coo(&cooCt, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nB, cooCt_rows.size(), cooCt_rows.data(), cooCt_cols.data(), cooCt_vals.data());
	status = mkl_sparse_z_create_coo(&cooCf, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nB, cooCf_rows.size(), cooCf_rows.data(), cooCf_cols.data(), cooCf_vals.data());
	status = mkl_sparse_z_create_coo(&cooYff, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYff_rows.size(), cooYff_rows.data(), cooYff_cols.data(), cooYff_vals.data());
	status = mkl_sparse_z_create_coo(&cooYft, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYft_rows.size(), cooYft_rows.data(), cooYft_cols.data(), cooYft_vals.data());
	status = mkl_sparse_z_create_coo(&cooYtf, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYtf_rows.size(), cooYtf_rows.data(), cooYtf_cols.data(), cooYtf_vals.data());
	status = mkl_sparse_z_create_coo(&cooYtt, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYtt_rows.size(), cooYtt_rows.data(), cooYtt_cols.data(), cooYtt_vals.data());
	status = mkl_sparse_z_create_coo(&cooYsh, SPARSE_INDEX_BASE_ZERO, sistema.nB, sistema.nB, cooYsh_rows.size(), cooYsh_rows.data(), cooYsh_cols.data(), cooYsh_vals.data());
	// mkl_sparse_d_create_csr(&csrA, SPARSE_INDEX_BASE_ZERO, sistema.nB, sistema.nB, nnz, rowIndex, colIndex, values);

	status = mkl_sparse_convert_csr(cooCt, SPARSE_OPERATION_NON_TRANSPOSE, &cooCt);
	status = mkl_sparse_convert_csr(cooCf, SPARSE_OPERATION_NON_TRANSPOSE, &cooCf);
	status = mkl_sparse_convert_csr(cooYff, SPARSE_OPERATION_NON_TRANSPOSE, &cooYff);
	status = mkl_sparse_convert_csr(cooYft, SPARSE_OPERATION_NON_TRANSPOSE, &cooYft);
	status = mkl_sparse_convert_csr(cooYtf, SPARSE_OPERATION_NON_TRANSPOSE, &cooYtf);
	status = mkl_sparse_convert_csr(cooYtt, SPARSE_OPERATION_NON_TRANSPOSE, &cooYtt);
	status = mkl_sparse_convert_csr(cooYsh, SPARSE_OPERATION_NON_TRANSPOSE, &cooYsh);
	
	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("csrCt:\n");
		printMKL_z_csr(&cooCt);
		printf("csrCf:\n");
		printMKL_z_csr(&cooCf);
		printf("csrYff:\n");
		printMKL_z_csr(&cooYff);
		printf("csrYft:\n");
		printMKL_z_csr(&cooYft);
		printf("csrYtf:\n");
		printMKL_z_csr(&cooYtf);
		printf("csrYtt:\n");
		printMKL_z_csr(&cooYtt);
		printf("csrYsh:\n");
		printMKL_z_csr(&cooYsh);
	}

	MKL_Complex16 alpha = std::complex<float_type>(1,0);

	// Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yf = Yff * Cf + Yft * Ct;
	// Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> Yt = Ytf * Cf + Ytt * Ct;

	// Yf = Yff * Cf + Yft * Ct;
	// sparse_status_t mkl_sparse_sp2m ( const sparse_operation_t transA , const struct matrix_descr descrA , const sparse_matrix_t A , const sparse_operation_t transB , const struct matrix_descr descrB , const sparse_matrix_t B , const sparse_request_t request , sparse_matrix_t *C );
	// status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYff, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooCf, SPARSE_STAGE_FULL_MULT_NO_VAL, &aux1);
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYff, cooCf, &aux1);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("aux1:\n");
		printMKL_z_csr(&aux1);
	}

	// status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYft, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooCt, SPARSE_STAGE_NNZ_COUNT, &aux2);
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYft, cooCt, &aux2);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("aux2:\n");
		printMKL_z_csr(&aux2);
	}

	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux1, alpha, aux2, &cooYf);
	
	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("cooYf:\n");
		printMKL_z_csr(&cooYf);
	}

	mkl_sparse_destroy(aux1);
	mkl_sparse_destroy(aux2);
	// Yt = Ytf * Cf + Ytt * Ct;
	// status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYtf, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooCf, SPARSE_STAGE_NNZ_COUNT, &aux1);
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYtf, cooCf, &aux1);
	// status = mkl_sparse_sp2m(SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYtt, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooCt, SPARSE_STAGE_NNZ_COUNT, &aux2);
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYtt, cooCt, &aux2);
	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux1, alpha, aux2, &cooYt);
	mkl_sparse_destroy(aux1);
	mkl_sparse_destroy(aux2);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("cooYt:\n");
		printMKL_z_csr(&cooYt);
	}
	// auto aux = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Cf.transpose());
	// auto aux1 = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Ct.transpose());
	// auto aux2 = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>((aux * Yf) + (aux1 * Yt));
	
	// aux2 = (aux * Yf) + (aux1 * Yt)
	// status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, descr, cooCf, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYf, SPARSE_STAGE_NNZ_COUNT, &aux1);
	status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, cooCf, cooYf, &aux1);
	// status = mkl_sparse_sp2m(SPARSE_OPERATION_TRANSPOSE, descr, cooCt, SPARSE_OPERATION_NON_TRANSPOSE, descr, cooYt, SPARSE_STAGE_NNZ_COUNT, &aux2);
	status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, cooCt, cooYt, &aux2);
	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux1, alpha, aux2, &aux3);
	mkl_sparse_destroy(aux1);
	mkl_sparse_destroy(aux2);

	// Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> ans = aux2 + Ysh;
	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux3, alpha, cooYsh, &aux1);

	status = mkl_sparse_order(aux1);
	
	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("spY:\n");
		printMKL_z_csr(&aux1);
	}
	
	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

	// *(sistema.spY) = ans;

	/////////////////////////////////////////////////////////////////
	// Leitura de Ybus //////////////////////////////////////////////

    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rows_start, *rows_end, *col_indx;
    MKL_Complex16 *values;

    status = mkl_sparse_z_export_csr(aux1, &indexing, &rows, &cols, &rows_start,
        &rows_end, &col_indx, &values);

    if(status != SPARSE_STATUS_SUCCESS) {
		printf("Error!\n");
	}

	// int outerIndexPtr[cols+1];
	// int innerIndices[nnz];
	// double values[nnz];
	// sistema.spY =
	auto x = new Eigen::Map< const Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> >(
		rows, cols, rows_end[rows - 1] - indexing, rows_start, col_indx, values
	);
	
	// sistema.spY = x.sum(Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB));

	// for (int i = 0; i < rows; i++)
	// {
	// 	for (int j = rows_start[i]; j < rows_start[i+1]; j++)
	// 	{
	// 		i, col_indx[j], values[j]
	// 	}
	// }
	
	/////////////////////////////////////////////////////////////////

	// Cópia da outra função

	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	sistema.nnzY = rows_end[rows - 1] - indexing;
	sistema.spYvalE = values;

    // conversão de sistema.spYvalE de std::complex<float_type> para complex_type
	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	for (int i = 0; i < sistema.nnzY; i++) {
		sistema.spYval[i].x = real(sistema.spYvalE[i]);
		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	}

	sistema.csrColIndY = col_indx; // (int*)malloc(sistema.spY->nonZeros() * sizeof(int));

	sistema.csrRowPtrY = rows_start;// (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));

	// funcional /////////////////////////////////////////////////////////////////////////////////////////

	// // cria matriz

    // std::vector<Eigen::Triplet<std::complex<float_type>>> triplY;
	// triplY.reserve(rows_end[rows - 1] - indexing);

	// for (int i = 0; i < rows; i++)
	// {
	// 	for (int j = rows_start[i]; j < rows_start[i+1]; j++)
	// 	{
	// 		triplY.emplace_back(Eigen::Triplet<std::complex<float_type>>(i, col_indx[j],
	// 			std::complex<float_type>(std::real(values[j]), std::imag(values[j]))));
	// 	}
	// }
	
	// sistema.spY->setFromTriplets(triplY.begin(), triplY.end());

	// // Cópia da outra função 2

	// // https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	// sistema.nnzY = sistema.spY->nonZeros();
	// sistema.spYvalE = sistema.spY->valuePtr();
	// //sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
	// //sistema.csrColIndY = sistema.spY->outerIndexPtr();

    // // conversão de sistema.spYvalE de std::complex<float_type> para complex_type
	// sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	// for (int i = 0; i < sistema.nnzY; i++) {
	// 	sistema.spYval[i].x = real(sistema.spYvalE[i]);
	// 	sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	// }

	// sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	// for (int j = 0; j < sistema.spY->nonZeros(); j++) {
	// 	sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j];
	// }

	// sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	// for (int j = 0; j <= sistema.spY->outerSize(); j++) {
	// 	sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j];
	// }

	// sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	// {
	// 	int acc = 0;
	// 	// csr to coo
	// 	for (int j = 0; j < sistema.spY->outerSize(); j++) {
	// 		for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
	// 			sistema.cooRowIndY[acc] = j;
	// 		}
	// 	}
	// }

	// mkl_sparse_destroy(aux1);
	//          /////////////////////////////////////////////////////////////////////////////////////////

	mkl_sparse_destroy(aux3);
}