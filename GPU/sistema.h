#pragma once
#include <cuComplex.h>

#include <string>
#include <fstream>
#include <iostream>
#include <math.h>

#include <vector>

#include "sistInfo.h"

#include <complex>

#include "Eigen/Sparse"

enum subMatJ {
	H = 1,
	L = 2,
	M = 3,
	N = 4
};

struct sistema {
	int* barrasPV;
	int nPV;
	int* barrasPQ;
	int nPQ;
	int barraVO;
	int nB;
	int nL;
	float_type baseMVA;

	complex_type* Y;

	Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>* spY;
	complex_type* spYval;
	std::complex<float_type>* spYvalE; 
	int* csrRowPtrY;
	int* csrColIndY;
	int* cooRowIndY;
	int nnzY;

	float_type* limQinf; 
	float_type* limQsup; 

	float_type* VfixadaPV; 
};

struct h_sparse {
	std::vector<float_type> spJval;
	 std::vector<subMatJ> spJsubmatType;
	 std::vector<int> cooColIndSubMatJ;
	 std::vector<int> cooRowIndSubMatJ;
	 std::vector<int> cooColIndJ;
	 std::vector<int> cooRowIndJ;
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

struct barra {
	int* id; 

	float_type* V;
	float_type* theta;

	float_type* Pliq;
	float_type* Qliq;
	float_type* Pload;
	float_type* Qload;
	float_type* Pg;
	float_type* Qg;

	float_type* Vbase;

	float_type* gsh;
	float_type* bsh;

	float_type* phi;
};

struct ramo {
	complex_type* z;
	float_type* bsh;
	complex_type* tap;

	float_type* phi;
	Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>* eigen_phi;
	int* d_csrColIndPhi; 
	int* d_csrRowPtrPhi; 
	float_type* phiVal;               
	int nnzPhi;          

	float_type* Pdp;
	float_type* Ppd;
	float_type* Qdp;
	float_type* Qpd;

	int* de;
	int* para;
};

struct iterativo {
	float_type* Pcalc;
	float_type* Qcalc;

	int iteracao, noMax;

	float_type* gLim;
	int ngLim;
	float_type* Jlim;

	bool* limQ; 
				
	int* barrasPVlim; 
	int nPVlim; 
	int* barrasPQlim; 
	int nPQlim; 

	float_type* QliqLim;       
	bool flgMudancaLimInjReat;
};

void InitCsrPhi(sistema& sistema, ramo& ramo) {
	std::vector<Eigen::Triplet<float_type>> valores;
	for (int i = 0; i < sistema.nL; i++) {
		valores.push_back(Eigen::Triplet<float_type>(IDX1F(ramo.de[i]), IDX1F(ramo.para[i]), ramo.phi[i]));
	}

	ramo.eigen_phi = new Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>;

	ramo.eigen_phi->resize(sistema.nL, sistema.nL);

	ramo.eigen_phi->setFromTriplets(valores.begin(), valores.end());

	ramo.nnzPhi = ramo.eigen_phi->nonZeros();

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
	sistema.barraVO = 0;

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

	iterativo.iteracao = 0;
	iterativo.noMax = global::no_max_iter;
	for (int i = 0; i < sistema.nB - 1; i++) {
		iterativo.Pcalc[i] = 0.;
	}
	for (int i = 0; i < sistema.nPQ; i++) {
		iterativo.Qcalc[i] = 0.;
	}

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

		std::getline(CDF, line); std::getline(CDF, line); 
		while (std::getline(CDF, line)) 
		{
			if (line.find("-999") == std::string::npos) {
				sistema.nB++;	
				switch (atoi(line.substr(24, 2).c_str())) { 
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
				int auxde = atoi(line.substr(0, 4).c_str()); 
				int auxpara = atoi(line.substr(4, 5).c_str()); 
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

int id2i(int id, sistema& sistema, barra& barra) {
	for (int i = 0; i <= sistema.nB; i++) {
		if (barra.id[i] == id) {
			return i + 1;
		}
	}
	throw - 1;
}

bool readCDF(std::string cdfFile, sistema& sistema, barra& barra, ramo& ramo) {
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

				barra.id[IDX1F(i)] = atoi(line.substr(0, 4).c_str()); 

				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  

				switch (atoi(line.substr(24, 2).c_str())) { 
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

					barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());
					sistema.VfixadaPV[IDX1F(sistema.nPV)] = barra.V[IDX1F(i)]; 

					sistema.limQsup[IDX1F(sistema.nPV)] = atof(line.substr(90, 7).c_str()) / sistema.baseMVA;  
					sistema.limQinf[IDX1F(sistema.nPV)] = atof(line.substr(98, 7).c_str()) / sistema.baseMVA;  

					break;
				case 3:
					if (sistema.barraVO == 0) {
						sistema.barraVO = i;

						barra.V[IDX1F(i)] = atof(line.substr(27, 6).c_str());

						barra.theta[IDX1F(i)] = atof(line.substr(33, 7).c_str()) * 3.14159265358979323846 / 180.;  
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
				
				barra.Pload[IDX1F(i)] = atof(line.substr(40, 9).c_str()) / sistema.baseMVA;  
				barra.Qload[IDX1F(i)] = atof(line.substr(49, 10).c_str()) / sistema.baseMVA; 
				barra.Pg[IDX1F(i)] = atof(line.substr(59, 8).c_str()) / sistema.baseMVA;  
				barra.Qg[IDX1F(i)] = atof(line.substr(67, 8).c_str()) / sistema.baseMVA;  
				barra.Vbase[IDX1F(i)] = atof(line.substr(76, 7).c_str());  
				barra.bsh[IDX1F(i)] = atof(line.substr(114, 8).c_str()); 
				barra.gsh[IDX1F(i)] = atof(line.substr(106, 8).c_str()); 

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
				
				int auxde = id2i(atoi(line.substr(0, 4).c_str()), sistema, barra); 
				int auxpara = id2i(atoi(line.substr(4, 5).c_str()), sistema, barra); 
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

					ramo.z[IDX1F(sistema.nL)] = _mkComplex(atof(line.substr(19, 10).c_str()),  
						atof(line.substr(29, 11).c_str())); 

					float_type auxA = atof(line.substr(76, 6).c_str()); 
					ramo.phi[IDX1F(sistema.nL)] = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.; 
					if (auxA != 0.) {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(auxA * cos(ramo.phi[IDX1F(sistema.nL)]), auxA * sin(ramo.phi[IDX1F(sistema.nL)])); 
					}
					else {
						ramo.tap[IDX1F(sistema.nL)] = _mkComplex(cos(ramo.phi[IDX1F(sistema.nL)]), sin(ramo.phi[IDX1F(sistema.nL)]));
					}

					ramo.bsh[IDX1F(sistema.nL)] = atof(line.substr(40, 10).c_str()); 
				}
				else { 
					if (!global::laconic_mode) {
						printf("\nRamo %d-%d declarado duas vezes. Criando equivalente paralelo...", ramo.de[IDX1F(sistema.nL)], ramo.para[IDX1F(sistema.nL)]);
					}

					float_type auxA = atof(line.substr(76, 6).c_str()); 
					float_type auxPhi = atof(line.substr(83, 7).c_str()) * 3.14159265358979323846 / 180.;
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

					aux = _mkComplex(atof(line.substr(19, 10).c_str()),  
						atof(line.substr(29, 11).c_str())); 

					ramo.z[i] = _cuDiv(_cuMul(aux, ramo.z[i]), _cuAdd(ramo.z[i], aux)); 

					float_type aux2 = atof(line.substr(40, 10).c_str());   
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

void calcYbus2(sistema& sistema, barra& barra, ramo& ramo) {
	for (int i = 1; i <= sistema.nL; i++) {
		complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); 
		aux = _cuMul(_mkComplex(-1., 0.), aux); 
		complex_type aux1 = _cuDiv(aux, ramo.z[IDX1F(i)]);	
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux1; 

		aux = _cuMul(ramo.z[IDX1F(i)], aux);
		aux1 = _cuDiv(_mkComplex(1., 0.), aux); 
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux1;

		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); 
		dAux = dAux * dAux; 

		aux = _mkComplex(dAux, 0.); 
		aux = _cuDiv(aux, ramo.z[IDX1F(i)]); 
		aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); 
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
		sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

		aux = _mkComplex(1 / dAux, 0.); 
		
		aux = _cuDiv(aux, ramo.z[IDX1F(i)]); 
		aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); 
		aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
		sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
	}

	for (int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(0., barra.bsh[IDX1F(i)])); 
	}
}

void calcYbus(sistema& sistema, barra& barra, ramo& ramo) {
		for (int i = 1; i <= sistema.nL; i++) {
			complex_type aux = _cuCon(ramo.tap[IDX1F(i)]); 
			aux = _cuDiv(_mkComplex(-1., 0.), aux); 

			complex_type aux2 = _cuDiv(_mkComplex(-1., 0.), ramo.tap[IDX1F(i)]); 
			aux2 = _cuDiv(aux2, ramo.z[IDX1F(i)]); 

			sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = _cuDiv(aux, ramo.z[IDX1F(i)]); 
			sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux2;

			float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]); 
			dAux = 1 / (dAux * dAux); 

			aux = _mkComplex(dAux, 0.); 
			aux = _cuDiv(aux, ramo.z[IDX1F(i)]); 
			aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); 
			aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
			sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;

			aux = _cuDiv(_mkComplex(1., 0.), ramo.z[IDX1F(i)]);
			aux = _cuAdd(sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)], aux); 
			aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); 
			sistema.Y[IDX2F(ramo.para[IDX1F(i)], ramo.para[IDX1F(i)], sistema.nB)] = aux;
		}

	for (int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])); 
	}

}

void calcYbusSp_eficinte(sistema& sistema, barra& barra, ramo& ramo) {
	std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	for (int i = 1; i <= sistema.nL; i++) {
		std::complex<float_type> aux = std::complex <float_type>(1., 0.) / (std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), -aux / std::conj(std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y))));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), -aux / std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y)));

		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]);
		dAux = 1. / (dAux * dAux);

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux * aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));
	}

	for (int i = 1; i <= sistema.nB; i++) {
		
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])));
	}

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);
	sistema.spY->setFromTriplets(valores.begin(), valores.end());

	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();

	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	for (int i = 0; i < sistema.nnzY; i++) {
		sistema.spYval[i].x = real(sistema.spYvalE[i]);
		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	}

	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] ;
	}

	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] ;
	}

	sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	{
		int acc = 0;
		
		for (int j = 0; j < sistema.spY->outerSize(); j++) {
			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
				sistema.cooRowIndY[acc] = j ;
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

	for (int i = 1; i <= sistema.nL; i++) {
		
		triplCf.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.de[IDX1F(i)]), 1));
		triplCt.emplace_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.para[IDX1F(i)]), 1));

		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		auto daux = _cuAbs(ramo.tap[IDX1F(i)]); 
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

	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();

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
	std::vector<MKL_INT>                  cooCt_rows, cooCf_rows, cooYff_rows, cooYft_rows, cooYtf_rows, cooYtt_rows, cooYsh_rows,
	                                      cooCt_cols, cooCf_cols, cooYff_cols, cooYft_cols, cooYtf_cols, cooYtt_cols, cooYsh_cols;
	std::vector<std::complex<float_type>> cooCt_vals, cooCf_vals, cooYff_vals, cooYft_vals, cooYtf_vals, cooYtt_vals, cooYsh_vals;

	for (int i = 1; i <= sistema.nL; i++) {
		
		cooCf_rows.push_back(IDX1F(i));
		cooCf_cols.push_back(IDX1F(ramo.de[IDX1F(i)]));
		cooCf_vals.push_back(1);

		cooCt_rows.push_back(IDX1F(i));
		cooCt_cols.push_back(IDX1F(ramo.para[IDX1F(i)]));
		cooCt_vals.push_back(1);

		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));

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

	for (int i = 1; i <= sistema.nB; i++) {
		
		cooYsh_rows.push_back(IDX1F(i));
		cooYsh_cols.push_back(IDX1F(i));
		cooYsh_vals.emplace_back(std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]));
	}

	sparse_status_t status;
	sparse_matrix_t cooCt, cooCf, cooYff, cooYft, cooYtf, cooYtt, cooYsh, cooYf, cooYt, aux1, aux2, aux3;
	status = mkl_sparse_z_create_coo(&cooCt, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nB, cooCt_rows.size(), cooCt_rows.data(), cooCt_cols.data(), cooCt_vals.data());
	status = mkl_sparse_z_create_coo(&cooCf, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nB, cooCf_rows.size(), cooCf_rows.data(), cooCf_cols.data(), cooCf_vals.data());
	status = mkl_sparse_z_create_coo(&cooYff, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYff_rows.size(), cooYff_rows.data(), cooYff_cols.data(), cooYff_vals.data());
	status = mkl_sparse_z_create_coo(&cooYft, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYft_rows.size(), cooYft_rows.data(), cooYft_cols.data(), cooYft_vals.data());
	status = mkl_sparse_z_create_coo(&cooYtf, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYtf_rows.size(), cooYtf_rows.data(), cooYtf_cols.data(), cooYtf_vals.data());
	status = mkl_sparse_z_create_coo(&cooYtt, SPARSE_INDEX_BASE_ZERO, sistema.nL, sistema.nL, cooYtt_rows.size(), cooYtt_rows.data(), cooYtt_cols.data(), cooYtt_vals.data());
	status = mkl_sparse_z_create_coo(&cooYsh, SPARSE_INDEX_BASE_ZERO, sistema.nB, sistema.nB, cooYsh_rows.size(), cooYsh_rows.data(), cooYsh_cols.data(), cooYsh_vals.data());

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

	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYff, cooCf, &aux1);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("aux1:\n");
		printMKL_z_csr(&aux1);
	}

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

	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYtf, cooCf, &aux1);

	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, cooYtt, cooCt, &aux2);
	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux1, alpha, aux2, &cooYt);
	mkl_sparse_destroy(aux1);
	mkl_sparse_destroy(aux2);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("cooYt:\n");
		printMKL_z_csr(&cooYt);
	}

	status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, cooCf, cooYf, &aux1);

	status = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, cooCt, cooYt, &aux2);
	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux1, alpha, aux2, &aux3);
	mkl_sparse_destroy(aux1);
	mkl_sparse_destroy(aux2);

	status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, aux3, alpha, cooYsh, &aux1);

	status = mkl_sparse_order(aux1);

	if (global::verbose_mode && !global::laconic_mode)
	{
		printf("spY:\n");
		printMKL_z_csr(&aux1);
	}

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

    sparse_index_base_t indexing;
    MKL_INT rows, cols, *rows_start, *rows_end, *col_indx;
    MKL_Complex16 *values;

    status = mkl_sparse_z_export_csr(aux1, &indexing, &rows, &cols, &rows_start,
        &rows_end, &col_indx, &values);

    if(status != SPARSE_STATUS_SUCCESS) {
		printf("Error!\n");
	}

	auto x = new Eigen::Map< const Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor> >(
		rows, cols, rows_end[rows - 1] - indexing, rows_start, col_indx, values
	);

	sistema.nnzY = rows_end[rows - 1] - indexing;
	sistema.spYvalE = values;

	sistema.spYval = (complex_type*)malloc(sistema.nnzY * sizeof(complex_type));
	for (int i = 0; i < sistema.nnzY; i++) {
		sistema.spYval[i].x = real(sistema.spYvalE[i]);
		sistema.spYval[i].y = imag(sistema.spYvalE[i]);
	}

	sistema.csrColIndY = col_indx; 

	sistema.csrRowPtrY = rows_start;

	mkl_sparse_destroy(aux3);
}