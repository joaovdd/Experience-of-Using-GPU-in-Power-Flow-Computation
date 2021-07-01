#include "ramo.h"
#include "../idx.h"

// popula estrutura csr com os dados de defasagem angular dos ramos do sistema e inicializa vetores
void InitCsrPhi(sistemaType& sistema, ramoType& ramo) {

	std::vector<Eigen::Triplet<float_type>> valores;
	for (unsigned int i = 0; i < sistema.nL; i++) {
		valores.push_back(Eigen::Triplet<float_type>(IDX1F(ramo.de[i]), IDX1F(ramo.para[i]), ramo.phi[i]));
	}

	// ramo.eigen_phi = new Eigen::SparseMatrix<float_type, Eigen::StorageOptions::RowMajor>(sistema.nL, sistema.nL);

	ramo.eigen_phi.resize(sistema.nL, sistema.nL);

	ramo.eigen_phi.setFromTriplets(valores.begin(), valores.end());

	//ramo.csrColIndPhi = (int*)malloc(ramo.phi->nonZeros() * sizeof(int));
	//for (int j = 0; j < ramo.eigen_phi->nonZeros(); j++) {
	//	ramo.csrColIndPhi[j] = ramo.eigen_phi->innerIndexPtr()[j] /*+ 1*/;
	//}

	//ramo.csrRowPtrPhi = (int*)malloc((ramo.eigen_phi->outerSize() + 1) * sizeof(int));
	//for (int j = 0; j <= ramo.eigen_phi->outerSize(); j++) {
	//	ramo.csrRowPtrPhi[j] = ramo.eigen_phi->outerIndexPtr()[j] /*+ 1*/;
	//}
}

void initBranch(sistemaType &sistema, ramoType &ramo){
	ramo.z   = (complex_type *)malloc(sistema.nL * sizeof(complex_type));
	ramo.bsh = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.tap = (complex_type *)malloc(sistema.nL * sizeof(complex_type));
	ramo.phi  = (float_type *)malloc(sistema.nL * sizeof(float_type));
	
	ramo.de	  = (unsigned int *)malloc(sistema.nL * sizeof(float_type));
	ramo.para = (unsigned int *)malloc(sistema.nL * sizeof(float_type));

	ramo.Pdp = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Qdp = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Ppd = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Qpd = (float_type *)malloc(sistema.nL * sizeof(float_type));


	for (unsigned int i = 0; i < sistema.nL; i++){
		ramo.z[i] = _mkComplex(0., 0.);
	}
	for (unsigned int i = 0; i < sistema.nL; i++){
		ramo.bsh[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nL; i++){
		ramo.tap[i] = _mkComplex(1., 0.);
	}
	for (unsigned int i = 0; i < sistema.nL; i++){
		ramo.de[i] = 0;
	}
	for (unsigned int i = 0; i < sistema.nL; i++){
		ramo.para[i] = 0;
	}
	for (unsigned int i = 0; i < sistema.nL; i++) {
		ramo.Pdp[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nL; i++) {
		ramo.Qdp[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nL; i++) {
		ramo.Ppd[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nL; i++) {
		ramo.Qpd[i] = 0.;
	}
}

void finBranch(ramoType &ramo) {
	if (ramo.z == nullptr) { free(ramo.z); }
	if (ramo.bsh == nullptr) { free(ramo.bsh); }
	if (ramo.tap == nullptr) { free(ramo.tap); }
	if (ramo.phi == nullptr) { free(ramo.phi); }
	if (ramo.Pdp == nullptr) { free(ramo.Pdp); }
	if (ramo.Qdp == nullptr) { free(ramo.Qdp); }
	if (ramo.Ppd == nullptr) { free(ramo.Ppd); }
	if (ramo.Qpd == nullptr) { free(ramo.Qpd); }

	if (ramo.de == nullptr) { free(ramo.de); }
	if (ramo.para == nullptr) { free(ramo.para); }
}