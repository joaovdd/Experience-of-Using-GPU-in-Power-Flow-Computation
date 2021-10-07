#include "ramo.h"
#include "../idx.h"

void InitCsrPhi(sistemaType& sistema, ramoType& ramo) {
	std::vector<Eigen::Triplet<float_type>> valores;
	for (int i = 0; i < sistema.nL; i++) {
		valores.push_back(Eigen::Triplet<float_type>(IDX1F(ramo.de[i]), IDX1F(ramo.para[i]), ramo.phi[i]));
	}

	ramo.eigen_phi.resize(sistema.nL, sistema.nL);

	ramo.eigen_phi.setFromTriplets(valores.begin(), valores.end());

}

void initBranch(sistemaType &sistema, ramoType &ramo){
	ramo.z   = (complex_type *)malloc(sistema.nL * sizeof(complex_type));
	ramo.bsh = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.tap = (complex_type *)malloc(sistema.nL * sizeof(complex_type));
	ramo.phi  = (float_type *)malloc(sistema.nL * sizeof(float_type));

	ramo.de	  = (int *)malloc(sistema.nL * sizeof(float_type));
	ramo.para = (int *)malloc(sistema.nL * sizeof(float_type));

	ramo.Pdp = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Qdp = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Ppd = (float_type *)malloc(sistema.nL * sizeof(float_type));
	ramo.Qpd = (float_type *)malloc(sistema.nL * sizeof(float_type));

	for (int i = 0; i < sistema.nL; i++){
		ramo.z[i] = _mkComplex(0., 0.);
	}
	for (int i = 0; i < sistema.nL; i++){
		ramo.bsh[i] = 0.;
	}
	for (int i = 0; i < sistema.nL; i++){
		ramo.tap[i] = _mkComplex(1., 0.);
	}
	for (int i = 0; i < sistema.nL; i++){
		ramo.de[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++){
		ramo.para[i] = 0;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Pdp[i] = 0.;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Qdp[i] = 0.;
	}
	for (int i = 0; i < sistema.nL; i++) {
		ramo.Ppd[i] = 0.;
	}
	for (int i = 0; i < sistema.nL; i++) {
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