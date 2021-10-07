#include "spCalcYbus.h"
#include <iostream>

void calcYbusSp_Matpower(sistemaType& sistema, barraType& barra, ramoType& ramo) {
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplCt, triplCf;
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplYff, triplYft, triplYtf, triplYtt, triplYsh;

	for (int i = 1; i <= sistema.nL; i++) {
		
		triplCf.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.de[IDX1F(i)]), 1));
		triplCt.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.para[IDX1F(i)]), 1));

		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		auto daux = _cuAbs(ramo.tap[IDX1F(i)]); 
		aux = _cuMul(aux, _mkComplex(1/(daux*daux), 0));
		triplYff.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), _cuCon(ramo.tap[IDX1F(i)])), ramo.z[IDX1F(i)]);
		triplYft.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuDiv(_cuDiv(_mkComplex(-1, 0), ramo.tap[IDX1F(i)]), ramo.z[IDX1F(i)]);
		triplYtf.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

		aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		triplYtt.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
			std::complex<float_type>(_cuReal(aux), _cuImag(aux))));

	}

	for (int i = 1; i <= sistema.nB; i++) {
		auto aux = _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]);
		triplYsh.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
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

	auto ans = aux2 + Ysh;

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

	*(sistema.spY) = ans;

	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();

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

void calcYbusSp_eficinte(sistemaType& sistema, barraType& barra, ramoType& ramo) {
	std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	for (int i = 1; i <= sistema.nL; i++) {
		std::complex<float_type> aux = std::complex < float_type>(1.,0) / (std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), - aux / std::conj(std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y))));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), - aux / std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y)));

		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]);
		dAux = 1. / (dAux * dAux);

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux * aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux  + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));
	}

	for (int i = 1; i <= sistema.nB; i++) {
		
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])));
	}

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

	sistema.spY->setFromTriplets(valores.begin(), valores.end());

	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();

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

