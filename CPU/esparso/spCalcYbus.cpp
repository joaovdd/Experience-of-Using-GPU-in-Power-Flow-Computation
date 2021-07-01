#include "spCalcYbus.h"
#include <iostream>

void calcYbusSp_Matpower(sistemaType& sistema, barraType& barra, ramoType& ramo) {
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplCt, triplCf;
	std::vector<Eigen::Triplet<std::complex<float_type>>> triplYff, triplYft, triplYtf, triplYtt, triplYsh;

	//percorre ramos
	//using namespace std::complex_literals;
	for (int i = 1; i <= sistema.nL; i++) {
		// construção das matrizes de conexao
		triplCf.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.de[IDX1F(i)]), 1));
		triplCt.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(ramo.para[IDX1F(i)]), 1));

		// construcao das matrizes diagonais
		auto aux = _cuAdd(_cuDiv(_mkComplex(1, 0), ramo.z[IDX1F(i)]), _mkComplex(0, ramo.bsh[IDX1F(i)] / 2.));
		auto daux = _cuAbs(ramo.tap[IDX1F(i)]); // cuCmul(_cuCon(ramo.tap[IDX1F(i)]), ramo.tap[IDX1F(i)]);
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

		//aux = _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)]);
		//triplYsh.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i),
		//	std::complex<float_type>(_cuReal(aux), _cuImag(aux))));
	}

	//percorre barras
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

	//auto diff = ans - *(sistema.spY);
	//std::cout << diff << std::endl;

	//std::cout << Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Cf.transpose()) << std::endl;
	//std::cout << Cf << std::endl;
	//std::cout << Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Ct.transpose()) << std::endl;
	//std::cout << Ct << std::endl;


	/*std::cout << aux.outerSize() << "; " << aux.innerSize() << "\n";
	std::cout << aux1.outerSize() << "; " << aux1.innerSize() << "\n";*/
	//std::cout << ans << "\n\n";
	//std::cout << *sistema.spY << "\n\n";
	//std::cout << Cf << "\n\n" << Ct << "\n\n" << Ysh << "\n\n";
	//std::cout << Yff << "\n\n" << Yft << "\n\n" << Ytf << "\n\n" << Ytt;
	// std::cout << Yf << "\n\n"; //<< Yf << "\n\n" << aux << "\n\n" << aux1;
	// std::cout << ans << "\n"; //<< aux2 << "\n" << Ysh;

	//auto aux2 = aux + aux1 + Ysh;

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);

	//*(sistema.spY) = Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Cf.transpose()) * Yf + Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(Ct.transpose()) * Yt + Ysh;
	*(sistema.spY) = ans;

	// Cópia da outra função

	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();
	//sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
	//sistema.csrColIndY = sistema.spY->outerIndexPtr();

	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] /*+ 1*/;
		// std::cout << sistema.csrColIndY[j] << ' ';
	}
	// std::cout << std::endl;

	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] /*+ 1*/;
		// std::cout << sistema.csrRowPtrY[j] << ' ';
	}

	sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	{
		int acc = 0;
		// csr to coo
		for (int j = 0; j < sistema.spY->outerSize(); j++) {
			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
				sistema.cooRowIndY[acc] = j /*+ 1*/;
				// std::cout << j /*+ 1*/ << ' ';
			}
		}
		//std::cout << std::endl;
	}
}

//void calcYbusSp(sistemaType &sistema, barraType &barra, ramoType &ramo) {
//	//inicializar Y
//
//	//using namespace std::complex_literals;
//
//	Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);
//
//	//percorre ramos
//	for (int i = 1; i <= sistema.nL; i++) {
//		
//		std::complex<float_type> aux = std::complex<float_type>(-1., 0.) / (std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y)) * std::conj(std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y));
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
//		dAux = 1/(dAux * dAux); // 1/|t_ft|^2
//
//		dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = dAux/std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y) + dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.);
//
//		// aux = _mkComplex(dAux, 0.); // 1/|t_ft|^2
//		// aux = _cuDiv(aux, ramo.z[IDX1F(i)]); // y_ft/|t_ft|^2
//		// aux = _cuAdd(sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)], aux); // Y_ff + |t_ft|^2*y_ft
//		// aux = _cuAdd(_mkComplex(0., ramo.bsh[IDX1F(i)] / 2.), aux); // Y_ff + |t_ft|^2*y_ft + j*bsh_ft/2
//		// // sistema.Y[IDX2F(ramo.de[IDX1F(i)], ramo.de[IDX1F(i)], sistema.nB)] = aux;
//		// dnY(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)])) = aux.x + aux.y*1i;
//
//		dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) = 1./(ramo.z[IDX1F(i)].x + ramo.z[IDX1F(i)].y*1i) + dnY(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)])) + (ramo.bsh[IDX1F(i)] * 1i) / 2.;
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
//		dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)]*1i;
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
//	sistema.csrColIndY= (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
//	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
//		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] /*+ 1*/;
//		// std::cout << sistema.csrColIndY[j] << ' ';
//	}
//	// std::cout << std::endl;
//
//	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
//	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
//		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] /*+ 1*/;
//		// std::cout << sistema.csrRowPtrY[j] << ' ';
//	}
//	// std::cout << std::endl;
//
//
//
//	//if(verbose_mode){
//	//	printf("\n\nType punning: std::complex<float_type> -> complex_type?\n\n");
//	//}
//
//	sistema.cooRowIndY = (int*)malloc(dnY.nonZeros() * sizeof(int));
//	{
//		int acc = 0;
//		// csr to coo
//		for (int j = 0; j < sistema.spY->outerSize(); j++){
//			for (; acc < sistema.spY->outerIndexPtr()[j+1]; acc++){
//				sistema.cooRowIndY[acc] = j /*+ 1*/;
//				// std::cout << j /*+ 1*/ << ' ';
//			}
//		}
//		//std::cout << std::endl;
//
//	}
//}

// dp_k: |t|^2 * y_km + jB_sh^km/2
// dp_m:         y_km + jB_sh^km/2
// km  : -t* * y_km
// mk  : -t  * y_km

void calcYbusSp_eficinte(sistemaType& sistema, barraType& barra, ramoType& ramo) {
	//inicializar Y - (check 2020)

	//using namespace std::complex_literals;

	std::vector<Eigen::Triplet<std::complex<float_type>>> valores;

	//Eigen::MatrixXcd dnY = Eigen::MatrixXcd::Zero(sistema.nB, sistema.nB);

	//percorre ramos
	for (int i = 1; i <= sistema.nL; i++) {

		std::complex<float_type> aux = std::complex < float_type>(1.,0) / (std::complex<float_type>(ramo.z[IDX1F(i)].x, ramo.z[IDX1F(i)].y));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), - aux / std::conj(std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y))));
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), - aux / std::complex<float_type>(ramo.tap[IDX1F(i)].x, ramo.tap[IDX1F(i)].y)));


		float_type dAux = _cuAbs(ramo.tap[IDX1F(i)]);
		dAux = 1. / (dAux * dAux);

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.de[IDX1F(i)]), IDX1F(ramo.de[IDX1F(i)]), dAux * aux + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));

		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(ramo.para[IDX1F(i)]), IDX1F(ramo.para[IDX1F(i)]), aux  + std::complex<float_type>(0., ramo.bsh[IDX1F(i)] / 2.)));
	}

	//percorre barras
	for (int i = 1; i <= sistema.nB; i++) {
		//dnY(IDX1F(i), IDX1F(i)) = dnY(IDX1F(i), IDX1F(i)) + barra.bsh[IDX1F(i)] * 1i;
		valores.push_back(Eigen::Triplet<std::complex<float_type>>(IDX1F(i), IDX1F(i), std::complex<float_type>(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])));
	}

	sistema.spY = new Eigen::SparseMatrix<std::complex<float_type>, Eigen::StorageOptions::RowMajor>(sistema.nB, sistema.nB);
	//*sistema.spY = dnY.sparseView(); // converte para esparsa;
	sistema.spY->setFromTriplets(valores.begin(), valores.end());
	
	// https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html#ac2684952b14b5c9b0f68ae3bb8c517a6
	sistema.nnzY = sistema.spY->nonZeros();
	sistema.spYvalE = sistema.spY->valuePtr();
	//sistema.csrRowPtrY = sistema.spY->innerIndexPtr();
	//sistema.csrColIndY = sistema.spY->outerIndexPtr();

	sistema.csrColIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	for (int j = 0; j < sistema.spY->nonZeros(); j++) {
		sistema.csrColIndY[j] = sistema.spY->innerIndexPtr()[j] /*+ 1*/;
		// std::cout << sistema.csrColIndY[j] << ' ';
	}
	// std::cout << std::endl;

	sistema.csrRowPtrY = (int*)malloc((sistema.spY->outerSize() + 1) * sizeof(int));
	for (int j = 0; j <= sistema.spY->outerSize(); j++) {
		sistema.csrRowPtrY[j] = sistema.spY->outerIndexPtr()[j] /*+ 1*/;
		// std::cout << sistema.csrRowPtrY[j] << ' ';
	}

	sistema.cooRowIndY = (int*)malloc(sistema.spY->nonZeros() * sizeof(int));
	{
		int acc = 0;
		// csr to coo
		for (int j = 0; j < sistema.spY->outerSize(); j++) {
			for (; acc < sistema.spY->outerIndexPtr()[j + 1]; acc++) {
				sistema.cooRowIndY[acc] = j /*+ 1*/;
				// std::cout << j /*+ 1*/ << ' ';
			}
		}
		//std::cout << std::endl;
	}
}

