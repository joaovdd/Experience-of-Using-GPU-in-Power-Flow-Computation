#include "dnCalcYbus.h"

// dp_k: |t|^2 * y_km + jB_sh^km/2
// dp_m:         y_km + jB_sh^km/2
// km  : -t* * y_km
// mk  : -t  * y_km

void calcYbus(sistemaType &sistema, barraType &barra, ramoType &ramo) {
	//inicializar Y - (check 2020)

	//percorre ramos
	for (unsigned int i = 1; i <= sistema.nL; i++) {
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
	for (unsigned int i = 1; i <= sistema.nB; i++) {
		sistema.Y[IDX2F(i, i, sistema.nB)] = _cuAdd(sistema.Y[IDX2F(i, i, sistema.nB)], _mkComplex(barra.gsh[IDX1F(i)], barra.bsh[IDX1F(i)])); // bsh_k + SUM[|t_km|^2*y_km + bsh_km/2]
	}
}