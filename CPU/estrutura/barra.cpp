#include "barra.h"

void initBus(sistemaType &sistema, barraType &barra){
	barra.id = (unsigned int *)malloc(sistema.nB * sizeof(unsigned int));
	barra.V		= (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.theta = (float_type *)malloc(sistema.nB * sizeof(float_type));

	barra.Pliq	= (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.Qliq	= (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.Pload = (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.Qload = (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.Pg	= (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.Qg	= (float_type *)malloc(sistema.nB * sizeof(float_type));

	barra.Vbase = (float_type *)malloc(sistema.nB * sizeof(float_type));

	barra.gsh	= (float_type *)malloc(sistema.nB * sizeof(float_type));
	barra.bsh	= (float_type *)malloc(sistema.nB * sizeof(float_type));

	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.V[i] = global::v_inicial;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.theta[i] = global::theta_inicial;
	}

	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Pliq[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Qliq[i] = 0.;
	}

	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Pload[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Qload[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Pg[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Qg[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.Vbase[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.gsh[i] = 0.;
	}
	for (unsigned int i = 0; i < sistema.nB; i++){
		barra.bsh[i] = 0.;
	}
}

void finBus(barraType &barra) {
	if (barra.id == nullptr) { free(barra.id); }
	if (barra.V == nullptr) { free(barra.V); }
	if (barra.theta == nullptr) { free(barra.theta); }
												 
	if (barra.Pliq == nullptr) { free(barra.Pliq); }
	if (barra.Qliq == nullptr) { free(barra.Qliq); }
	if (barra.Pload == nullptr) { free(barra.Pload); }
	if (barra.Qload == nullptr) { free(barra.Qload); }
	if (barra.Pg == nullptr) { free(barra.Pg); }
	if (barra.Qg == nullptr) { free(barra.Qg); }
												 
	if (barra.Vbase == nullptr) { free(barra.Vbase); }
												 
	if (barra.gsh == nullptr) { free(barra.gsh);	}
	if (barra.bsh == nullptr) { free(barra.bsh);	}
}