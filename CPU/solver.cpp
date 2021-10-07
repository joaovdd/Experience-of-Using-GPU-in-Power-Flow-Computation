#include "solver.h"

bool solver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	switch (global::metodo)
	{
	case denso:
		dnEigenSolver(sistema, barra, ramo, iterativo);
		break;
	case denso_LAPACKE:
		
		break;
	case esparsoSimples:
		spEigenSolverNaive(sistema, barra, ramo, iterativo);
		break;
	case esparso:
		return spSolveMKL(sistema, barra, ramo, iterativo);
		
		break;
	default:
		printf("ERRO [solver] metodo inválido!\n");
		break;
	}
	return 0; 
}