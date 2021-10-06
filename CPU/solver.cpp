#include "solver.h"

// soluciona o sistema linear. Se bem sucedida retorna 0 (apenas esparsa).
bool solver(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	switch (global::metodo)
	{
	case denso:
		dnEigenSolver(sistema, barra, ramo, iterativo);
		break;
	case denso_LAPACKE:
		//LAPACKE_solver(sistema, barra, ramo, iterativo);
		break;
	case esparsoSimples:
		spEigenSolverNaive(sistema, barra, ramo, iterativo);
		break;
	case esparso:
		return spSolveMKL(sistema, barra, ramo, iterativo);
		// return spSolve(sistema, barra, ramo, iterativo);
		break;
	default:
		printf("ERRO [solver] metodo inválido!\n");
		break;
	}
	return 0; // implementar, ponto boleano para a densa também caso necessário
}