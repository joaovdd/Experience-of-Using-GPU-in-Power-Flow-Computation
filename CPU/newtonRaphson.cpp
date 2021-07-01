#include "newtonRaphson.h"

using namespace std;
using namespace Eigen;


// /* Auxiliary routine: printing a matrix */
// void print_matrix( char* desc, int m, int n, float_type* a, int lda ) {
//         int i, j;
//         printf( "\n %s\n", desc );
//         for( i = 0; i < m; i++ ) {
//                 for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
//                 printf( "\n" );
//         }
// }

// /* Auxiliary routine: printing a vector of integers */
// void print_int_vector( char* desc, int n, int* a ) {
//         int j;
//         printf( "\n %s\n", desc );
//         for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
//         printf( "\n" );
// }

float_type maxi(float_type* vet, unsigned int dim){
	float_type aux = 0;
	for (unsigned int i = 0; i < dim; i++){
		if (abs(vet[i]) > aux) {
			aux = abs(vet[i]);
			//printf("oi %d\n", aux);
		}
	}
	return aux;
}

// Calcula o resíduo
void calcRes(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo){
	// deltaPs
	for (unsigned int i = 1; i < sistema->barraVO; i++){ //antes da swing; i percorre barras
		// deltaP = |Pcalc - Pesp|
		iterativo->g[IDX1F(i)] = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
		//printf("g[%d] = ", i-1);
		//printf("%e = |%e - %e| ", iterativo->g[IDX1F(i)], iterativo->Pcalc[IDX1F(i)], barra->Pliq[IDX1F(i)]);
		//printf("g[%d] = %e\n", i, iterativo->g[IDX1F(i)]);
	}

	for (unsigned int i = sistema->barraVO+1; i <= sistema->nB; i++){ // depois da swing; i percorre barras
		// deltaP = |Pcalc - Pesp|
		float_type aux = barra->Pliq[IDX1F(i)] - iterativo->Pcalc[IDX1F(i)];
		iterativo->g[IDX1F(i) - 1] = aux;
		//printf("g[%d] = ", i-1-1);
		//printf("%e = |%14.6e - %14.6e| ", iterativo->g[IDX1F(i)-1], iterativo->Pcalc[IDX1F(i)], barra->Pliq[IDX1F(i)]);
		//printf("g[%d] = %e\n", i-1, iterativo->g[IDX1F(i) - 1]);
	}

	// deltaQs
	const unsigned int offset = sistema->nB-1;
	for (unsigned int i = 0; i < sistema->nPQ; i++){ // i � iterador das barras PQ
		// deltaQ = |Qcalc - Qesp|
		iterativo->g[i+offset] = barra->Qliq[IDX1F(sistema->barrasPQ[i])] - iterativo->Qcalc[IDX1F(sistema->barrasPQ[i])];
		//printf("g[%d] = ", i+sistema->nPQ+sistema->nPV);
		//printf("%e = |%14.6e - %14.6e| ", iterativo->g[i+sistema->nPQ+sistema->nPV], iterativo->Qcalc[IDX1F(sistema->barrasPQ[i])], barra->Qliq[IDX1F(sistema->barrasPQ[i])]);
		//printf("g[%d] = %14.6e\n", i+sistema->nPQ+sistema->nPV, iterativo->g[i+sistema->nPQ+sistema->nPV]);
	}
}

void attVO(sistemaType* sistema, barraType* barra, iterativoType* iterativo){
	for (unsigned int i = 1; i < sistema->barraVO; i++){ // i percorre as barras diretamente e pula a barra swing
		barra->theta[IDX1F(i)] += iterativo->g[IDX1F(i)];  // atualiza valor de theta; g � igual a -deltaX!
		//printf("attVO: theta[%d] = %f\n", i, barra->theta[IDX1F(i)]);
	}
	for (unsigned int i = sistema->barraVO+1; i <= sistema->nB; i++) { // i percorre as barras diretamente e pula a barra swing
		barra->theta[IDX1F(i)] += iterativo->g[IDX1F(i)-1];              // atualiza valor de theta; g � igual a -deltaX!
		//printf("attVO: theta[%d] = %f\n", i, barra->theta[IDX1F(i)]);
	}

	const unsigned int offset = sistema->nB - 1;
	for (unsigned int i = 0; i < sistema->nPQ; i++) {                    // i � iterador das barras PQ
		barra->V[IDX1F(sistema->barrasPQ[i])] += iterativo->g[i + offset]; // atualiza valor de V; g � igual a -deltaX!
		//printf("attVO: V[%d] = %f\n", sistema->barrasPQ[i], barra->V[IDX1F(sistema->barrasPQ[i])]);
	}

	if (global::verbose_mode) {
		printf("\nV^(%d):\n", iterativo->iteracao);
		showVec(barra->V, sistema->nB, 6);
		printf("\ntheta^(%d):\n", iterativo->iteracao);
		showVec(barra->theta, sistema->nB, 6);
	}

}

void evalLimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	chkLimQ(sistema, barra, ramo, iterativo);   // atualiza sistema->limQ
	geraVetPV(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPVlim e iterativo->nPQlim
	geraVetPQ(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPQlim e iterativo->nPQlim  [[ver openmp while]]
	iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; // atualiza o tamanho de g com os limites
	//calculePQ(sistema, barra, ramo, iterativo); // aualiza potências

#if JACOBIANO_ESPARSO_STENCIL
	if (global::metodo == esparso) {
		// verificar se ocorreu mudanca no chkLimQ

		{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
			//zera Jac
			sistema->spJval.clear();
			sistema->spJsubmatType.clear();
			sistema->cooColIndSubMatJ.clear();
			sistema->cooRowIndSubMatJ.clear();
			sistema->cooColIndJ.clear();
			sistema->cooRowIndJ.clear();
			//sistema->csrRowPtrJ.clear();

			Jstencil0based(sistema, barra, ramo, iterativo);
		}
	}
#endif
}

void do_LimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	do_chkLimQ(sistema, barra, ramo, iterativo);   // atualiza sistema->limQ
	geraVetPV(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPVlim e iterativo->nPQlim
	geraVetPQ(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPQlim e iterativo->nPQlim  [[ver openmp while]]
	iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; // atualiza o tamanho de g com os limites
	//calculePQ(sistema, barra, ramo, iterativo); // aualiza potências

#if JACOBIANO_ESPARSO_STENCIL
	//#ifdef ESPARSO
	if (global::metodo == esparso) {
		// verificar se ocorreu mudanca no chkLimQ

		{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
			//zera Jac
			sistema->spJval.clear();
			sistema->spJsubmatType.clear();
			sistema->cooColIndSubMatJ.clear();
			sistema->cooRowIndSubMatJ.clear();
			sistema->cooColIndJ.clear();
			sistema->cooRowIndJ.clear();
			//sistema->csrRowPtrJ.clear();

			Jstencil0based(sistema, barra, ramo, iterativo);
		}

	}
#endif
}

void undo_LimInjReat(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	if (undo_chkLimQ(sistema, barra, ramo, iterativo)) {
		geraVetPV(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPVlim e iterativo->nPQlim
		geraVetPQ(sistema, barra, ramo, iterativo); // atualiza sistema->barrasPQlim e iterativo->nPQlim  [[ver openmp while]]
		iterativo->ngLim = iterativo->nPQlim + iterativo->nPQlim + iterativo->nPVlim; // atualiza o tamanho de g com os limites
		//calculePQ(sistema, barra, ramo, iterativo); // aualiza potências

#if JACOBIANO_ESPARSO_STENCIL
		if (global::metodo == esparso) {
			// verificar se ocorreu mudanca no chkLimQ

			{ BENCHMARK_JACOBIANOSTENCIL_REBUILD
				//zera Jac
				sistema->spJval.clear();
				sistema->spJsubmatType.clear();
				sistema->cooColIndSubMatJ.clear();
				sistema->cooRowIndSubMatJ.clear();
				sistema->cooColIndJ.clear();
				sistema->cooRowIndJ.clear();
				//sistema->csrRowPtrJ.clear();

				Jstencil0based(sistema, barra, ramo, iterativo);
			}

		}
#endif
	}
}

// retorna falso caso haja violação de dimite de injeção de reativo, verdadeiro, caso contrário
bool simpleChkLimQ(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo) {
	//printf("\n");
	auto folga = global::tol_limInjReat; //0.01; //global::tol
//#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 0; i < sistema->nPV; i++) { // percorre barras PV
		// Qg = Qliq - Qload
		float_type aux = (iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])]);
		if (aux < (sistema->limQinf[i] - folga)) {
			return false;
		}
		else if (aux > (sistema->limQsup[i] + folga)) {
			return false;
		}
	}
	return true;
}

bool hasTrue(bool* arr, unsigned int size) {
//#pragma omp parallel for if (global::openmp)
	for (unsigned int i = 0; i < size; i++) { // percorre barras PV
		if (arr[i]) {
			return true;
		}
	}
	return false;
}

float_type accumulate(float_type* arr, unsigned int size) {
	float_type ans = 0;
	for (unsigned int i = 0; i < size; i++) { // percorre barras PV
		ans += abs(arr[i]);
	}
	return ans;
}

void nR(sistemaType* sistema, barraType* barra, ramoType* ramo, iterativoType* iterativo){
	//bool flgLimInjReat = false;
	
	int n = sistema->nPQ + sistema->nPV + sistema->nPQ;

	initLim(sistema, barra, ramo, iterativo);
	if (global::verbose_mode && global::lim_inj_reat) {
		printAlLim(*sistema, *barra, *ramo, *iterativo);
	}

	{ BENCHMARK_PROCESSOITERATIVO
		for (iterativo->iteracao = 0; iterativo->iteracao < global::no_max_iter; iterativo->iteracao++) {
			if (global::verbose_mode) {
				printf("\n****************************************************************** ITERACAO %d ******************************************************************\n", iterativo->iteracao);
			}

			{ BENCHMARK_CALCPQ
				switch (global::metodo)
				{
				case denso:
				case esparsoSimples:
					dnCalculePQ(sistema, barra, ramo, iterativo);
					break;
				case esparso:
					spCalcP(*sistema, *barra, *ramo, *iterativo);
					spCalcQ(*sistema, *barra, *ramo, *iterativo);
					break;
				default:
					printf("\n\n[ERRO] calcPQ: metodo inválido!\n\n");
					break;
				}
			}

			//calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g

			//// controle (lim inj reat)
			if ((iterativo->iteracao >= 1) && global::lim_inj_reat) {
				if (maxi(&(iterativo->gLim[0]), iterativo->ngLim) < /*global::tol * 100.*/0.1) { // pula a primeira e segunda iterações e deltaX deve ser inferior a um centésimo da tolerancia
					evalLimInjReat(sistema, barra, ramo, iterativo);
					//calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g
				}
				else if (hasTrue(iterativo->limQ, sistema->nPV) /*&& maxi(&(iterativo->gLim[0]), iterativo->ngLim) < 1 /*global::tol  /*&& 0*/) {
					// se o limInjReat de alguma barra foi violado
					// checar se ainda está violado
					undo_LimInjReat(sistema, barra, ramo, iterativo);
					//calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g
				}
				/*
				if (0) {
					do_LimInjReat(sistema, barra, ramo, iterativo);
				}*/
			}

			calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g

			//if (maxi(&(iterativo->gLim[0]), iterativo->ngLim) > 1000) {
			//	for (int i = 0; i <= sistema->nPQ; i++) {
			//		barra->V[IDX1F(iterativo->barrasPQlim[i])] = 1.;
			//		barra->theta[IDX1F(iterativo->barrasPQlim[i])] = 0.;
			//	}
			//}

			//// controle
			//if ((iterativo->iteracao >= 1) && global::lim_inj_reat) {
			//	if (accumulate(iterativo->gLim, iterativo->ngLim) < .1) { // pula a primeira e segunda iterações e deltaX deve ser inferior a um centésimo da tolerancia
			//		evalLimInjReat(sistema, barra, ramo, iterativo);
			//		calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g
			//	}
			//	else if (hasTrue(iterativo->limQ, sistema->nPV) && accumulate(iterativo->gLim, iterativo->ngLim) < .1) {
			//		// se o limInjReat de alguma barra foi violado
			//		// checar se ainda está violado
			//		undo_LimInjReat(sistema, barra, ramo, iterativo);
			//		calcResLim(sistema, barra, ramo, iterativo); // resíduo = iterativo->g
			//	}
			//}

			if (global::verbose_mode) {
				printf("Pcalc^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->Pcalc, sistema->nB, 4);
				printf("\nQcalc^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->Qcalc, sistema->nB, 4);
				printf("\ndeltaP^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->gLim, sistema->nB - 1, 4);
				printf("\ndeltaQ^(%d) =\n", iterativo->iteracao);
				showVec(iterativo->gLim + sistema->nB - 1, iterativo->nPQlim, 4);

				printf("\nQg_barrasPV^(%d) =\n", iterativo->iteracao);
				cout << '[' << iterativo->Qcalc[IDX1F(sistema->barrasPV[0])] + barra->Qload[IDX1F(sistema->barrasPV[0])] << endl;
				for (unsigned int i = 1; i < (sistema->nPV - 1); i++) {
					cout << ' ' << iterativo->Qcalc[IDX1F(sistema->barrasPV[i])] + barra->Qload[IDX1F(sistema->barrasPV[i])] << endl;
				}
				cout << ' ' << iterativo->Qcalc[IDX1F(sistema->barrasPV[sistema->nPV - 1])] + barra->Qload[IDX1F(sistema->barrasPV[sistema->nPV - 1])] << ']' << endl;

				cout << "\nlimQ:" << endl;
				showVec(iterativo->limQ, sistema->nPV, 4);
				cout << "\nbarrasPVlim:" << endl;
				showVec(iterativo->barrasPVlim, iterativo->nPVlim, 4);
				cout << "\nbarrasPQlim:" << endl;
				showVec(iterativo->barrasPQlim, iterativo->nPQlim, 4);
				cout << "\nlimQinf:" << endl;
				showVec(sistema->limQinf, sistema->nPV, 4);
				cout << "\nlimQsup:" << endl;
				showVec(sistema->limQsup, sistema->nPV, 4);
				cout << "\nQliqLim:" << endl;
				showVec(iterativo->QliqLim, sistema->nB, 4);
			}

			float_type maxG = maxi(iterativo->gLim, iterativo->ngLim); // retorna valor absoluto do elto de maior módulo em g 
			if (global::output_processo_iterativo && !global::laconic_mode) {
				printf("max(gLim) = %e; sum(gLim) = %e\n", maxG, accumulate(iterativo->gLim, iterativo->ngLim));
			}
			if (maxG <= global::tol) {

				//checa Lim inj reat ou, caso esteja deslicado pula
				if (global::lim_inj_reat && !simpleChkLimQ(sistema, barra, ramo, iterativo)) {
					// continue;
					evalLimInjReat(sistema, barra, ramo, iterativo);
					//undo_LimInjReat(sistema, barra, ramo, iterativo);
				}
				else {
					if (global::verbose_mode) {
						printf("\n\nbreak!\nmax(g) = %e\n\n", maxi(iterativo->gLim, iterativo->ngLim));
					}
					return; //break; // terminou!!!
				}

			}
			if (global::verbose_mode) {
				printf("     ...     continua:\n\n");
			}

			switch (global::metodo) {
				case esparso:
	#if JACOBIANO_ESPARSO_STENCIL
				{ BENCHMARK_JACOBIANOSTENCIL_FILL
					// fillJstencil(sistema, barra, ramo, iterativo);
					fillJstencil0based(sistema, barra, ramo, iterativo);
				}
	#else
				{ BENCHMARK_JACOBIANO
					//spCalcJ(sistema, barra, ramo, iterativo);
					spCalcJ_eficiente(sistema, barra, ramo, iterativo);
					//calcJacLim(sistema, barra, ramo, iterativo);
					
					// imprime jacobino COO
					//for (int i = 0; i < sistema->nnzJ; i++) {
					//	std::cout << sistema->cooRowIndJ[i] + 1 << " " << sistema->cooColIndJ[i] + 1 << " " << sistema->spJval[i] << "\n";
					//}


				}
	#endif
					break;
				default:
				{ BENCHMARK_JACOBIANO
					calcJacLim(sistema, barra, ramo, iterativo);
				}
					break;
			}

			{ BENCHMARK_SISTEMALINEAR
				if (solver(sistema, barra, ramo, iterativo)) {
					iterativo->iteracao = global::no_max_iter; // apenas para AIRE
					break; // sistema divergiu
				}
			}

			if(global::verbose_mode){
				printf("\ndeltaX =\n");
				showVec(iterativo->gLim, sistema->nB - 1 + iterativo->nPQlim, 5);
			}
		
			// atualiza theta e V
			attVOlim(sistema, barra, iterativo);

		}
	}

	if (maxi(&(iterativo->gLim[0]), iterativo->nPQlim + iterativo->nPVlim + iterativo->nPQlim) > global::tol) {
		//if (verbose_mode) {
			printf("\n\nO processo iterativo nao convergiu!\nmax(g) = %e\n\n", maxi(iterativo->gLim, iterativo->nPQlim + iterativo->nPVlim + iterativo->nPQlim));
		//}
		//return; //break; // terminou!!!
	}

	// imprimir solucao
}