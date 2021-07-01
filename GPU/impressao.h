#pragma once

#include <stdio.h>
#include "sistema.h"
#include "benchmarks.h"
#include <string>
#include <fstream>
#include <iterator>



// file output
//#include <fstream>
#include <stdio.h>
#include <math.h> 

void impressao(sistema &sistema, barra &barra, ramo &ramo, iterativo &iterativo) {
	printf("\n\n                                 FLUXOS DE POTENCIA\n");
            
	printf("        BUS      V         ANG         P          Q        TO       PIJ          QIJ\n");
	for (unsigned short i = 0; i < sistema.nB; i++)
	{
		printf("      %4d%11.5f%11.5f%11.5f%11.5f  ", i + 1, barra.V[i], barra.theta[i], iterativo.Pcalc[i], iterativo.Qcalc[i]); // barra.Pliq[i], barra.Qliq[i]);
		//Todo: fazer um return com erro antes do while extrapolar os limites dos vetores
		//'**ERROR** UNABLE TO CALCULATE POWER FLOWS'
		unsigned short j;
		for (j = 0; j < sistema.nL; j++)
		{
			if (ramo.de[j] == (i + 1)) {
				printf("%4d  %11.5f  %11.5f\n", ramo.para[j], ramo.Pdp[j], ramo.Qdp[j]);
				break;
			}
			else if (ramo.para[j] == (i + 1)) {
				printf("%4d  %11.5f  %11.5f\n", ramo.de[j], ramo.Ppd[j], ramo.Qpd[j]);
				break;
			}
		}
		j++;
		for (; j < sistema.nL; j++)
		{
			if (ramo.de[j] == (i + 1)) {
				printf("                                                        %4d  %11.5f  %11.5f\n", ramo.para[j], ramo.Pdp[j], ramo.Qdp[j]);
			}
			else if (ramo.para[j] == (i + 1)) {
				printf("                                                        %4d  %11.5f  %11.5f\n", ramo.de[j], ramo.Ppd[j], ramo.Qpd[j]);
			}
		}

	}
}

// imprime Pg e Pl
void impressao2(sistema& sistema, barra& barra, ramo& ramo, iterativo& iterativo) {
	printf("\n\n                                 FLUXOS DE POTENCIA\n");

	printf("        BUS      V         ANG        Pli        Qli        Pge        Qge      TO        PIJ          QIJ\n");
	for (unsigned int i = 0; i < sistema.nB; i++)
	{
		//printf("      %4d%11.5f%11.5f%11.5f%11.5f  ", i + 1, barra.V[i], barra.theta[i], barra.Pliq[i], barra.Qliq[i]);
		printf("      %4d%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f  ", barra.id[i], barra.V[i], barra.theta[i], iterativo.Pcalc[i], iterativo.Qcalc[i], iterativo.Pcalc[i] + barra.Pload[i], iterativo.Qcalc[i] + barra.Qload[i]);
		//Todo: fazer um return com erro antes do while extrapolar os limites dos vetores
		//'**ERROR** UNABLE TO CALCULATE POWER FLOWS'
		unsigned int j;
		for (j = 0; j < sistema.nL; j++)
		{
			if (ramo.de[j] == (i + 1)) {
				printf("%4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.para[j])], ramo.Pdp[j], ramo.Qdp[j]);
				break;
			}
			else if (ramo.para[j] == (i + 1)) {
				printf("%4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.de[j])], ramo.Ppd[j], ramo.Qpd[j]);
				break;
			}
		}
		j++;
		for (; j < sistema.nL; j++)
		{
			if (ramo.de[j] == (i + 1)) {
				printf("                                                                              %4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.para[j])], ramo.Pdp[j], ramo.Qdp[j]);
			}
			else if (ramo.para[j] == (i + 1)) {
				printf("                                                                              %4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.de[j])], ramo.Ppd[j], ramo.Qpd[j]);
			}
		}

	}
}

std::string met2str(metodo met) {
	switch (met)
	{
	case metodo::denso:
		return "denso";
		break;
	case metodo::esparso:
		return "esparso";
		break;
	case metodo::hibridoA:
		return "hibridoA";
		break;
	case metodo::hibridoB:
		return "hibridoB";
		break;
	//case esparsoSimples:
	//	return "esparsosimples";
	//	break;
	default:
		return "nda";
		break;
	}
}

void benchmarkModePrint(iterativo& iterPon, std::chrono::duration<double, std::milli>& duracao) {
	std::string aux("out_"+ met2str(global::metodo) + "-" + global::arq_entrada + ".csv");
	
	FILE* pFile;

	pFile = fopen(aux.c_str(), "a");

	// fprintf(pFile, "%d; iteracoes.; %f; ms.;\n", iterPon.iteracao, duracao.count());
	double auxFrac, auxInt;
	auxFrac = modf(duracao.count(), &auxInt);

	fprintf(pFile, "%d; iteracoes.; %d,%d; ms.;\n", iterPon.iteracao, (int) auxInt, (int) round(auxFrac*1000000));

	fclose(pFile);
}

void printDoubleToFile(std::ofstream& outfile, double numero) {
	//double auxFrac, auxInt;
	//auxFrac = modf(numero, &auxInt);
	//outfile << (int)auxInt << ',' << (int)round(auxFrac * 1000000) << ';';
	std::string aux = to_string(numero);
	aux[aux.find('.')] = ',';
	outfile << aux << ';';
}

// para benchmarks.h

// rotinas de impress�o para os valores de benchmark coletados
void benchmarksPrint(iterativo& iterPon) {
	//std::string aux("out_" + met2str(global::metodo) + "-" + global::arq_entrada.substr(9) + ".csv");

	double accJacobiano = 0, accCudaMemcpy = 0, accJacobianoStencil = 0, accCalcPQ = 0, accSistemaLinear = 0, accInitGPU = 0/*, accInitLib = 0*/;
	bool /*flgProcessoIterativo*/flgIterativo = 0, flgCudaMemcpy = 0, flgGeral = 0, flgAdmitancia = 0, flgJacobiano = 0, flgJacobianoStencil = 0, flgCalcPQ = 0, flgSistemaLinear = 0, flgInitGPU = 0, flgFluxo = 0, flgInitLib = 0;

	// acumula valores coletados no processo iterativo
	for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		switch (std::get<0>(i)) {
		case benchmarkType::geral:
			if (flgGeral) {
				throw - 1;
			}
			flgGeral = 1;
			break;
		case benchmarkType::admitancia:
			if (flgAdmitancia) {
				throw - 1;
			}
			flgAdmitancia = 1;
			break;
		case benchmarkType::init_GPU:
			accInitGPU += std::get<2>(i);
			flgInitGPU = 1;
			break;
		case benchmarkType::init_lib:
			//accInitLib += std::get<2>(i);
			flgInitLib = 1;
			break;
		case benchmarkType::calcPQ:
			accCalcPQ += std::get<2>(i);
			flgCalcPQ = 1;
			break;
		case benchmarkType::jacobiano:
			accJacobiano += std::get<2>(i);
			flgJacobiano = 1;
			break;
		case benchmarkType::jacobianoStencil_fill:
		case benchmarkType::jacobianoStencil_build:
		case benchmarkType::jacobianoStencil_rebuild:
			accJacobianoStencil += std::get<2>(i);
			flgJacobianoStencil = 1;
			break;
		case benchmarkType::sistemaLinear:
			accSistemaLinear += std::get<2>(i);
			flgSistemaLinear = 1;
			break;
		case benchmarkType::processoIterativo:
			if (flgIterativo) {
				throw - 1;
			}
			flgIterativo = 1;
			break;
		case benchmarkType::fluxo:
			if (flgFluxo) {
				throw - 1;
			}
			flgFluxo = 1;
			break;
		case benchmarkType::cudaMemcpy:
			flgCudaMemcpy = 1;
			accCudaMemcpy += std::get<2>(i);
			break;
		}
	}

	printf("\n Relatorio de execucao:\n\n");
	if (flgGeral) {
		printf("*Geral***********************************\n\n" AVISO_ERRO_TEMPO_TOTAL);
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::geral)
			{
				printf(" Tempo total de execucao: %f ms.\n Numero de iteracoes: %d.\n\n", std::get<2>(i), iterPon.iteracao);
			}
		}
	}
	if (flgAdmitancia) {
		printf("*Ybus************************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::admitancia)
			{
				printf(" Tempo total: %f ms.\n\n", std::get<2>(i));
			}
		}
	}
	if (flgInitGPU) {
		printf("*Init GPU********************************\n\n");
		
		bool first = 1;

		printf(" acao                        -  tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::init_GPU)
			{
				if (first) {
					printf(" alocacao de memoria inicial - %14f ms (alem do inicio da comunicacao com a GPU)\n", std::get<2>(i));
					first = 0;
				}
				else {
					printf(" alocacao de varia. esparsas - %14f ms\n", std::get<2>(i));
				}
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n\n", accInitGPU);
	}
	if (flgCudaMemcpy) {
		printf("*Transferencias de memoria***************\n\n");

		int cnt = 0;

		printf(" acao                              -  tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::cudaMemcpy)
			{
				cnt++;
				printf(" %d-esima transferencia de memoria - %14f ms\n", cnt, std::get<2>(i));
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n\n", accInitGPU);
	}
	if (flgInitLib) {
		printf("*Inicializacao das Bibliotecas***********\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::init_lib)
			{
				printf(" Tempo total: %f ms.\n\n", std::get<2>(i));
			}
		}
	}
	if (flgCalcPQ) {
		printf("*CalcPQ**********************************\n\n");
		printf("*Tempo nao inclui transferencias de dados para o host!\n\n");
		printf(" iteracao - tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::calcPQ)
			{
				printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n", accCalcPQ);
		printf(" Tempo medio: %f ms.\n\n", accCalcPQ / (double)iterPon.iteracao);
	}
	if (flgJacobiano) {
		printf("*Jacobiano*******************************\n\n");
		printf(" iteracao - tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::jacobiano)
			{
				printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n", accJacobiano);
		printf(" Tempo medio: %f ms.\n\n", accJacobiano / (double)iterPon.iteracao);
	}
	if (flgJacobianoStencil) {
		printf("*Jacobiano Stencil***********************\n");
		printf("  acao  - iteracao - tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::jacobianoStencil_fill)
			{
				printf(" preenc - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
			else if (std::get<0>(i) == benchmarkType::jacobianoStencil_build) {
				printf(" constr - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
			else if (std::get<0>(i) == benchmarkType::jacobianoStencil_rebuild) {
				printf(" recons - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n", accJacobianoStencil);
		printf(" Tempo medio: %f ms (por iteracao).\n\n", accJacobianoStencil / (double)iterPon.iteracao);
	}
	if (flgSistemaLinear) {
		printf("*Sistema Linear**************************\n\n");
		printf(" iteracao - tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::sistemaLinear)
			{
				printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
		}
		printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n", accSistemaLinear);
		printf(" Tempo medio: %f ms.\n\n", accSistemaLinear / (double)iterPon.iteracao);
	}
	if (flgIterativo) {
		printf("*Processo Iterativo**********************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::processoIterativo)
			{
				printf(" Tempo total: %f ms.\n\n", std::get<2>(i));
			}
		}
	}
	if (flgFluxo) {
		printf("*Fluxo***********************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::fluxo)
			{
				printf(" Tempo total: %f ms.\n\n", std::get<2>(i));
				//break;
			}
		}

		printf("*****************************************\n");
	}
}

// rotinas de impress�o para os valores de benchmark coletados
void benchmarksPrintFile(iterativo& iterPon) {
	std::string aux;
	std::string ptFl = DOUBLE_MODE ? "double_" : "float_";
	if (global::CPUsolverFlg) {
		aux = std::string("out_" + ptFl + std::string("hibridoC") + "-" + /*global::arq_entrada*/ global::arq_entrada.substr(0, global::arq_entrada.find_last_of('.')) + ".csv");
	}
	else {
		aux = std::string("out_" + ptFl + met2str(global::metodo) + "-" + /*global::arq_entrada*/ global::arq_entrada.substr(0, global::arq_entrada.find_last_of('.')) + ".csv");
	}

	unsigned noDoTeste = 1;
	// verifica se o cabe�alho ja foi escrito
	std::ifstream infile(aux);
	if (infile.good()) {
		// conta n�mero de linhas
		// new lines will be skipped unless we stop it from happening:    
		infile.unsetf(std::ios_base::skipws);

		// conta linhas:
		noDoTeste = std::count(
			std::istream_iterator<char>(infile),
			std::istream_iterator<char>(),
			'\n');
		infile.close();
		noDoTeste--;
	}
	else {
		// arquivo percisa ser criado e cabecalho inserido
		std::ofstream outfile(aux);
		outfile << "Execucao;Aviso;ms total;no iteracoes;ms Ybus;ms Init GPU;;;ms Inicializacao;ms calcPQ;;ms Jacobiano;;ms Sistema Linear;;ms Processo;ms Fluxo;\n;;;;;total; Aloca��o Vari�veis iniciais; Aloca��o de vari�veis esparsas; Biblioteca; total; medio por iteracao; total; medio por iteracao; total; medio por iteracao; Iterativo;\n";
		outfile.close();
	}
	std::ofstream outfile(aux, std::ios_base::app);
	outfile << noDoTeste << ";"; // no da execucao
	//outfile.close();



	//std::ofstream outfile;
	//outfile.open(aux, std::ios_base::app); // append
	//double auxFrac, auxInt;
	//auxFrac = modf(duracao/*.count()*/, &auxInt);
	//outfile << iterPon.iteracao << "; iteracoes.; " << (int)auxInt << ',' << (int)round(auxFrac * 1000000) << "; ms.;\n";
	/*printDoubleToFile(outfile, 10.5);*/

	double accJacobiano = 0, accJacobianoStencil = 0, accCalcPQ = 0, accSistemaLinear = 0, accInitGPU = 0/*, accInitLib = 0*/;
	bool /*flgProcessoIterativo*/flgIterativo = 0, flgGeral = 0, flgAdmitancia = 0, flgJacobiano = 0, flgJacobianoStencil = 0, flgCalcPQ = 0, flgSistemaLinear = 0, flgInitGPU = 0, flgFluxo = 0, flgInitLib = 0;

	// acumula valores coletados no processo iterativo
	for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		switch (std::get<0>(i)) {
		case benchmarkType::geral:
			if (flgGeral) {
				throw - 1;
			}
			flgGeral = 1;
			break;
		case benchmarkType::admitancia:
			if (flgAdmitancia) {
				throw - 1;
			}
			flgAdmitancia = 1;
			break;
		case benchmarkType::init_GPU:
			accInitGPU += std::get<2>(i);
			flgInitGPU = 1;
			break;
		case benchmarkType::init_lib:
			//accInitLib += std::get<2>(i);
			flgInitLib = 1;
			break;
		case benchmarkType::calcPQ:
			accCalcPQ += std::get<2>(i);
			flgCalcPQ = 1;
			break;
		case benchmarkType::jacobiano:
			accJacobiano += std::get<2>(i);
			flgJacobiano = 1;
			break;
		case benchmarkType::jacobianoStencil_fill:
		case benchmarkType::jacobianoStencil_build:
		case benchmarkType::jacobianoStencil_rebuild:
			accJacobianoStencil += std::get<2>(i);
			flgJacobianoStencil = 1;
			break;
		case benchmarkType::sistemaLinear:
			accSistemaLinear += std::get<2>(i);
			flgSistemaLinear = 1;
			break;
		case benchmarkType::processoIterativo:
			if (flgIterativo) {
				throw - 1;
			}
			flgIterativo = 1;
			break;
		case benchmarkType::fluxo:
			if (flgFluxo) {
				throw - 1;
			}
			flgFluxo = 1;
			break;
		}
	}

	//printf("\n Relatorio de execucao:\n\n");
	if (flgGeral) {
		//printf("*Geral***********************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::geral)
			{
				outfile << AVISO_ERRO_TEMPO_TOTAL_FILE << "; ";
				printDoubleToFile(outfile, std::get<2>(i)); // Tempo total de execucao
				outfile << iterPon.iteracao << ';';         // Numero de iteracoes
			}
		}
	}
	else {
		outfile << ';' << ';';
	}
	if (flgAdmitancia) {
		//printf("*Ybus************************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::admitancia)
			{
				printDoubleToFile(outfile, std::get<2>(i)); // Tempo de execucao
			}
		}
	}
	else {
		outfile << ';';
	}
	if (flgInitGPU) {
		printDoubleToFile(outfile, accInitGPU);
		
		bool first = 1;
		int aux = 1;
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::init_GPU)
			{
				if (first) {
					printDoubleToFile(outfile, std::get<2>(i));
					first = 0;
					aux++;
				}
				else {
					printDoubleToFile(outfile, std::get<2>(i));
					first = 1;
					aux++;
				}
			}
		}
		if (/*!first*/ aux <= 2) {
			// preenche celula em branco
			outfile << ' ' << ';';
		}
	}
	else {
		outfile << ';' << ';' << ';';
	}
	if (flgInitLib) {
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::init_lib)
			{
				printDoubleToFile(outfile, std::get<2>(i));
			}
		}
	}
	else {
		outfile << ';' << ';';
	}
	if (flgCalcPQ) {
		//printf("*CalcPQ**********************************\n\n");
		//printf(" iteracao - tempo de execucao\n");
		//for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		//	if (std::get<0>(i) == calcPQ)
		//	{
		//		printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//}
		/*printf("-----------------------------\n");
		printf(" Tempo total: %f ms.\n", accCalcPQ);
		printf(" Tempo medio: %f ms.\n\n", accCalcPQ / (double)iterPon.iteracao);*/
		printDoubleToFile(outfile, accCalcPQ);                            // Tempo de total
		printDoubleToFile(outfile, accCalcPQ / (double)iterPon.iteracao); // Tempo medio
	}
	else {
		outfile << ';' << ';';
	}
	if (flgJacobiano) {
		//printf("*Jacobiano*******************************\n\n");
		//printf(" iteracao - tempo de execucao\n");
		//for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		//	if (std::get<0>(i) == jacobiano)
		//	{
		//		printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//}
		//printf("-----------------------------\n");
		//printf(" Tempo total: %f ms.\n", accJacobiano);
		//printf(" Tempo medio: %f ms.\n\n", accJacobiano / (double)iterPon.iteracao);
		printDoubleToFile(outfile, accJacobiano);                            // Tempo de total
		printDoubleToFile(outfile, accJacobiano / (double)iterPon.iteracao); // Tempo medio
	}

	if (flgJacobianoStencil) {
		//printf("*Jacobiano Stencil***********************n\n");
		//printf("  acao  - iteracao - tempo de execucao\n");
		//for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		//	if (std::get<0>(i) == jacobianoStencil_fill)
		//	{
		//		printf(" preenc - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//	else if (std::get<0>(i) == jacobianoStencil_build) {
		//		printf(" constr - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//	else if (std::get<0>(i) == jacobianoStencil_rebuild) {
		//		printf(" recons - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//}
		//printf("-----------------------------\n");
		//printf(" Tempo total: %f ms.\n", accJacobianoStencil);
		//printf(" Tempo medio: %f ms (por iteracao).\n\n", accJacobianoStencil / (double)iterPon.iteracao);
		printDoubleToFile(outfile, accJacobianoStencil);                            // Tempo de total
		printDoubleToFile(outfile, accJacobianoStencil / (double)iterPon.iteracao); // Tempo medio
	}
	else if (!flgJacobiano) {
		// se nenhum jacobiano...
		outfile << ';' << ';';
	}
	if (flgSistemaLinear) {
		//printf("*Sistema Linear**************************\n\n");
		//printf(" iteracao - tempo de execucao\n");
		//for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		//	if (std::get<0>(i) == sistemaLinear)
		//	{
		//		printf(" %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
		//	}
		//}
		//printf("-----------------------------\n");
		//printf(" Tempo total: %f ms.\n", accSistemaLinear);
		//printf(" Tempo medio: %f ms.\n\n", accSistemaLinear / (double)iterPon.iteracao);
		printDoubleToFile(outfile, accSistemaLinear);                            // Tempo de total
		printDoubleToFile(outfile, accSistemaLinear / (double)iterPon.iteracao); // Tempo medio
	}
	else {
		outfile << ';' << ';';
	}
	if (flgIterativo) {
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::processoIterativo)
			{
				printDoubleToFile(outfile, std::get<2>(i));
			}
		}
	}
	else {
		outfile << ';';
	}
	if (flgFluxo) {
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == benchmarkType::fluxo)
			{
				printDoubleToFile(outfile, std::get<2>(i));
				//break;
			}
		}
	}
	else {
		outfile << ';';
	}

	outfile << '\n';

	outfile.close();

	//std::ofstream outfile;
	//outfile.open(aux, std::ios_base::app); // append
	//double auxFrac, auxInt;
	//auxFrac = modf(duracao/*.count()*/, &auxInt);
	//outfile << iterPon.iteracao << "; iteracoes.; " << (int)auxInt << ',' << (int)round(auxFrac * 1000000) << "; ms.;\n";
}