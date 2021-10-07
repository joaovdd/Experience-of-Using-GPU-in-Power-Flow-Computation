#include "impressao.h"

void impressao2(sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo) {
	printf("\n\n                                 FLUXOS DE POTENCIA\n");

	printf("        BUS      V         ANG        Pli        Qli        Pge        Qge      TO        PIJ          QIJ\n");
	for (int i = 0; i < sistema.nB; i++)
	{
		
		printf("      %4d%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f  ", barra.id[i], barra.V[i], barra.theta[i], iterativo.Pcalc[i], iterativo.Qcalc[i], iterativo.Pcalc[i] + barra.Pload[i], iterativo.Qcalc[i] + barra.Qload[i]);
		
		int j;
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

void impressao(sistemaType &sistema, barraType &barra, ramoType &ramo, iterativoType &iterativo) {
	printf("\n\n                                 FLUXOS DE POTENCIA\n");
            
	printf("        BUS      V         ANG         P          Q        TO       PIJ          QIJ\n");
	for (int i = 0; i < sistema.nB; i++)
	{
		
		printf("      %4d%11.5f%11.5f%11.5f%11.5f  ", barra.id[i], barra.V[i], barra.theta[i], iterativo.Pcalc[i], iterativo.Qcalc[i]);
		
		int j;
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
				printf("                                                        %4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.para[j])], ramo.Pdp[j], ramo.Qdp[j]);
			}
			else if (ramo.para[j] == (i + 1)) {
				printf("                                                        %4d  %11.5f  %11.5f\n", barra.id[IDX1F(ramo.de[j])], ramo.Ppd[j], ramo.Qpd[j]);
			}
		}

	}
}

std::string met2str(metodoType met) {
	switch (met)
	{
	case denso:
		return "denso";
		break;
	case esparso:
		return "esparso";
		break;
	case esparsoSimples:
		return "esparsosimples";
		break;
	default:
		return "nda";
		break;
	}
}

void benchmarkModePrint(iterativoType& iterPon, double& duracao) {
	std::string aux("out_" + met2str(global::metodo) + "-" + global::arq_entrada.substr(9) + ".csv");

	std::ofstream outfile;
	outfile.open(aux, std::ios_base::app); 
	double auxFrac, auxInt;
	auxFrac = modf(duracao, &auxInt);
	outfile << iterPon.iteracao << "; iteracoes.; " << (int)auxInt << ',' << (int)round(auxFrac * 1000000) << "; ms.;\n";
}

void benchmarksPrint(iterativoType& iterPon) {
	double accJacobiano = 0, accJacobianoStencil = 0, accCalcPQ = 0, accSistemaLinear = 0;
	bool flgFluxo = 0, flgIterativo = 0, flgGeral = 0, flgAdmitancia = 0, flgJacobiano = 0, flgJacobianoStencil = 0, flgCalcPQ = 0, flgSistemaLinear = 0;

	for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		switch (std::get<0>(i)) {
		case geral:
			if (flgGeral) {
				throw - 1;
			}
			flgGeral = 1;
			break;
		case admitancia:
			if (flgAdmitancia) {
				throw - 1;
			}
			flgAdmitancia = 1;
			break;
		case calcPQ:
			accCalcPQ += std::get<2>(i);
			flgCalcPQ = 1;
			break;
		case jacobiano:
			accJacobiano += std::get<2>(i);
			flgJacobiano = 1;
			break;
		case jacobianoStencil_fill:
		case jacobianoStencil_build:
		case jacobianoStencil_rebuild:
			accJacobianoStencil += std::get<2>(i);
			flgJacobianoStencil = 1;
			break;
		case sistemaLinear:
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

	printf("\n Relatorio de execucao:\n\n");
	if (flgGeral) {
		printf("*Geral***********************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == geral)
			{
				printf(" Tempo total de execucao: %f ms.\n Numero de iteracoes: %d.\n\n", std::get<2>(i), iterPon.iteracao);
			}
		}
	}
	if (flgAdmitancia) {
		printf("*Ybus************************************\n\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == admitancia)
			{
				printf(" Tempo total: %f ms.\n\n", std::get<2>(i));
			}
		}
	}
	if (flgCalcPQ) {
		printf("*CalcPQ**********************************\n\n");
		printf(" iteracao - tempo de execucao\n");
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == calcPQ)
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
			if (std::get<0>(i) == jacobiano)
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
			if (std::get<0>(i) == jacobianoStencil_fill)
			{
				printf(" preenc - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
			else if (std::get<0>(i) == jacobianoStencil_build) {
				printf(" constr - %8d - %14f ms\n", std::get<1>(i), std::get<2>(i));
			}
			else if (std::get<0>(i) == jacobianoStencil_rebuild) {
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
			if (std::get<0>(i) == sistemaLinear)
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
				
			}
		}

		printf("*****************************************\n");
	}
}

void printDoubleToFile(std::ofstream& outfile, double numero) {
	double auxFrac, auxInt;
	auxFrac = modf(numero, &auxInt);
	outfile << (int)auxInt << ',' << (int)round(auxFrac * 1000000) << ';';
}

void benchmarksPrintFile(iterativoType& iterPon) {
	std::string ptFl = DOUBLE_MODE ? "double_" : "float_";

	std::string aux("out_" + ptFl + met2str(global::metodo) + "-" +  global::arq_entrada.substr(global::arq_entrada.find_last_of('/') + 1, global::arq_entrada.find_last_of('.') - 1 - global::arq_entrada.find_last_of('/')) + ".csv"); 

	unsigned noDoTeste = 1;

	std::ifstream infile(aux);
	if (infile.good()) {
		
		infile.unsetf(std::ios_base::skipws);

		noDoTeste = std::count(
			std::istream_iterator<char>(infile),
			std::istream_iterator<char>(),
			'\n');
		infile.close();
		noDoTeste--;
	}
	else {
		
		std::ofstream outfile(aux);
		outfile << "Execucao;ms total;no iteracoes;ms Ybus;ms calcPQ;;ms Jacobiano;;ms Sistema Linear;;ms Processo;ms Fluxo;\n;;;; total; medio por iteracao; total; medio por iteracao; total; medio por iteracao;\n";
		outfile.close();
	}
	std::ofstream outfile(aux, std::ios_base::app);
	outfile << noDoTeste << ";"; 

	double accJacobiano = 0, accJacobianoStencil = 0, accCalcPQ = 0, accSistemaLinear = 0;
	bool flgFluxo = 0, flgIterativo = 0, flgGeral = 0, flgAdmitancia = 0, flgJacobiano = 0, flgJacobianoStencil = 0, flgCalcPQ = 0, flgSistemaLinear = 0;

	for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
		switch (std::get<0>(i)) {
		case geral:
			if (flgGeral) {
				throw - 1;
			}
			flgGeral = 1;
			break;
		case admitancia:
			if (flgAdmitancia) {
				throw - 1;
			}
			flgAdmitancia = 1;
			break;
		case calcPQ:
			accCalcPQ += std::get<2>(i);
			flgCalcPQ = 1;
			break;
		case jacobiano:
			accJacobiano += std::get<2>(i);
			flgJacobiano = 1;
			break;
		case jacobianoStencil_fill:
		case jacobianoStencil_build:
		case jacobianoStencil_rebuild:
			accJacobianoStencil += std::get<2>(i);
			flgJacobianoStencil = 1;
			break;
		case sistemaLinear:
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

	if (flgGeral) {
		
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == geral)
			{
				printDoubleToFile(outfile, std::get<2>(i)); 
				outfile << iterPon.iteracao << ';';         
			}
		}
	}
	else {
		outfile << ';' << ';';
	}
	if (flgAdmitancia) {
		
		for (std::tuple<benchmarkType, int, double> i : global::tracker.benchmarkTable) {
			if (std::get<0>(i) == admitancia)
			{
				printDoubleToFile(outfile, std::get<2>(i)); 
			}
		}
	}
	else {
		outfile << ';';
	}
	if (flgCalcPQ) {
		
		printDoubleToFile(outfile, accCalcPQ);                            
		printDoubleToFile(outfile, accCalcPQ / (double)iterPon.iteracao); 
	}
	else {
		outfile << ';' << ';';
	}
	if (flgJacobiano) {
		
		printDoubleToFile(outfile, accJacobiano);                            
		printDoubleToFile(outfile, accJacobiano / (double)iterPon.iteracao); 
	}

	if (flgJacobianoStencil) {
		
		printDoubleToFile(outfile, accJacobianoStencil);                            
		printDoubleToFile(outfile, accJacobianoStencil / (double)iterPon.iteracao); 
	}
	else if (!flgJacobiano) {
		
		outfile << ';' << ';';
	}
	if (flgSistemaLinear) {
		
		printDoubleToFile(outfile, accSistemaLinear);                            
		printDoubleToFile(outfile, accSistemaLinear / (double)iterPon.iteracao); 
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
				
			}
		}
	}
	else {
		outfile << ';';
	}

	outfile << '\n';

	outfile.close();

}