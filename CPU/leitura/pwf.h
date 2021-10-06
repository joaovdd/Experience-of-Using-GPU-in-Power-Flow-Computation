// Baseado no Manual do Programa de Análise de Redes (Anarede) V10.02.03
//
// (C) 2020 João Victor Daher Daibes


// Manual: file:///C:/Cepel/Anarede/V100203/Manuais/Manual-Anarede.pdf
// DBAR -> p82 (2.13.4)
// DLIN -> p153(2.54.4)

// DAGR -> p73 (2.9.4)				Leitura dos dados de agregadores genéricos.
// DOPC -> p167(2.62.4)				Leitura e modificação dos dados das Opções de Controle de Execução padrão.
// DCTE -> p111(2.29.1)				Leitura e modificação dos dados de constantes utilizadas no programa.
// DCSC -> p109(2.28.4)				Leitura dos dados de CSC (Compensador Série Controlável).
// DBSH // FBAN -> p87  (2.15.4)	Leitura  dos  dados  de  bancos  de  capacitores  e/ou  reatores  individualizados  conectados  a  barras  CA  ou  linhas  de  transmissão. 
// DSHL -> p176(2.69.4)				Leitura dos dados de dispositivos shunt de circuito CA.
// DCER -> p100(2.22.4)				Leitura dos dados de compensador estático de reativos. 
// DCTR -> p123(2.31.4)				Leitura  dos  dados  complementares  de  transformadores.
// DGBT -> p135(2.40.4)				Leitura dos dados de grupos de base de tensão de barras CA. 

#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include <array>
#include <map>

#include "../estrutura/sistema.h"
#include "../estrutura/barra.h"
#include "../estrutura/ramo.h"
#include "../estrutura/iterativo.h"
#include "id2i.h"
#include "../idx.h"
#include "../global.h"

struct DBARtype
{
	//int numBarra = 0;
	std::vector<int> numero;
	std::vector<char> operacao;
	std::vector<char> estado;
	std::vector<char> tipo;
	//std::vector<std::pair <char, char>> grupoDeBaseDeTensao;
	std::vector <std::string> grupoDeBaseDeTensao;
	std::vector<std::string> nome;
	//std::vector<std::pair <char, char>> grupoDeLimiteDeTensao;
	std::vector<std::string> grupoDeLimiteDeTensao;
	std::vector<float> tensao;
	std::vector<float> angulo;
	std::vector<float> geracaoAtiva;
	std::vector<float> geracaoReativa;
	std::vector<float> geracaoReativaMinima;
	std::vector<float> geracaoReativaMaxima;
	std::vector<int> barraControlada;
	std::vector<float> cargaAtiva;
	std::vector<float> cargaReativa;
	std::vector<float> capacitorReator;
	std::vector<int> area;
	std::vector<float> tensaoParaDefinicaoDeCarga;
	std::vector<float> modoDeVisualizacao;
	std::vector<std::vector<char>> agregadores;
};

struct DLINtype
{
	std::vector<int> daBarra;
	std::vector<char> aberturaDaBarra;
	std::vector<char> operacao;
	std::vector<char> aberturaParaABarra;
	std::vector<int> paraBarra;
	std::vector<int> circuito;
	std::vector<char> estado;
	std::vector<char> proprietario;
	std::vector<float> resistencia; // %
	std::vector<float> reatancia; // %
	std::vector<float> susceptancia; // %
	std::vector<float> tap; // p.u.
	std::vector<float> tapMinimo;
	std::vector<float> tapMaximo;
	std::vector<float> defasagem; // graus
	std::vector<int> barraControlada; // número
	std::vector<float> capacidadeNormal; // MVA
	std::vector<float> capacidadeEmEmergencia; // MVA
	std::vector<int> numeroDeTaps; // número
	std::vector<float> capacidadeDeEquipamento; // Capacidade de carregamento do equipamento com menor capacidade de carregamento conectado ao circuito. 
	std::vector<std::vector<char>> agregadores;
};

struct DCTEtype
{
	float TEPA = 0.1f;   // MW
	float TEPR = 0.1f;   // MVAr
	float TLPR = 0.1f;   // MVAr
	float TLVC = 0.5f;   // %
	float TLTC = 0.01f;  // %
	float TETP = 5.0f;   // MVAr
	float TBPA = 5.0f;   // MVAr
	float TSFR = 0.01f;  // %
	float TUDC = 0.001f; // %
	float TADC = 0.01f;  // %
	float BASE = 100.0f; // MVA
	float DASE = 100.0f; // MW
	float ZMAX = 500.0f; // %
	int ACIT = 30;
	int LPIT = 50;
	int LFLP = 10;
	int LFIT = 10;
	int DCIT = 10;
	int VSIT = 10;
	int LCRT = 23;
	int LPRT = 60;
	int LFCV = 1;
	int TPST = 2 * TEPR;
	int QLST = 4 * TEPR;
	int EXST = 4 * TEPA;
	float TLPP = 1.0;
	float TSBZ = 0.01f;  // MW
	float TSBA = 1.0f;   // %
	float PGER = 0.01f;  // MW
	float VDVN = 5.0f;   // MW
	float VDVM = 30.0f;  // %
	float ASTP = 40.0f;  // rd
	float VSTP = 5.0f;   // %
	float CSTP = 5.0f;   // %
	float VFLD = 70.f;   // %
	int HIST;
	float ZMIN = 0.0001f;// %
	int PDIT = 10;
	int ICIT = 30;
	float FDIV = 2.0f;   // 
	int DMAX = 5;   // 
	float ICMN = 0.05f;  // %
	float VART = 5.0f;   // %
	int TSTP = 33;   // 
	float TSDC = 0.01f;   // 
	float ASDC = 1.f;     // º
	float ICMV = 0.5f;   // %
	int APAS = 90;   // %
	int CPAR = 70;   // %
	float VAVT = 2.0f;   // %
	float VAVF = 5.0f;   // %
	float VMVF = 15.0f;  // MW
	float VPVT = 2.0f;   // %
	float VPVF = 5.0f;   // %
	float VPMF = 10.0f;  // MW
	float VSVF = 20.0f;  // %
	float VINF = 1.0f;   // 
	float VSUP = 1.0f;   // 
	float TLSI = 0.0f;   // 
	float NDIR = 20.f;   // 
	float STIR = 1.f;   // 
	float STTR = 5.f;    // %
	float TRPT = 100.f;  // %
	float BFPO = 1.f;    // Mvar
	float LFPO = 0.1f;  // MW
};

struct DCSCtype
{
	std::vector<int> daBarra;
	std::vector<char> operacao;
	std::vector<int> paraBarra;
	std::vector<int> circuito;
	std::vector<char> estado;
	std::vector<char> proprietario;
	std::vector<char> bypass;
	std::vector<float> valorMinimo; // %
	std::vector<float> valorMaximo; // %
	std::vector<float> valorInicial; // %
	std::vector<char> modoDeControle; // = 'X';
	std::vector<float> valorEspecificado; // Reatância do CSC, em %, se o modo de controle especificado é Reatância Constante (X). 
	std::vector<int> extremidadeDeMedicao; // Número da barra terminal do CSC na qual a potência ou a corrente é medida, como definido no campo Número do Código de Execução DBAR.
	std::vector<int> numeroDeEstagios; // Número de estágios do CSC discreto (TSSC - Thyristor Switched Series Capacitor). O valor default é para o CSC que opera de modo contínuo (TCSC - Thyristor Controlled Series Capacitor). 
	std::vector<int> agregador;
};

struct DGBTtype
{
	std::map<std::string, float> mapa;
};

bool readPWF(sistemaType& sistema, barraType& barra, ramoType& ramo, DCTEtype& DCTE, DBARtype& DBAR, DLINtype& DLIN, DGBTtype& DGBT, DCSCtype& DCSC);

std::string ifBlanckSetDefault(std::string str, std::string defaultValue);

float parseFloatPtoDecImp(std::string str, int casaUnitariaDoPontoDecimalImplicito);

void parseDBARline(DBARtype& DBAR, std::string& str_pwf, unsigned long long& inicio);

void parseDLINline(DLINtype& DLIN, std::string& str_pwf, unsigned long long& inicio);

void parseDCTEline(DCTEtype& DCTE, std::string& str_pwf, unsigned long long& inicio); // inacabado

void parseDCSCline(DCSCtype& DCSC, std::string& str_pwf, unsigned long long& inicio);

void parseDGBTline(DGBTtype& DGBT, std::string& str_pwf, unsigned long long& inicio);

void parse(DBARtype& estrutura, std::string& str_pwf);

void parse(DLINtype& estrutura, std::string& str_pwf);

void parse(DCTEtype& estrutura, std::string& str_pwf);

void parse(DGBTtype& estrutura, std::string& str_pwf);

void parse(DCSCtype& estrutura, std::string& str_pwf);

void getPWF(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo);

void lerPWFEAlocarMemoria(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo);