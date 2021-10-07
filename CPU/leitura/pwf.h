
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
	std::vector<int> numero;
	std::vector<char> operacao;
	std::vector<char> estado;
	std::vector<char> tipo;

	std::vector <std::string> grupoDeBaseDeTensao;
	std::vector<std::string> nome;

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
	std::vector<float> resistencia; 
	std::vector<float> reatancia; 
	std::vector<float> susceptancia; 
	std::vector<float> tap; 
	std::vector<float> tapMinimo;
	std::vector<float> tapMaximo;
	std::vector<float> defasagem; 
	std::vector<int> barraControlada; 
	std::vector<float> capacidadeNormal; 
	std::vector<float> capacidadeEmEmergencia; 
	std::vector<int> numeroDeTaps; 
	std::vector<float> capacidadeDeEquipamento; 
	std::vector<std::vector<char>> agregadores;
};

struct DCTEtype
{
	float TEPA = 0.1f;   
	float TEPR = 0.1f;   
	float TLPR = 0.1f;   
	float TLVC = 0.5f;   
	float TLTC = 0.01f;  
	float TETP = 5.0f;   
	float TBPA = 5.0f;   
	float TSFR = 0.01f;  
	float TUDC = 0.001f; 
	float TADC = 0.01f;  
	float BASE = 100.0f; 
	float DASE = 100.0f; 
	float ZMAX = 500.0f; 
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
	float TSBZ = 0.01f;  
	float TSBA = 1.0f;   
	float PGER = 0.01f;  
	float VDVN = 5.0f;   
	float VDVM = 30.0f;  
	float ASTP = 40.0f;  
	float VSTP = 5.0f;   
	float CSTP = 5.0f;   
	float VFLD = 70.f;   
	int HIST;
	float ZMIN = 0.0001f;
	int PDIT = 10;
	int ICIT = 30;
	float FDIV = 2.0f;   
	int DMAX = 5;   
	float ICMN = 0.05f;  
	float VART = 5.0f;   
	int TSTP = 33;   
	float TSDC = 0.01f;   
	float ASDC = 1.f;     
	float ICMV = 0.5f;   
	int APAS = 90;   
	int CPAR = 70;   
	float VAVT = 2.0f;   
	float VAVF = 5.0f;   
	float VMVF = 15.0f;  
	float VPVT = 2.0f;   
	float VPVF = 5.0f;   
	float VPMF = 10.0f;  
	float VSVF = 20.0f;  
	float VINF = 1.0f;   
	float VSUP = 1.0f;   
	float TLSI = 0.0f;   
	float NDIR = 20.f;   
	float STIR = 1.f;   
	float STTR = 5.f;    
	float TRPT = 100.f;  
	float BFPO = 1.f;    
	float LFPO = 0.1f;  
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
	std::vector<float> valorMinimo; 
	std::vector<float> valorMaximo; 
	std::vector<float> valorInicial; 
	std::vector<char> modoDeControle; 
	std::vector<float> valorEspecificado; 
	std::vector<int> extremidadeDeMedicao; 
	std::vector<int> numeroDeEstagios; 
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

void parseDCTEline(DCTEtype& DCTE, std::string& str_pwf, unsigned long long& inicio); 

void parseDCSCline(DCSCtype& DCSC, std::string& str_pwf, unsigned long long& inicio);

void parseDGBTline(DGBTtype& DGBT, std::string& str_pwf, unsigned long long& inicio);

void parse(DBARtype& estrutura, std::string& str_pwf);

void parse(DLINtype& estrutura, std::string& str_pwf);

void parse(DCTEtype& estrutura, std::string& str_pwf);

void parse(DGBTtype& estrutura, std::string& str_pwf);

void parse(DCSCtype& estrutura, std::string& str_pwf);

void getPWF(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo);

void lerPWFEAlocarMemoria(std::string pwfFile, sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo);