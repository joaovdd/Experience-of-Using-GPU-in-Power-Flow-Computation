#pragma once
#include <iostream>
#include <fstream>

#include "../idx.h"
#include "../estruturas.h"
#include "../global.h"
#include "id2i.h"

#include "config.h"

struct mpBusType
{
    std::vector<int> BUS_I;      
    std::vector<int> BUS_TYPE;   
    std::vector<double> PD;      
    std::vector<double> QD;      
    std::vector<double> GS;      
    std::vector<double> BS;      
    std::vector<int> BUS_AREA;   
    std::vector<double> VM;      
    std::vector<double> VA;      
    std::vector<double> BASE_KV; 
    std::vector<int> ZONE;       
    std::vector<double> VMAX;    
    std::vector<double> VMIN;    
    std::vector<double> LAM_P;   
    std::vector<double> LAM_Q;   
    std::vector<double> MU_VMAX; 
    std::vector<double> MU_VMIN; 
};

struct mpGenType
{
    std::vector<int> GEN_BUS;       
    std::vector<double> PG;         
    std::vector<double> QG;         
    std::vector<double> QMAX;       
    std::vector<double> QMIN;       
    std::vector<double> VG;         
    std::vector<double> MBASE;      
    std::vector<double> GEN_STATUS; 
    std::vector<double> PMAX;       
    std::vector<double> PMIN;       
    std::vector<double> PC1;        
    std::vector<double> PC2;        
    std::vector<double> QC1MIN;     
    std::vector<double> QC1MAX;     
    std::vector<double> QC2MIN;     
    std::vector<double> QC2MAX;     
    std::vector<double> RAMP_AGC;   
    std::vector<double> RAMP_10;    
    std::vector<double> RAMP_30;    
    std::vector<double> RAMP_Q;     
    std::vector<double> APF;        
    std::vector<double> MU_PMAX;    
    std::vector<double> MU_PMIN;    
    std::vector<double> MU_QMAX;    
    std::vector<double> MU_QMIN;    
};

struct mpBranchType
{
    std::vector<double> F_BUS;     
    std::vector<double> T_BUS;     
    std::vector<double> BR_R;      
    std::vector<double> BR_X;      
    std::vector<double> BR_B;      
    std::vector<double> RATE_A;    
    std::vector<double> RATE_B;    
    std::vector<double> RATE_C;    
    std::vector<double> TAP;       
    std::vector<double> SHIFT;     
    std::vector<double> BR_STATUS; 
    std::vector<double> ANGMIN;    
    std::vector<double> ANGMAX;    
    std::vector<double> PF;        
    std::vector<double> QF;        
    std::vector<double> PT;        
    std::vector<double> QT;        
    std::vector<double> MU_SF;     
    std::vector<double> MU_ST;     
    std::vector<double> MU_ANGMIN; 
    std::vector<double> MU_ANGMAX; 
};

struct matPowerDataType
{
    mpBusType bus;
    mpGenType gen;
    mpBranchType branch;
    double baseMVA;
    double version;
};

matPowerDataType lerArquivoMatPOWER(std::string filePath);
matPowerDataType lerArquivoMatPOWEReficiente(const std::string& filePath);

matPowerDataType lerMatPowerEAlocarMemoria(std::string mFile, sistemaType& sistema, barraType& barra, ramoType& ramo, iterativoType& iterativo);

bool lerTamanhosMatPOWER(matPowerDataType data, sistemaType& sistema);

bool storeMatPOWER(sistemaType& sistema, barraType& barra, ramoType& ramo, matPowerDataType& mpData);