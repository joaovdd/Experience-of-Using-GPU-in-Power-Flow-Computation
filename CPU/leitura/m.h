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
    std::vector<int> BUS_I;      //   bus number (positive integer)
    std::vector<int> BUS_TYPE;   //   bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)
    std::vector<double> PD;      //   real power demand (MW)
    std::vector<double> QD;      //   reactive power demand (MVAr)
    std::vector<double> GS;      //   shunt conductance (MW demanded at V = 1.0 p.u.)
    std::vector<double> BS;      //   shunt susceptance (MVAr injected at V = 1.0 p.u.)
    std::vector<int> BUS_AREA;   //   area number (positive integer)
    std::vector<double> VM;      //   voltage magnitude (p.u.)
    std::vector<double> VA;      //   voltage angle (degrees)
    std::vector<double> BASE_KV; //   base voltage (kV)
    std::vector<int> ZONE;       //   loss zone (positive integer)
    std::vector<double> VMAX;    // † maximum voltage magnitude (p.u.)
    std::vector<double> VMIN;    // † minimum voltage magnitude (p.u.)
    std::vector<double> LAM_P;   // † Lagrange multiplier on real power mismatch (u/MW)
    std::vector<double> LAM_Q;   // † Lagrange multiplier on reactive power mismatch (u/MVAr)
    std::vector<double> MU_VMAX; // † Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
    std::vector<double> MU_VMIN; // † Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)
};

struct mpGenType
{
    std::vector<int> GEN_BUS;       //   1 bus number
    std::vector<double> PG;         //   2 real power output (MW)
    std::vector<double> QG;         //   3 reactive power output (MVAr)
    std::vector<double> QMAX;       //   4 maximum reactive power output (MVAr)
    std::vector<double> QMIN;       //   5 minimum reactive power output (MVAr)
    std::vector<double> VG;         // ‡ 6 voltage magnitude setpoint (p.u.)
    std::vector<double> MBASE;      //   7 total MVA base of machine, defaults to baseMVA
    std::vector<double> GEN_STATUS; //   8 machine status, > 0 = machine in-service, <= 0 = machine out-of-service
    std::vector<double> PMAX;       //   9 maximum real power output (MW)
    std::vector<double> PMIN;       //   10 minimum real power output (MW)
    std::vector<double> PC1;        // * 11 lower real power output of PQ capability curve (MW)
    std::vector<double> PC2;        // * 12 upper real power output of PQ capability curve (MW)
    std::vector<double> QC1MIN;     // * 13 minimum reactive power output at PC1 (MVAr)
    std::vector<double> QC1MAX;     // * 14 maximum reactive power output at PC1 (MVAr)
    std::vector<double> QC2MIN;     // * 15 minimum reactive power output at PC2 (MVAr)
    std::vector<double> QC2MAX;     // * 16 maximum reactive power output at PC2 (MVAr)
    std::vector<double> RAMP_AGC;   // * 17 ramp rate for load following/AGC (MW/min)
    std::vector<double> RAMP_10;    // * 18 ramp rate for 10 minute reserves (MW)
    std::vector<double> RAMP_30;    // * 19 ramp rate for 30 minute reserves (MW)
    std::vector<double> RAMP_Q;     // * 20 ramp rate for reactive power (2 sec timescale) (MVAr/min)
    std::vector<double> APF;        // * 21 area participation factor
    std::vector<double> MU_PMAX;    // † 22 Kuhn-Tucker multiplier on upper Pg limit (u/MW)
    std::vector<double> MU_PMIN;    // † 23 Kuhn-Tucker multiplier on lower Pg limit (u/MW)
    std::vector<double> MU_QMAX;    // † 24 Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
    std::vector<double> MU_QMIN;    // † 25 Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)
};

struct mpBranchType
{
    std::vector<double> F_BUS;     //   1 "from" bus number
    std::vector<double> T_BUS;     //   2 "to" bus number
    std::vector<double> BR_R;      //   3 resistance (p.u.)
    std::vector<double> BR_X;      //   4 reactance (p.u.)
    std::vector<double> BR_B;      //   5 total line charging susceptance (p.u.)
    std::vector<double> RATE_A;    // * 6 MVA rating A (long term rating), set to 0 for unlimited
    std::vector<double> RATE_B;    // * 7 MVA rating B (short term rating), set to 0 for unlimited
    std::vector<double> RATE_C;    // * 8 MVA rating C (emergency rating), set to 0 for unlimited
    std::vector<double> TAP;       //   9 transformer o nominal turns ratio, if non-zero (taps at \from" bus, impedance at \to" bus, i.e. if r = x = b = 0, tap = jVf j jVtj ; tap = 0 used to indicate transmission line rather than transformer, i.e. mathematically equivalent to transformer with tap = 1)
    std::vector<double> SHIFT;     //   10 transformer phase shift angle (degrees), positive ) delay
    std::vector<double> BR_STATUS; //   11 initial branch status, 1 = in-service, 0 = out-of-service
    std::vector<double> ANGMIN;    // † 12 minimum angle dierence, f 􀀀 t (degrees)
    std::vector<double> ANGMAX;    // † 13 maximum angle dierence, f 􀀀 t (degrees)
    std::vector<double> PF;        // ‡ 14 real power injected at \from" bus end (MW)
    std::vector<double> QF;        // ‡ 15 reactive power injected at \from" bus end (MVAr)
    std::vector<double> PT;        // ‡ 16 real power injected at \to" bus end (MW)
    std::vector<double> QT;        // ‡ 17 reactive power injected at \to" bus end (MVAr)
    std::vector<double> MU_SF;     // § 18 Kuhn-Tucker multiplier on MVA limit at \from" bus (u/MVA)
    std::vector<double> MU_ST;     // § 19 Kuhn-Tucker multiplier on MVA limit at \to" bus (u/MVA)
    std::vector<double> MU_ANGMIN; // § 20 Kuhn-Tucker multiplier lower angle dierence limit (u/degree)
    std::vector<double> MU_ANGMAX; // § 21 Kuhn-Tucker multiplier upper angle dierence limit (u/degree)
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

// bool readMatPOWER(matPowerDataType data, sistemaType& sistema, barraType& barra, ramoType& ramo);

bool storeMatPOWER(sistemaType& sistema, barraType& barra, ramoType& ramo, matPowerDataType& mpData);