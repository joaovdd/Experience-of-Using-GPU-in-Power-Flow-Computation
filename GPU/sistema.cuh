#pragma once
#include "sistema.h"
#include "dim.h"

#include <assert.h>
#include "cusparse.h"
#include "cusolverSp.h"

void inline d_showVecf(complex_type* show, int dim);
void d_showVecf(int* show, int dim);

sistema d_initSistema(sistema& h_sistema, sistema* d_sistema) {
	sistema aux; 

	aux.barraVO = 0;
	aux.nPV = 0;
	aux.nPQ = 0;
	aux.nB = 0;
	aux.nL = 0;
	aux.baseMVA = 0;

	switch (global::metodo) {
	case metodo::denso:

		aux.csrRowPtrY = nullptr;
		aux.csrColIndY = nullptr;
		aux.cooRowIndY = nullptr;
		checkCudaErrors(cudaMalloc(&(aux.Y), h_sistema.nB * h_sistema.nB * sizeof(complex_type))); 	
		checkCudaErrors(cudaMemset(aux.Y, 0, h_sistema.nB * h_sistema.nB * sizeof(complex_type)));
		break;
	case metodo::esparso:
	case metodo::hibridoB:
		aux.Y = nullptr;
		checkCudaErrors(cudaMalloc(&(aux.csrColIndY), h_sistema.nnzY * sizeof(int)));
		checkCudaErrors(cudaMalloc(&(aux.csrRowPtrY), (h_sistema.nB + 1) * sizeof(int)));
		checkCudaErrors(cudaMalloc(&(aux.cooRowIndY), h_sistema.nnzY * sizeof(int)));

		checkCudaErrors(cudaMalloc(&(aux.spYval), h_sistema.nnzY * sizeof(complex_type)));
		break;
	}

	checkCudaErrors(cudaMalloc(&(aux.barrasPV), h_sistema.nPV * sizeof(int))); 					
	checkCudaErrors(cudaMemset(aux.barrasPV, 0, h_sistema.nPV * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.barrasPQ), h_sistema.nPQ * sizeof(int))); 					
	checkCudaErrors(cudaMemset(aux.barrasPQ, 0, h_sistema.nPQ * sizeof(int)));

	checkCudaErrors(cudaMemcpy(d_sistema, &aux, sizeof(sistema), cudaMemcpyHostToDevice));

	return aux;
}

void d_finSistema(sistema* d_sistema, sistema aux) {
	if (aux.barrasPV) { checkCudaErrors(cudaFree(aux.barrasPV)); }
	if (aux.barrasPQ) { checkCudaErrors(cudaFree(aux.barrasPQ)); }
	if (aux.Y != nullptr) { checkCudaErrors(cudaFree(aux.Y)); }

	if (aux.spYval) { cudaFree(aux.spYval); }
	if (aux.csrRowPtrY) { cudaFree(aux.csrRowPtrY); }
	if (aux.csrColIndY) { cudaFree(aux.csrColIndY); }
	if (aux.cooRowIndY) { cudaFree(aux.cooRowIndY); }

	if (d_sistema) { checkCudaErrors(cudaFree(d_sistema)); }
}

void sistemacpyH2D(sistema& h_sistema, sistema* d_sistema, sistema& aux) {
	aux.barraVO = h_sistema.barraVO;
	aux.nPV = h_sistema.nPV;
	aux.nPQ = h_sistema.nPQ;
	aux.nB = h_sistema.nB;
	aux.nL = h_sistema.nL;

	checkCudaErrors(cudaMemcpy(aux.barrasPV, h_sistema.barrasPV, h_sistema.nPV * sizeof(int), cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(aux.barrasPQ, h_sistema.barrasPQ, h_sistema.nPQ * sizeof(int), cudaMemcpyHostToDevice));

	switch (global::metodoDeCalculoDeYbus) {
	case metodoDeCalculoDeYbus::dnCPU:

		checkCudaErrors(cudaMemcpy(aux.Y, h_sistema.Y, h_sistema.nB * h_sistema.nB * sizeof(complex_type), cudaMemcpyHostToDevice));
		break;

	case metodoDeCalculoDeYbus::spCPU:
		checkCudaErrors(cudaMemcpy(aux.csrColIndY, h_sistema.csrColIndY, h_sistema.nnzY * sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(aux.csrRowPtrY, h_sistema.csrRowPtrY, (h_sistema.nB + 1) * sizeof(int), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(aux.cooRowIndY, h_sistema.cooRowIndY, h_sistema.nnzY * sizeof(int), cudaMemcpyHostToDevice)); 

		checkCudaErrors(cudaMemcpy(aux.spYval, h_sistema.spYval, h_sistema.nnzY * sizeof(complex_type), cudaMemcpyHostToDevice));

		aux.nnzY = h_sistema.nnzY;

		if (global::verbose_mode) {
			printf("sistPon.nnzY =\n");
			printf("%d\n", aux.nnzY);
			printf("sistPon.spYval =\n");
			d_showVecf(aux.spYval, aux.nnzY);
			printf("sistPon.csrRowPtrY =\n");
			d_showVecf(aux.csrRowPtrY, aux.nB + 1);
			printf("sistPon.csrColIndY =\n");
			d_showVecf(aux.csrColIndY, aux.nnzY);
			printf("sistPon.cooRowIndY =\n");
			d_showVecf(aux.cooRowIndY, aux.nnzY);
		}
		break;
	default:
		cout << "[sistemacpyH2D] ERRO metodo invalido de Ybus.\n" << endl;
	}

	checkCudaErrors(cudaMemcpy(d_sistema, &aux, sizeof(sistema), cudaMemcpyHostToDevice));
}

void sistemacpyD2H(sistema &h_sistema, sistema *d_sistema, sistema &aux){
	 
	checkCudaErrors(cudaMemcpy(&(aux), d_sistema , sizeof(sistema), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_sistema.Y, aux.Y, h_sistema.nB*h_sistema.nB * sizeof(complex_type), cudaMemcpyDeviceToHost));

	cudaDeviceSynchronize();

	h_sistema.barraVO = aux.barraVO;
	h_sistema.nPV     = aux.nPV;
	h_sistema.nPQ     = aux.nPQ;
	h_sistema.nB      = aux.nB;
	h_sistema.nL      = aux.nL;

}

__global__ void printAddress(barra* d_barra){
	printf("d_initBarra: d_barra->V     = %p\n", d_barra->V);
	printf("d_initBarra: d_barra->theta = %p\n", d_barra->theta);
	printf("d_initBarra: d_barra->Pliq  = %p\n", d_barra->Pliq);
	printf("d_initBarra: d_barra->Qliq  = %p\n", d_barra->Qliq);
	printf("d_initBarra: d_barra->Pload = %p\n", d_barra->Pload);
	printf("d_initBarra: d_barra->Qload = %p\n", d_barra->Qload);
	printf("d_initBarra: d_barra->Pg    = %p\n", d_barra->Pg);
	printf("d_initBarra: d_barra->Qg    = %p\n", d_barra->Qg);
	printf("d_initBarra: d_barra->Vbase = %p\n", d_barra->Vbase);
	printf("d_initBarra: d_barra->gsh   = %p\n", d_barra->gsh);
	printf("d_initBarra: d_barra->bsh   = %p\n\n", d_barra->bsh);
}

__global__ void printAddress(ramo* d_ramo){
	printf("d_initRamo: d_ramo->z    = %p\n", d_ramo->z);
	printf("d_initRamo: d_ramo->bsh  = %p\n", d_ramo->bsh);
	printf("d_initRamo: d_ramo->Pdp  = %p\n", d_ramo->Pdp);
	printf("d_initRamo: d_ramo->Pdp  = %p\n", d_ramo->Ppd);
	printf("d_initRamo: d_ramo->Qdp  = %p\n", d_ramo->Qdp);
	printf("d_initRamo: d_ramo->Qdp  = %p\n", d_ramo->Qpd);
	printf("d_initRamo: d_ramo->tap  = %p\n", d_ramo->tap);
	printf("d_initRamo: d_ramo->de   = %p\n", d_ramo->de);
	printf("d_initRamo: d_ramo->para = %p\n\n", d_ramo->para);
}

barra d_initBarra(sistema &h_sistema, barra* d_barra){
	barra aux; 

	checkCudaErrors(cudaMalloc(&(aux.V), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.V, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.theta), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.theta, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Pliq), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Pliq, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Qliq), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Qliq, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Pload), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Pload, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Qload), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Qload, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Pg), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Pg, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Qg), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Qg, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Vbase), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Vbase, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.gsh), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.gsh, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.bsh), h_sistema.nB * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.bsh, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMemcpy(d_barra, &aux, sizeof(barra), cudaMemcpyHostToDevice));

	return aux;
}

void d_finBarra(barra* d_barra, barra aux){
	if (aux.V) { checkCudaErrors(cudaFree(aux.V)); }
	if (aux.theta) { checkCudaErrors(cudaFree(aux.theta)); }
	if (aux.Pliq) { checkCudaErrors(cudaFree(aux.Pliq)); }
	if (aux.Qliq) { checkCudaErrors(cudaFree(aux.Qliq)); }
	if (aux.Pload) { checkCudaErrors(cudaFree(aux.Pload)); }
	if (aux.Qload) { checkCudaErrors(cudaFree(aux.Qload)); }
	if (aux.Pg) { checkCudaErrors(cudaFree(aux.Pg)); }
	if (aux.Qg) { checkCudaErrors(cudaFree(aux.Qg)); }
	if (aux.Vbase) { checkCudaErrors(cudaFree(aux.Vbase)); }
	if (aux.gsh) { checkCudaErrors(cudaFree(aux.gsh)); }
	if (aux.bsh) { checkCudaErrors(cudaFree(aux.bsh)); }

	if (d_barra) { checkCudaErrors(cudaFree(d_barra)); }

}

void barracpyH2D(sistema &h_sistema, barra &h_barra, barra* d_barra, barra aux){ 
	checkCudaErrors(cudaMemcpy(aux.V, h_barra.V, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(aux.theta, h_barra.theta, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Pliq,  h_barra.Pliq,  h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Qliq,  h_barra.Qliq,  h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Pload, h_barra.Pload, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Qload, h_barra.Qload, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Pg,    h_barra.Pg,    h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Qg,    h_barra.Qg,    h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Vbase, h_barra.Vbase, h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.gsh,   h_barra.gsh,   h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.bsh,   h_barra.bsh,   h_sistema.nB * sizeof(float_type), cudaMemcpyHostToDevice));
}

void barracpyD2H(sistema &h_sistema, barra &h_barra, barra* d_barra, barra aux){ 
	checkCudaErrors(cudaMemcpy(aux.V,     h_barra.V,     h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.theta, h_barra.theta, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Pliq,  h_barra.Pliq,  h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Qliq,  h_barra.Qliq,  h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Pload, h_barra.Pload, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Qload, h_barra.Qload, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Pg,    h_barra.Pg,    h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Qg,    h_barra.Qg,    h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Vbase, h_barra.Vbase, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.gsh,   h_barra.gsh,   h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.bsh,   h_barra.bsh,   h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
}

ramo d_initRamo(sistema &h_sistema, ramo& h_ramo, ramo* d_ramo){
	ramo aux; 

	checkCudaErrors(cudaMalloc(&(aux.z), h_sistema.nL * sizeof(complex_type)));	
	checkCudaErrors(cudaMemset(aux.z, 0, h_sistema.nL * sizeof(complex_type)));	

	checkCudaErrors(cudaMalloc(&(aux.bsh), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.bsh, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.tap), h_sistema.nL * sizeof(complex_type)));
	checkCudaErrors(cudaMemset(aux.tap, 0, h_sistema.nL * sizeof(complex_type)));

	checkCudaErrors(cudaMalloc(&(aux.phi), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.phi, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.d_csrColIndPhi), h_ramo.eigen_phi->innerSize() * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.d_csrColIndPhi, 0, h_ramo.eigen_phi->innerSize() * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.d_csrRowPtrPhi), h_ramo.eigen_phi->outerSize() * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.d_csrRowPtrPhi, 0, h_ramo.eigen_phi->outerSize() * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.phiVal), h_ramo.eigen_phi->nonZeros() * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.phiVal, 0, h_ramo.eigen_phi->nonZeros() * sizeof(float_type)));

	aux.nnzPhi = h_ramo.eigen_phi->nonZeros();

	checkCudaErrors(cudaMalloc(&(aux.Pdp), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Pdp, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Ppd), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Ppd, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Qdp), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Qdp, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.Qpd), h_sistema.nL * sizeof(float_type)));
	checkCudaErrors(cudaMemset(aux.Qpd, 0, h_sistema.nL * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.de), h_sistema.nL * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.de, 0, h_sistema.nL * sizeof(int)));

	checkCudaErrors(cudaMalloc(&(aux.para), h_sistema.nL * sizeof(int)));
	checkCudaErrors(cudaMemset(aux.para, 0, h_sistema.nL * sizeof(int)));

	checkCudaErrors(cudaMemcpy(d_ramo, &aux, sizeof(ramo), cudaMemcpyHostToDevice)); 

	return aux;
}

void d_finRamo(ramo *d_ramo, ramo aux){
	if (aux.z)    { checkCudaErrors(cudaFree(aux.z));   }
	if (aux.bsh)  { checkCudaErrors(cudaFree(aux.bsh)); }
	if (aux.tap)  { checkCudaErrors(cudaFree(aux.tap)); }
	if (aux.phi)  { checkCudaErrors(cudaFree(aux.phi)); }
	if (aux.d_csrColIndPhi)  { checkCudaErrors(cudaFree(aux.d_csrColIndPhi)); }
	if (aux.d_csrRowPtrPhi)  { checkCudaErrors(cudaFree(aux.d_csrRowPtrPhi)); }
	if (aux.phiVal)  { checkCudaErrors(cudaFree(aux.phiVal)); }

	if (aux.Pdp)  { checkCudaErrors(cudaFree(aux.Pdp)); }
	if (aux.Ppd)  { checkCudaErrors(cudaFree(aux.Ppd)); }
	if (aux.Qdp)  { checkCudaErrors(cudaFree(aux.Qdp)); }
	if (aux.Qpd)  { checkCudaErrors(cudaFree(aux.Qpd)); }
	if (aux.de)   { checkCudaErrors(cudaFree(aux.de));  }
	if (aux.para) { checkCudaErrors(cudaFree(aux.para));}

	if (d_ramo)   { checkCudaErrors(cudaFree(d_ramo)); }

}

void ramocpyH2D(sistema &h_sistema, ramo &h_ramo, ramo *d_ramo, ramo aux){ 

	checkCudaErrors(cudaMemcpy(aux.z,   h_ramo.z,   h_sistema.nL * sizeof(complex_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.bsh, h_ramo.bsh, h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.tap, h_ramo.tap, h_sistema.nL * sizeof(complex_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.phi, h_ramo.phi, h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(aux.phi, h_ramo.phi, h_sistema.nL * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.d_csrColIndPhi, h_ramo.eigen_phi->innerIndexPtr(), h_ramo.eigen_phi->innerSize() * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.d_csrRowPtrPhi, h_ramo.eigen_phi->outerIndexPtr(), h_ramo.eigen_phi->outerSize() * sizeof(float_type), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.phiVal, h_ramo.eigen_phi->valuePtr(), h_ramo.eigen_phi->nonZeros() * sizeof(float_type), cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(aux.Ppd,   h_ramo.Ppd,   h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Pdp,   h_ramo.Pdp,   h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Qpd,   h_ramo.Qpd,   h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.Qdp,   h_ramo.Qdp,   h_sistema.nL * sizeof(float_type),          cudaMemcpyHostToDevice));

	checkCudaErrors(cudaMemcpy(aux.de,   h_ramo.de,   h_sistema.nL * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(aux.para, h_ramo.para, h_sistema.nL * sizeof(int), cudaMemcpyHostToDevice));
}

void ramocpyD2H(sistema &h_sistema, ramo &h_ramo, ramo *d_ramo, ramo aux){ 

	checkCudaErrors(cudaMemcpy(h_ramo.z,   aux.z,   h_sistema.nL * sizeof(complex_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.bsh, aux.bsh, h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.tap, aux.tap, h_sistema.nL * sizeof(complex_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.phi, aux.phi, h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaMemcpy(aux.Ppd,   h_ramo.Ppd,   h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Pdp,   h_ramo.Pdp,   h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Qpd,   h_ramo.Qpd,   h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(aux.Qdp,   h_ramo.Qdp,   h_sistema.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaMemcpy(h_ramo.de,   aux.de,   h_sistema.nL * sizeof(int), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.para, aux.para, h_sistema.nL * sizeof(int), cudaMemcpyDeviceToHost));
}

iterativo d_initIter(sistema &h_sistema, iterativo *d_iterativo){
	iterativo aux; 

	checkCudaErrors(cudaMalloc(&(aux.Pcalc), h_sistema.nB * sizeof(float_type))); 
	checkCudaErrors(cudaMemset(aux.Pcalc, 0, h_sistema.nB * sizeof(float_type))); 

	checkCudaErrors(cudaMalloc(&(aux.Qcalc), h_sistema.nB * sizeof(float_type))); 
	checkCudaErrors(cudaMemset(aux.Qcalc, 0, h_sistema.nB * sizeof(float_type))); 

	if (global::metodo == metodo::denso || (global::metodo == metodo::hibridoB && !global::CPUsolverFlg)) {
		cudaMalloc(&(aux.Jlim), (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * sizeof(float_type)); 
		cudaMemset(aux.Jlim, 0, (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * (h_sistema.nPV + h_sistema.nPV + h_sistema.nPQ + h_sistema.nPQ) * sizeof(float_type)); 
	}
	else {
		aux.Jlim = nullptr;
	}

	aux.iteracao = 0;
	aux.noMax = 0;

	checkCudaErrors(cudaMalloc(&(aux.barrasPQlim), (h_sistema.nPV + h_sistema.nPQ) * sizeof(int))); 
	checkCudaErrors(cudaMalloc(&(aux.barrasPVlim), h_sistema.nPV * sizeof(int))); 

	checkCudaErrors(cudaMalloc(&(aux.gLim), (h_sistema.nB - 1) * 2 * sizeof(float_type))); 
	checkCudaErrors(cudaMemset(aux.gLim, 0, (h_sistema.nB - 1) * 2 * sizeof(float_type)));

	checkCudaErrors(cudaMalloc(&(aux.QliqLim), h_sistema.nB * sizeof(float_type))); 
	checkCudaErrors(cudaMemset(aux.QliqLim, 0, h_sistema.nB * sizeof(float_type)));

	checkCudaErrors(cudaMemcpy(d_iterativo, &aux, sizeof(iterativo), cudaMemcpyHostToDevice));						 

	return aux;
}

void d_finIter(iterativo* d_iterativo, iterativo aux){
	if (aux.Pcalc)   { checkCudaErrors(cudaFree(aux.Pcalc)); }
	if (aux.Qcalc)   { checkCudaErrors(cudaFree(aux.Qcalc)); }
	if (aux.Jlim != nullptr) { checkCudaErrors(cudaFree(aux.Jlim)); }

	if (aux.barrasPQlim) { checkCudaErrors(cudaFree(aux.barrasPQlim)); }
	if (aux.barrasPVlim) { checkCudaErrors(cudaFree(aux.barrasPVlim)); }
	if (aux.gLim) { checkCudaErrors(cudaFree(aux.gLim)); }
	if (aux.QliqLim) { checkCudaErrors(cudaFree(aux.QliqLim)); }

	if (d_iterativo) { checkCudaErrors(cudaFree(d_iterativo)); }
}

void itercpyH2D(sistema &h_sistema, iterativo &h_iterativo, iterativo *d_iterativo, iterativo &aux){
	aux.iteracao = h_iterativo.iteracao;
	aux.noMax    = h_iterativo.noMax;

	checkCudaErrors(cudaMemcpy(d_iterativo, &aux, sizeof(iterativo), cudaMemcpyHostToDevice));

}

void itercpyD2H(sistema &h_sistema, iterativo &h_iterativo, iterativo *d_iterativo, iterativo &aux){
	checkCudaErrors(cudaMemcpy(&aux,  d_iterativo,  sizeof(iterativo), cudaMemcpyDeviceToHost));

	h_iterativo.iteracao = aux.iteracao;
	h_iterativo.noMax    = aux.noMax;

	checkCudaErrors(cudaMemcpy(h_iterativo.Pcalc, aux.Pcalc, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_iterativo.Qcalc, aux.Qcalc, h_sistema.nB * sizeof(float_type), cudaMemcpyDeviceToHost));

}

void printbarra(barra* d_barra) {
	printf("printBarra: &(d_barra->V)     = %p\n", &(d_barra->V));
	printf("printBarra: &(d_barra->theta) = %p\n", &(d_barra->theta));
	printf("printBarra: &(d_barra->Pliq)  = %p\n", &(d_barra->Pliq));
	printf("printBarra: &(d_barra->Qliq)  = %p\n", &(d_barra->Qliq));
	printf("printBarra: &(d_barra->Pload) = %p\n", &(d_barra->Pload));
	printf("printBarra: &(d_barra->Qload) = %p\n", &(d_barra->Qload));
	printf("printBarra: &(d_barra->Pg)    = %p\n", &(d_barra->Pg));
	printf("printBarra: &(d_barra->Qg)    = %p\n", &(d_barra->Qg));
	printf("printBarra: &(d_barra->Vbase) = %p\n", &(d_barra->Vbase));
	printf("printBarra: &(d_barra->gsh)   = %p\n", &(d_barra->gsh));
	printf("printBarra: &(d_barra->bsh)   = %p\n\n", &(d_barra->bsh));
}

cudaDeviceProp initGPU(){
	int deviceCount = 0;
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

	if (error_id != cudaSuccess) {
		printf("cudaGetDeviceCount returned %d\n-> %s\n",
			static_cast<int>(error_id), cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		exit(EXIT_FAILURE);
	}

	printf("\n\n");

	if (deviceCount == 0) {
		printf("There are no available device(s) that support CUDA\n");
	} else {
		printf("Detected %d CUDA Capable device(s)\n", deviceCount);
	}

    int devID = gpuGetMaxGflopsDeviceId();
    checkCudaErrors(cudaSetDevice(devID));
	cudaDeviceProp deviceProp;
    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));

    if(!global::laconic_mode){
		printf("GPU Device %d: \"%s\" with compute capability %d.%d\n\n", devID, deviceProp.name, deviceProp.major, deviceProp.minor);

		printf("  Warp size:                                    %d\n",
			   deviceProp.warpSize);
		printf("  Maximum number of threads per multiprocessor: %d\n",
			   deviceProp.maxThreadsPerMultiProcessor);
		printf("  Maximum number of threads per block:          %d\n",
			   deviceProp.maxThreadsPerBlock);
		printf("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
			   deviceProp.multiProcessorCount,
			   _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
			   _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
				   deviceProp.multiProcessorCount);
    }

	return deviceProp;

}

void d_criarYesparso(sistema& sistPon){
	cusparseHandle_t cusparseHandle = NULL;
	cusparseMatDescr_t descr = NULL;
	int nnzTotal = 0;
	int* d_nnzPerRowColumn;
	checkCudaErrors(cudaMalloc(&(d_nnzPerRowColumn), (sistPon.nB) * sizeof(int)));

	cusparseStatus_t cuSPARSEstatus = cusparseCreate(&cusparseHandle);
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	cuSPARSEstatus = cusparseCreateMatDescr(&descr);
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	cuSPARSEstatus = cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL); 
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	cuSPARSEstatus = cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ONE);
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	cuSPARSEstatus = _cusparseNnz(cusparseHandle, CUSPARSE_DIRECTION_ROW,
    			                  sistPon.nB, sistPon.nB, descr, sistPon.Y, sistPon.nB,
								  d_nnzPerRowColumn, &nnzTotal);

	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	int* d_csrRowPtr;
	int* d_csrColInd;
	complex_type* d_csrVal;

	checkCudaErrors(cudaMalloc((void **)&d_csrRowPtr, sizeof(int)*(sistPon.nB+1)));
	checkCudaErrors(cudaMalloc((void **)&d_csrColInd, sizeof(int)*(nnzTotal)));
	checkCudaErrors(cudaMalloc((void **)&d_csrVal, sizeof(complex_type)*(nnzTotal)));

	cuSPARSEstatus = _cusparseDense2csr(cusparseHandle, sistPon.nB, sistPon.nB, descr, sistPon.Y, sistPon.nB,
										d_nnzPerRowColumn, d_csrVal, d_csrRowPtr, d_csrColInd);
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);									

	int* d_cooRowInd;
	checkCudaErrors(cudaMalloc((void **)&d_cooRowInd, sizeof(int)*(nnzTotal)));

	cuSPARSEstatus = cusparseXcsr2coo(cusparseHandle, d_csrRowPtr, nnzTotal, sistPon.nB,	d_cooRowInd, CUSPARSE_INDEX_BASE_ONE);
	assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);

	if (cusparseHandle != NULL) {
		cuSPARSEstatus = cusparseDestroy(cusparseHandle);
		assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
	}
	if (descr != NULL) {
		cuSPARSEstatus = cusparseDestroyMatDescr(descr);
		assert(CUSPARSE_STATUS_SUCCESS == cuSPARSEstatus);
	}

	sistPon.spYval = d_csrVal;
	sistPon.csrRowPtrY = d_csrRowPtr;
	sistPon.csrColIndY = d_csrColInd;
	sistPon.cooRowIndY = d_cooRowInd;
	sistPon.nnzY = nnzTotal;

	if(d_nnzPerRowColumn) {checkCudaErrors(cudaFree(d_nnzPerRowColumn)); }
}