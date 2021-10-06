// numero de threads: deve ser potência de 2
// nl é o número de LTs conectadas à barra k (card[k])
int noThreadsPQ(int nl){
	int aux = 1;
	while (nl > aux){
		aux <<= 1;
		//printf("noThrads: %d\n", aux);
	}
	return aux;
}

// cardinalidade do conjunto K-1: a quantas barras está ligada cada barra
int* cardK(sistema &sistema, ramo &ramo){
	int* card = (int *)malloc(sistema.nB * sizeof(int));
	for (int i = 0; i < sistema.nB; i++){
		card[i] = 0;
	}
	for (int i = 0; i < sistema.nB; i++){
		card[IDX1F(ramo.de[i])]++;
		card[IDX1F(ramo.para[i])]++;
	}
	return card;
}

// retorna a posição do elemento de índices (lin, col) no vetor de valores da matriz
// armazenada da forma csr (overload para indexadores int)
__host__ __device__ int coeffPos(int lin, int col, int* csrRowPtrY, int* csrColIndY, int nnzY) {
	for (int i = csrRowPtrY[/*IDX1F(*/lin/*)*/]; i < csrRowPtrY[/*IDX1F(*/lin + 1/*)*/]; i++) {
		if (csrColIndY[i] == /*IDX1F(*/col/*)*/) {
			return i;
		}
	}
	return -1;
}

// retorna valor da defasagem angular entre as barras a e b
// pode se tornar mais eficiente (será?) ao se utilizar a função coeff (O(log(nnz_j)).
__device__ float_type phif(int a, int b, const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	auto aux = coeffPos(IDX1F(a), IDX1F(b), ramoPon.d_csrRowPtrPhi, ramoPon.d_csrColIndPhi, ramoPon.nnzPhi);
	if (aux == -1) {
		aux = coeffPos(IDX1F(b), IDX1F(a), ramoPon.d_csrRowPtrPhi, ramoPon.d_csrColIndPhi, ramoPon.nnzPhi);
		if (aux == -1) {
			return 0.;
		}
		else {
			return -ramoPon.phiVal[aux];
		}
	}
	else {
		return ramoPon.phiVal[aux];
	}
}

// retorna valor da defasagem angular entre as barras a e b
__host__ __device__ float_type phif(int a, int b, sistema* sistema, ramo* ramo) {
#ifdef __CUDA_ARCH__
	auto aux = coeffPos(IDX1F(a), IDX1F(b), ramo->d_csrRowPtrPhi, ramo->d_csrColIndPhi, ramo->nnzPhi);
	if (aux == -1) {
		aux = coeffPos(IDX1F(b), IDX1F(a), ramo->d_csrRowPtrPhi, ramo->d_csrColIndPhi, ramo->nnzPhi);
		if (aux == -1) {
			return 0.;
		}
		else {
			return -ramo->phiVal[aux];
}
	}
	else {
		return ramo->phiVal[aux];
	}
#else
	// return ramo->eigen_phi->coeff(IDX1F(a), IDX1F(b));
	// return ramo->eigen_phi->coeff(IDX1F(row), IDX1F(col));
	auto aux = coeffPos(IDX1F(a), IDX1F(b), ramo->eigen_phi->outerIndexPtr(), ramo->eigen_phi->innerIndexPtr(), ramo->eigen_phi->nonZeros());
	if (aux == -1) {
		aux = coeffPos(IDX1F(b), IDX1F(a), ramo->eigen_phi->outerIndexPtr(), ramo->eigen_phi->innerIndexPtr(), ramo->eigen_phi->nonZeros());
		if (aux == -1) {
			return 0.;
		}
		else {
			return -ramo->eigen_phi->valuePtr()[aux];
		}
	}
	else {
		return ramo->eigen_phi->valuePtr()[aux];
	}
#endif
}

//__device__ float_type phif(int a, int b, const sistema sistPon, const barra barraPon, const ramo ramoPon) {
//	for (int i = 0; i < sistPon.nL; i++) { // i percorre ramos
//		if (ramoPon.de[i] == a) { // se k está ligada à para[i] =>
//			if (ramoPon.para[i] == b) { // se k está ligada à para[i] =>
//				return (ramoPon.phi[i]);
//			}
//		}
//		else {
//			if (ramoPon.para[i] == a) { // se k está ligada à para[i] =>
//				if (ramoPon.de[i] == b) { // se k está ligada à para[i] =>
//					return (-(ramoPon.phi[i]));
//				}
//			}
//		}
//	}
//	return(0.);
//}
//
//__host__ __device__ float_type phif(int a, int b, sistema* sistema, ramo* ramo) {
//	for (int i = 0; i < sistema->nL; i++) { // i percorre ramos
//		if (ramo->de[i] == a) { // se k está ligada à para[i] =>
//			if (ramo->para[i] == b) { // se k está ligada à para[i] =>
//				return (ramo->phi[i]);
//			}
//		}
//		else {
//			if (ramo->para[i] == a) { // se k está ligada à para[i] =>
//				if (ramo->de[i] == b) { // se k está ligada à para[i] =>
//					return (-(ramo->phi[i]));
//				}
//			}
//		}
//	}
//	return(0.);
//}

__host__ __device__ float_type bshf(int a, int b, const sistema &sistPon, const ramo &ramoPon) {
	for (int i = 0; i < sistPon.nL; i++) { // i percorre ramos
		if (ramoPon.de[i] == a) { // se k está ligada à para[i] =>
			if (ramoPon.para[i] == b) { // se k está ligada à para[i] =>
				return (ramoPon.bsh[i]);
			}
		}
		else {
			if (ramoPon.para[i] == a) { // se k está ligada à para[i] =>
				if (ramoPon.de[i] == b) { // se k está ligada à para[i] =>
					return (ramoPon.bsh[i]);
				}
			}
		}
	}
	return(0.);
}

// // debugado! WORKING
// __global__ void calcPeficiente(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon){
// 	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
//  	if(idx <= sistPon.nB){
//  	//calculo de P_idx	
//  		float_type aux = 0;
// 		float_type acc = 0;
// 		//extern __shared__ int s[]; //martiz Y

// 		// idx itera sistPon.cooRowIndY
// 		// i é a coluna // i   itera sistPon.csrColIndY
// 		for(int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++){
// 			aux = phif(idx, sistPon.csrColIndY[IDX1F(i)], d_sistema, d_ramo); // defasagem do transformador
// 			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
// 			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
// 			// printf("***\n");
// 			aux += barraPon.theta[IDX1F(sistPon.cooRowIndY[IDX1F(i)])] - barraPon.theta[IDX1F(sistPon.csrColIndY[IDX1F(i)])]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
// 			aux *= barraPon.V[IDX1F(sistPon.csrColIndY[IDX1F(i)])];
// 			acc += aux;
// 			aux = 0;
// 		}
//  		acc *= barraPon.V[IDX1F(idx)];
// 		iterPon.Pcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
// 		//printf("P(%d) = %f\n", idx, acc);
//  		//printf("P(%d) = \n", idx);
//  	}
// }

// debugado!
__global__ void calcPeficiente(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
 	if(idx <= sistPon.nB){
 	//calculo de P_idx	
 		float_type aux = 0;
		float_type acc = 0;
		//extern __shared__ int s[]; //martiz Y

		// idx itera sistPon.cooRowIndY
		// i é a coluna // i   itera sistPon.csrColIndY
		for(int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++){
			//aux = phif(idx, sistPon.csrColIndY[IDX1F(i)], d_sistema, d_ramo); // defasagem do transformador
			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
			// printf("***\n");
			aux = barraPon.theta[IDX1F(sistPon.cooRowIndY[IDX1F(i)])] - barraPon.theta[IDX1F(sistPon.csrColIndY[IDX1F(i)])]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
			aux = _cuReal(sistPon.spYval[IDX1F(i)]) * cos(aux) + _cuImag(sistPon.spYval[IDX1F(i)]) * sin(aux);
			aux *= barraPon.V[IDX1F(sistPon.csrColIndY[IDX1F(i)])];
			acc += aux;
			aux = 0;
		}
 		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
		//printf("P(%d) = %f\n", idx, acc);
 		//printf("P(%d) = \n", idx);
 	}
}

__global__ void calcQeficiente(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
 	if(idx <= sistPon.nB){
 	//calculo de P_idx	
 		float_type aux = 0;
		float_type acc = 0;
		//extern __shared__ int s[]; //martiz Y

		// idx itera sistPon.cooRowIndY
		// i é a coluna // i   itera sistPon.csrColIndY
		for(int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++){
			//aux = phif(idx, sistPon.csrColIndY[IDX1F(i)], d_sistema, d_ramo); // defasagem do transformador
			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
			// printf("***\n");
			aux = barraPon.theta[IDX1F(sistPon.cooRowIndY[IDX1F(i)])] - barraPon.theta[IDX1F(sistPon.csrColIndY[IDX1F(i)])]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
			aux = _cuReal(sistPon.spYval[IDX1F(i)]) * sin(aux) - _cuImag(sistPon.spYval[IDX1F(i)]) * cos(aux);
			aux *= barraPon.V[IDX1F(sistPon.csrColIndY[IDX1F(i)])];
			acc += aux;
			aux = 0;
		}
 		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Qcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
		//printf("P(%d) = %f\n", idx, acc);
 		//printf("P(%d) = \n", idx);
 	}
}

__global__ void calcPeficiente_0based_sha(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;

	extern __shared__ int s_csrRowPtrY[];
	int* s_csrColIndY = s_csrRowPtrY + sistPon.nB;

		
	for (int i = threadIdx.x; i <= sistPon.nB; i+= blockDim.x)	{
		s_csrRowPtrY[i] = sistPon.csrRowPtrY[i];
	}

	for (int i = threadIdx.x; i < sistPon.nnzY; i += blockDim.x) {
		s_csrColIndY[i] = sistPon.csrColIndY[i];
	}
		

	if (idx <= sistPon.nB) {
		//calculo de P_idx	
		float_type aux = 0;
		float_type acc = 0;
		//extern __shared__ int s[]; //martiz Y

		// idx itera sistPon.cooRowIndY
		// i é a coluna // i   itera sistPon.csrColIndY
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			//aux = phif(idx, sistPon.csrColIndY[i] + 1, d_sistema, d_ramo); // defasagem do transformador
			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
			// printf("***\n");
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
			aux = _cuReal(sistPon.spYval[i]) * cos(aux) + _cuImag(sistPon.spYval[i]) * sin(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
		//printf("P(%d) = %f\n", idx, acc);
		//printf("P(%d) = \n", idx);
	}
}

__global__ void calcPeficiente_0based(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nB) {
		//calculo de P_idx	
		float_type aux = 0;
		float_type acc = 0;
		//extern __shared__ int s[]; //martiz Y

		// idx itera sistPon.cooRowIndY
		// i é a coluna // i   itera sistPon.csrColIndY
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			//aux = phif(idx, sistPon.csrColIndY[i] + 1, d_sistema, d_ramo); // defasagem do transformador
			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
			// printf("***\n");
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
			aux = _cuReal(sistPon.spYval[i]) * cos(aux) + _cuImag(sistPon.spYval[i]) * sin(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
		//printf("P(%d) = %f\n", idx, acc);
		//printf("P(%d) = \n", idx);
	}
}

__global__ void calcQeficiente_0based(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nB) {
		//calculo de P_idx	
		float_type aux = 0;
		float_type acc = 0;
		//extern __shared__ int s[]; //martiz Y

		// idx itera sistPon.cooRowIndY
		// i é a coluna // i   itera sistPon.csrColIndY
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			//aux = phif(idx, sistPon.csrColIndY[i] + 1, d_sistema, d_ramo); // defasagem do transformador
			// printf("sistPon.cooRowIndY(i) = %d\n", sistPon.cooRowIndY[IDX1F(i)]);
			// printf("sistPon.csrColIndY(i) = %d\n", sistPon.csrColIndY[IDX1F(i)]);
			// printf("***\n");
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; // theta_k para[i]
// 			aux = _cuReal(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * cos(aux) + _cuImag(sistPon.Y[IDX2F(sistPon.cooRowIndY[IDX1F(i)], sistPon.csrColIndY[IDX1F(i)], sistPon.nB)]) * sin(aux);
			aux = _cuReal(sistPon.spYval[i]) * sin(aux) - _cuImag(sistPon.spYval[i]) * cos(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Qcalc[IDX1F(idx)] = acc; // copia para a memoria global //ERRO!!!!!
		//printf("P(%d) = %f\n", idx, acc);
		//printf("P(%d) = \n", idx);
	}
}

__global__ void calcP_sha(const int k, sistema* d_sistema, barra* barra, ramo* ramo, float_type* d_Out) {
	extern __shared__ float_type s_data[];
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	int tid = threadIdx.x;
	if (idx <= d_sistema->nB) {
		//		printf("(%d)\n", idx);
		float_type aux = 0;

		//aux = phif(k, idx, d_sistema, ramo); // defasagem do transformador
//		printf("(%d); aux = %f\n", idx, aux);
		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; // theta_k para[i]
//		printf("(%d); aux = %f\n", idx, aux);
		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux) + _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux);
		//		printf("(%d); aux = %f\n", idx, aux);
		aux *= barra->V[IDX1F(idx)];
		//		printf("(%d); aux = %f\n", idx, aux);

		s_data[tid] = barra->V[IDX1F(k)] * aux;
	}
	else {
		s_data[tid] = 0.;
	}
	//printf("hello from (%d)\n", idx);

	//idx--;

	__syncthreads();

	//if (tid == 0) {
	//	for (int i = 0; i < 16; i++) {
	//		printf("s_data[%d] = %f\n", i, s_data[i]);
	//	}
	//}

	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			s_data[tid] += s_data[tid + s];
			//printf("idx = %d;\ts = %d\n", tid, s);
		}
		__syncthreads();
	}
	if (tid == 0) {
		d_Out[blockIdx.x] = s_data[0];
	}
}

//idx = 1:nB
__global__ void calcP(const int k, float_type* d_Pparcelas, sistema* d_sistema, barra* barra, ramo* ramo, int threads){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
	if(idx <= d_sistema->nB){
//		printf("(%d)\n", idx);
		float_type aux = 0;

		//aux = phif(k, idx, d_sistema, ramo); // defasagem do transformador
//		printf("(%d); aux = %f\n", idx, aux);
		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; // theta_k para[i]
//		printf("(%d); aux = %f\n", idx, aux);
		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux) + _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux);
//		printf("(%d); aux = %f\n", idx, aux);
		aux *= barra->V[IDX1F(idx)];
//		printf("(%d); aux = %f\n", idx, aux);

		d_Pparcelas[idx-1] = barra->V[IDX1F(k)]*aux;
	}
//	printf("hello from (%d)\n", idx);

	if (idx == 1) {
		for (int i = 0; i < 16; i++) {
			printf("s_data[%d] = %f\n", i, d_Pparcelas[i]);
		}
	}

	idx--;

		__syncthreads();

	for (int s = threads/2; s > 0; s >>= 1){
		if (idx < s){
			d_Pparcelas[idx] += d_Pparcelas[idx + s];
			//printf("idx = %d;\ts = %d\n", idx, s);
		}
		__syncthreads();
	}
}

//void d_calcPf(const int k, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, iterativo& iterPon, int* card, int cudaCap, int noSM) {
//	float_type* d_Pparcelas = NULL;
//	float_type aux = 0;
//
//	int /*threads         = noThreadsPQ(card[IDX1F(k)]), // potência de dois imediatamente superior ao número de barras ligadas à k
//				   threadsPerBlock = noThreadsPerBlock(threads, cudaCap, noSM),*/
//		noBlocks = noMaxBlocks(cudaCap, noSM);
//	checkCudaErrors(cudaMalloc(&d_Pparcelas, sizeof(float_type) * noThreadsPQ(h_sistema.nB)));
//	checkCudaErrors(cudaMemset(d_Pparcelas, 0, sizeof(float_type) * noThreadsPQ(h_sistema.nB)));
//
//	//	cudaDeviceSynchronize();
//		//printf("Threads = %d\n", noThreadsPQ(h_sistema.nB));
//
//		//printf("d_calcPf: tamanho = %d\n", tamanho);
//		//calcP<<<noBlocks, threadsPerBlock>>>(k, d_Pparcelas, d_sistema, barra, ramo, threads);
//	calcP << <1, noThreadsPQ(h_sistema.nB) >> > (k, d_Pparcelas, d_sistema, d_barra, d_ramo, noThreadsPQ(h_sistema.nB));
//	//printf("d_calcPf: calcPsci executada!\n");
//
////	cudaDeviceSynchronize();
//
//	//printf("oioioioioioioioioi\n\n");
//
//	checkCudaErrors(cudaMemcpy(iterPon.Pcalc + (k - 1), &(d_Pparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); // guarda valor de Pk no vetor d_P
//	if (d_Pparcelas) { checkCudaErrors(cudaFree(d_Pparcelas)); }
//
//	if (global::verbose_mode) {
//		// chequemos o resultado
//		checkCudaErrors(cudaMemcpy(&aux, iterPon.Pcalc + (k - 1), sizeof(float_type), cudaMemcpyDeviceToHost));
//		printf("P_%d\t= %f\n", k, aux);
//	}
//}

int min_(int x, int y) {
	return y ^ ((x ^ y) & -(x < y));
}

void d_dnCalcPf_sha(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) { 

	const int threads = noThreadsPQ(sistPon.nB);
	const int threadsPerBlock = min_(512, threads);
	const int blocks = ceil((float_type)threads / (float_type)threadsPerBlock);

	float_type *d_Out = nullptr, *h_Out = nullptr;
	checkCudaErrors(cudaMalloc(&d_Out, sizeof(float_type) * blocks));
	h_Out = (float_type*)malloc(sizeof(float_type) * blocks);
	//checkCudaErrors(cudaMemset(d_Pparcelas, 0, sizeof(float_type) * noThreadsPQ(sistPon.nB)));
	
	//checkCudaErrors(cudaDeviceSynchronize());
	calcP_sha<<<blocks, threadsPerBlock, sizeof(float_type)* threads>>> (k, d_sistema, d_barra, d_ramo, d_Out);
	//checkCudaErrors(cudaDeviceSynchronize());

	checkCudaErrors(cudaMemcpy(h_Out, d_Out, sizeof(float_type) * blocks, cudaMemcpyDeviceToHost));

	float_type acc = 0;
	for (int i = 0; i < blocks; i++) {
		acc += h_Out[i];
	}

	checkCudaErrors(cudaMemcpy(iterPon.Pcalc + (k - 1), &(acc), sizeof(float_type), cudaMemcpyHostToDevice)); // guarda valor de Pk no vetor d_P

	if (d_Out == nullptr) { checkCudaErrors(cudaFree(d_Out)); }
	if (h_Out == nullptr) { free(d_Out); }

	if (global::verbose_mode) {
		float_type aux = 0;
		// chequemos o resultado
		checkCudaErrors(cudaMemcpy(&aux, iterPon.Pcalc + (k - 1), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("P_%d\t= %f\n", k, aux);
	}
}

void d_dnCalcPf(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) { //(const int k, sistema& h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, iterativo& iterPon) {
	float_type* d_Pparcelas = NULL;
	float_type aux = 0;

	checkCudaErrors(cudaMalloc(&d_Pparcelas, sizeof(float_type) * noThreadsPQ(sistPon.nB)));
	checkCudaErrors(cudaMemset(d_Pparcelas, 0, sizeof(float_type) * noThreadsPQ(sistPon.nB)));


	calcP <<<1, noThreadsPQ(sistPon.nB)>>> (k, d_Pparcelas, d_sistema, d_barra, d_ramo, noThreadsPQ(sistPon.nB));

	checkCudaErrors(cudaMemcpy(iterPon.Pcalc + (k - 1), &(d_Pparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); // guarda valor de Pk no vetor d_P
	if (d_Pparcelas) { checkCudaErrors(cudaFree(d_Pparcelas)); }

	if (global::verbose_mode) {
		// chequemos o resultado
		checkCudaErrors(cudaMemcpy(&aux, iterPon.Pcalc + (k - 1), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("P_%d\t= %f\n", k, aux);
	}
}

//***************************************************************************************************************************************************************************************
//***************************************************************************************************************************************************************************************
//***************************************************************************************************************************************************************************************
//***************************************************************************************************************************************************************************************

//idx = 0:nL-1
__global__ void calcQ(const int k, float_type* d_Qparcelas, sistema* d_sistema, barra* d_barra, ramo* d_ramo, int threads){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
	if(idx <= d_sistema->nB){
//		printf("(%d)\n", idx);
		float_type aux = 0;

		//aux = phif(k, idx, d_sistema, d_ramo); // defasagem do transformador
		//printf("(%d); aux = %f\n", idx, aux);
		aux = d_barra->theta[IDX1F(k)] - d_barra->theta[IDX1F(idx)]; // theta_k para[i]
		//printf("(%d); aux = %f\n", idx, aux);
		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux) - _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux);
		//printf("(%d); aux = %f\n", idx, aux);
		aux *= d_barra->V[IDX1F(idx)];
//		printf("(%d); aux = %f\n", idx, aux);

		d_Qparcelas[idx-1] = d_barra->V[IDX1F(k)]*aux;
	}
//	printf("hello from (%d)\n", idx);

	idx--;

//		__syncthreads();
//
//		if(idx == 1){
//			float_type aux = 0;
//			printf("somando em %d\n", idx);
//			for (int i = 0; i < d_sistema->nB; i++){
//				printf("d_Qparcelas[%d] = %e\n", i, d_Qparcelas[i]);
//				aux += d_Qparcelas[i];
//			}
//
//			printf("ans = %e\n", aux);
//		}

	__syncthreads();

	// da pra ganhar eficiencia aqui ainda... ver memória shared e aqueles slides
	// o que acontece quando tiver mais de um bloco no mesmo P??????????
	for (int s = threads/2; s > 0; s >>= 1){
		if (idx < s){
			d_Qparcelas[idx] += d_Qparcelas[idx + s];
			//printf("idx = %d;\ts = %d\n", idx, s);
		}
		__syncthreads();
	}
//		if(idx == 0){
//			printf("ansRed[%d] = %e\n", k, d_Qparcelas[0]);
//		}

	// ATENÇÃO: MULTIPLICAR POR VK na cpu ou em outro kernel otimiza alguma coisa?????????????????????????????????????????????????????????????????????????????
}

__global__ void calcQ_sha(const int k, sistema* d_sistema, barra* barra, ramo* ramo, float_type* d_Out) {
	extern __shared__ float_type s_data[];
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	int tid = threadIdx.x;
	if (idx <= d_sistema->nB) {
		//		printf("(%d)\n", idx);
		float_type aux = 0;

		//aux = phif(k, idx, d_sistema, ramo); // defasagem do transformador
//		printf("(%d); aux = %f\n", idx, aux);
		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; // theta_k para[i]
//		printf("(%d); aux = %f\n", idx, aux);
		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux) - _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux);
		//		printf("(%d); aux = %f\n", idx, aux);
		aux *= barra->V[IDX1F(idx)];
		//		printf("(%d); aux = %f\n", idx, aux);

		s_data[tid] = barra->V[IDX1F(k)] * aux;
	}
	else {
		s_data[tid] = 0.;
	}
	//printf("hello from (%d)\n", idx);

	//idx--;

	__syncthreads();

	//if (tid == 0) {
	//	for (int i = 0; i < 16; i++) {
	//		printf("s_data[%d] = %f\n", i, s_data[i]);
	//	}
	//}

	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			s_data[tid] += s_data[tid + s];
			//printf("idx = %d;\ts = %d\n", tid, s);
		}
		__syncthreads();
	}
	if (tid == 0) {
		d_Out[blockIdx.x] = s_data[0];
	}
}

//void d_calcQf(const int k, sistema &h_sistema, sistema* d_sistema, barra* d_barra, ramo* d_ramo, iterativo* d_iterativo, iterativo& iterPon, int* card, int cudaCap, int noSM){
//    float_type* d_Qparcelas = NULL;
//    float_type aux = 0;
//    int /*threads         = noThreadsPQ(card[IDX1F(k)]),
//				   threadsPerBlock = noThreadsPerBlock(threads, cudaCap, noSM),*/
//				   noBlocks        = noMaxBlocks(cudaCap, noSM);
//
//    checkCudaErrors(cudaMalloc(&d_Qparcelas, sizeof(float_type)*noThreadsPQ(h_sistema.nB)));
//	checkCudaErrors(cudaMemset(d_Qparcelas, 0, sizeof(float_type)*noThreadsPQ(h_sistema.nB)));
//
//    //printf("[tamanho = %d\n", tamanho);
//	//calcQ<<<noBlocks, threadsPerBlock>>>(k, d_Qparcelas, d_sistema, d_barra, d_ramo, threads);
//	calcQ<<<1, noThreadsPQ(h_sistema.nB)>>>(k, d_Qparcelas, d_sistema, d_barra, d_ramo, noThreadsPQ(h_sistema.nB));
//
////	cudaDeviceSynchronize();
//
//	//printf("oioioioioioioioioi\n\n");
//
//	checkCudaErrors(cudaMemcpy(iterPon.Qcalc+(k-1), &(d_Qparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); // guarda valor de Qk no vetor d_Q
//	if (d_Qparcelas) { checkCudaErrors(cudaFree(d_Qparcelas)); }
//
//	// chequemos o resultado
//	if(global::verbose_mode){
//		checkCudaErrors(cudaMemcpy(&aux, (iterPon.Qcalc+(k-1)), sizeof(float_type), cudaMemcpyDeviceToHost));
//		printf("Q_%d\t= %f\n", k, aux);
//	}
//}

void d_dnCalcQf_sha(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) {

	const int threads = noThreadsPQ(sistPon.nB);
	const int threadsPerBlock = min_(512, threads);
	const int blocks = ceil((float_type)threads / (float_type)threadsPerBlock);

	float_type* d_Out = nullptr, * h_Out = nullptr;
	checkCudaErrors(cudaMalloc(&d_Out, sizeof(float_type) * blocks));
	h_Out = (float_type*)malloc(sizeof(float_type) * blocks);
	//checkCudaErrors(cudaMemset(d_Pparcelas, 0, sizeof(float_type) * noThreadsPQ(sistPon.nB)));

	//checkCudaErrors(cudaDeviceSynchronize());
	calcQ_sha <<<blocks, threadsPerBlock, sizeof(float_type)* threads >>> (k, d_sistema, d_barra, d_ramo, d_Out);
	//checkCudaErrors(cudaDeviceSynchronize());

	checkCudaErrors(cudaMemcpy(h_Out, d_Out, sizeof(float_type) * blocks, cudaMemcpyDeviceToHost));

	float_type acc = 0;
	for (int i = 0; i < blocks; i++) {
		acc += h_Out[i];
	}

	checkCudaErrors(cudaMemcpy(iterPon.Qcalc + (k - 1), &(acc), sizeof(float_type), cudaMemcpyHostToDevice)); // guarda valor de Pk no vetor d_P

	if (d_Out == nullptr) { checkCudaErrors(cudaFree(d_Out)); }
	if (h_Out == nullptr) { free(d_Out); }

	if (global::verbose_mode) {
		float_type aux = 0;
		// chequemos o resultado
		checkCudaErrors(cudaMemcpy(&aux, (iterPon.Qcalc + (k - 1)), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("Q_%d\t= %f\n", k, aux);
	}
}

void d_dnCalcQf(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) {
	float_type* d_Qparcelas = NULL;
	float_type aux = 0;

	checkCudaErrors(cudaMalloc(&d_Qparcelas, sizeof(float_type) * noThreadsPQ(sistPon.nB)));
	checkCudaErrors(cudaMemset(d_Qparcelas, 0, sizeof(float_type) * noThreadsPQ(sistPon.nB)));

	calcQ <<<1, noThreadsPQ(sistPon.nB) >>> (k, d_Qparcelas, d_sistema, d_barra, d_ramo, noThreadsPQ(sistPon.nB));

	checkCudaErrors(cudaMemcpy(iterPon.Qcalc + (k - 1), &(d_Qparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); // guarda valor de Qk no vetor d_Q
	if (d_Qparcelas) { checkCudaErrors(cudaFree(d_Qparcelas)); }

	// chequemos o resultado
	if (global::verbose_mode) {
		checkCudaErrors(cudaMemcpy(&aux, (iterPon.Qcalc + (k - 1)), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("Q_%d\t= %f\n", k, aux);
	}
}

void d_dnCalculePQ(sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, barra& barraPon, iterativo& iterPon) {
	// Pcalc e Qcalc s�o percorridos por iteradores an�logos a barrasPQ e barrasPV
	// <N�O!!!>Pcalc e Qcalc possuem nB entradas e a entrada k � relativa � k-�sima barra

	// Calcula os valores de Pk...
	for (int i = 1; i <= sistPon.nB; i++) {
		d_dnCalcPf_sha(i, d_sistema, d_barra, d_ramo, sistPon, iterPon);
	}
	// Calcula os valores de Qk...
	for (int i = 1; i <= sistPon.nB; i++) {
		d_dnCalcQf_sha(i, d_sistema, d_barra, d_ramo, sistPon, iterPon);
	}
}

void d_calculePQ(sistema* d_sistema, ramo* d_ramo, sistema &sistPon, barra &barraPon, iterativo &iterPon){
	checkCudaErrors(cudaGetLastError());

	int threadsPerBlock = 128;

	calcPeficiente<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);
	calcQeficiente<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);

	checkCudaErrors(cudaGetLastError());
}

void d_calculePQ_0based(sistema* d_sistema, ramo* d_ramo, sistema &sistPon, barra &barraPon, iterativo &iterPon, cudaStream_t* streams){
	checkCudaErrors(cudaGetLastError());

	int threadsPerBlock = 128;

	calcPeficiente_0based<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock, 0, streams[0]>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);
	//calcPeficiente_0based_sha<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock, (sistPon.nnzY + sistPon.nB + 1) * sizeof(int), streams[0]>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);
	calcQeficiente_0based<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock, 0, streams[1]>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);

	checkCudaErrors(cudaGetLastError());
}

// 1 thread para cada ramo
void __global__ d_P_Eficiente(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Pdp[IDX1F(idx)] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
		ramoPon.Ppd[IDX1F(idx)] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x*cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y*sin(-aux));
	}
}

// 1 thread para cada ramo
void __global__ d_Q_Eficiente(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Qdp[IDX1F(idx)] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*sin(aux));
		ramoPon.Qpd[IDX1F(idx)] = (-barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * (sistPon.Y[IDX2F(m, k, sistPon.nB)].y + bshf(m, k, sistPon, ramoPon)) + barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y*cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x*sin(-aux));	
	}
}

void calcFluxf_Eficiente(sistema &h_sistema, barra &h_barra, ramo &h_ramo, iterativo &h_iterativo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, const iterativo &iterPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nL;
	dim3 dimBlock(16*deviceprop.warpSize, 1);
	dim3 dimGrid(
			(int) ceil(
					((float) tamanho) / (float) (8*deviceprop.warpSize)), 1);

	// Calcula os valores de Pkm e Pmk
	d_P_Eficiente<<<dimGrid, dimBlock>>>(sistPon, barraPon, ramoPon);
	
	// Calcula os valores de Qkm e Qmk
	d_Q_Eficiente<<<dimGrid, dimBlock>>>(sistPon, barraPon, ramoPon);
}

// 1 thread para cada ramo
void __global__ d_P_Eficiente_Sp(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Pdp[IDX1F(idx)] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux)  - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		ramoPon.Ppd[IDX1F(idx)] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(-aux));
	}
}

// 1 thread para cada ramo
void __global__ d_Q_Eficiente_Sp(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Qdp[IDX1F(idx)] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y*cos(aux)  - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x*sin(aux));
		ramoPon.Qpd[IDX1F(idx)] = (-barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * (sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y + bshf(m, k, sistPon, ramoPon)) + barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y*cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x*sin(-aux));
	}
}

void calcFluxf_Eficiente_Sp(sistema &h_sistema, barra &h_barra, ramo &h_ramo, iterativo &h_iterativo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, const iterativo &iterPon, cudaDeviceProp deviceprop, cudaStream_t* streams) {
	int tamanho = sistPon.nL;
	dim3 dimBlock(16*deviceprop.warpSize, 1);
	dim3 dimGrid(
			(int) ceil(
					((float) tamanho) / (float) (8*deviceprop.warpSize)), 1);

	//printf(" tamanho = %d, threads/block = %d, blocks/grid = %d\n\n", tamanho, 16 * deviceprop.warpSize, (int)ceil(((float)tamanho) / (float)(16 * deviceprop.warpSize)));

	// Calcula os valores de Pkm e Pmk
	d_P_Eficiente_Sp<<<dimGrid, dimBlock, 0, streams[0]>>>(sistPon, barraPon, ramoPon);
	
	// Calcula os valores de Qkm e Qmk
	d_Q_Eficiente_Sp<<<dimGrid, dimBlock, 0, streams[1]>>>(sistPon, barraPon, ramoPon);
}

void __global__ d_P_dp(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
	ramoPon.Pdp[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
}

void __global__ d_Q_dp(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
	ramoPon.Qdp[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*sin(aux));
}

void __global__ d_P_pd(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
	ramoPon.Ppd[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
}

void __global__ d_Q_pd(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
	ramoPon.Qpd[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*sin(aux));
}

//muito ineficiente
void calcFluxf(sistema &h_sistema, barra &h_barra, ramo &h_ramo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nL; // nPQ
	dim3 dimBlock(3*deviceprop.warpSize, 1);
	dim3 dimGrid(
			(int) ceil(
					((float) tamanho) / (float) (3*deviceprop.warpSize)), 1);

	// Calcula os valores de Pkm...
	for (int i = 1; i <= sistPon.nL; i++) {
		d_P_dp<<<dimGrid, dimBlock>>>(h_ramo.de[IDX1F(i)], h_ramo.para[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}
	// Calcula os valores de Qkm...
	for (int i = 1; i <= sistPon.nL; i++) {
		d_Q_dp<<<dimGrid, dimBlock>>>(h_ramo.de[IDX1F(i)], h_ramo.para[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}
	// Calcula os valores de Pmk
	for (int i = 1; i <= sistPon.nL; i++) {
		d_P_pd<<<dimGrid, dimBlock>>>(h_ramo.para[IDX1F(i)], h_ramo.de[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}
	// Calcula os valores de Qmk...
	for (int i = 1; i <= sistPon.nL; i++) {
		d_Q_pd<<<dimGrid, dimBlock>>>(h_ramo.para[IDX1F(i)], h_ramo.de[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}
//	cudaDeviceSynchronize();

	checkCudaErrors(cudaMemcpy(h_ramo.Ppd,   ramoPon.Ppd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Pdp,   ramoPon.Pdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qpd,   ramoPon.Qpd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qdp,   ramoPon.Qdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

}

// *eficientes*********************************************
void __global__ d_P_ef(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int i = (threadIdx.x + blockDim.x * blockIdx.x);
	int idx = i + 1;

	if (i < sistPon.nL) {
		int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];

		float_type aux = 0;
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Pdp[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x * cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y * sin(aux));
		ramoPon.Ppd[i] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y * sin(-aux));
	}
	//printf("%d -> %d\n", threadIdx.x + blockDim.x * blockIdx.x + 1, idx);
}

void __global__ d_Q_ef(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int i = (threadIdx.x + blockDim.x * blockIdx.x);
	int idx = i + 1;

	if (i < sistPon.nL) {
		int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];

		float_type aux = 0;
		aux = phif(k, m, sistPon, barraPon, ramoPon); // defasagem do transformador
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; // theta_k para[i]
		ramoPon.Qdp[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y * cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x * sin(aux));
		ramoPon.Qpd[i] = (-barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * (sistPon.Y[IDX2F(m, k, sistPon.nB)].y + bshf(m, k, sistPon, ramoPon)) + barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x * sin(-aux));
	}
	//printf("%d -> %d\n", threadIdx.x + blockDim.x * blockIdx.x + 1, idx);
}

// ********************************************************

void calcFluxf_ef(sistema &h_sistema, barra &h_barra, ramo &h_ramo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nL;
	
	int threadsPerBlock = 128;
	int blocksPerGrid = (int)ceil((float)(tamanho) / (float)threadsPerBlock);

	//dim3 dimBlock(3*deviceprop.warpSize, 1);
	//dim3 dimGrid(
	//		(int) ceil(
	//				((float) tamanho) / (float) (3*deviceprop.warpSize)), 1);

	// Calcula os valores de P...
	d_P_ef<<< blocksPerGrid, threadsPerBlock >>>(sistPon, barraPon, ramoPon);
	d_Q_ef<<< blocksPerGrid, threadsPerBlock >>>(sistPon, barraPon, ramoPon);

	
	checkCudaErrors(cudaMemcpy(h_ramo.Ppd,   ramoPon.Ppd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Pdp,   ramoPon.Pdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qpd,   ramoPon.Qpd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qdp,   ramoPon.Qdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

}