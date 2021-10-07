
int noThreadsPQ(int nl){
	int aux = 1;
	while (nl > aux){
		aux <<= 1;
		
	}
	return aux;
}

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

__host__ __device__ int coeffPos(int lin, int col, int* csrRowPtrY, int* csrColIndY, int nnzY) {
	for (int i = csrRowPtrY[lin]; i < csrRowPtrY[lin + 1]; i++) {
		if (csrColIndY[i] == col) {
			return i;
		}
	}
	return -1;
}

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

__host__ __device__ float_type bshf(int a, int b, const sistema &sistPon, const ramo &ramoPon) {
	for (int i = 0; i < sistPon.nL; i++) { 
		if (ramoPon.de[i] == a) { 
			if (ramoPon.para[i] == b) { 
				return (ramoPon.bsh[i]);
			}
		}
		else {
			if (ramoPon.para[i] == a) { 
				if (ramoPon.de[i] == b) { 
					return (ramoPon.bsh[i]);
				}
			}
		}
	}
	return(0.);
}

__global__ void calcPeficiente(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
 	if(idx <= sistPon.nB){
 	
 		float_type aux = 0;
		float_type acc = 0;
		
		for(int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++){
			
			aux = barraPon.theta[IDX1F(sistPon.cooRowIndY[IDX1F(i)])] - barraPon.theta[IDX1F(sistPon.csrColIndY[IDX1F(i)])]; 

			aux = _cuReal(sistPon.spYval[IDX1F(i)]) * cos(aux) + _cuImag(sistPon.spYval[IDX1F(i)]) * sin(aux);
			aux *= barraPon.V[IDX1F(sistPon.csrColIndY[IDX1F(i)])];
			acc += aux;
			aux = 0;
		}
 		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; 
		
 	}
}

__global__ void calcQeficiente(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
 	if(idx <= sistPon.nB){
 	
 		float_type aux = 0;
		float_type acc = 0;
		
		for(int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++){
			
			aux = barraPon.theta[IDX1F(sistPon.cooRowIndY[IDX1F(i)])] - barraPon.theta[IDX1F(sistPon.csrColIndY[IDX1F(i)])]; 

			aux = _cuReal(sistPon.spYval[IDX1F(i)]) * sin(aux) - _cuImag(sistPon.spYval[IDX1F(i)]) * cos(aux);
			aux *= barraPon.V[IDX1F(sistPon.csrColIndY[IDX1F(i)])];
			acc += aux;
			aux = 0;
		}
 		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Qcalc[IDX1F(idx)] = acc; 
		
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
		
		float_type aux = 0;
		float_type acc = 0;
		
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; 

			aux = _cuReal(sistPon.spYval[i]) * cos(aux) + _cuImag(sistPon.spYval[i]) * sin(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; 
		
	}
}

__global__ void calcPeficiente_0based(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nB) {
		
		float_type aux = 0;
		float_type acc = 0;
		
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; 

			aux = _cuReal(sistPon.spYval[i]) * cos(aux) + _cuImag(sistPon.spYval[i]) * sin(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Pcalc[IDX1F(idx)] = acc; 
		
	}
}

__global__ void calcQeficiente_0based(sistema* d_sistema, ramo* d_ramo, sistema sistPon, barra barraPon, iterativo iterPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nB) {
		
		float_type aux = 0;
		float_type acc = 0;
		
		for (int i = sistPon.csrRowPtrY[IDX1F(idx)]; i < sistPon.csrRowPtrY[IDX1F(idx) + 1]; i++) {
			
			aux = barraPon.theta[sistPon.cooRowIndY[i]] - barraPon.theta[sistPon.csrColIndY[i]]; 

			aux = _cuReal(sistPon.spYval[i]) * sin(aux) - _cuImag(sistPon.spYval[i]) * cos(aux);
			aux *= barraPon.V[sistPon.csrColIndY[i]];
			acc += aux;
			aux = 0;
		}
		acc *= barraPon.V[IDX1F(idx)];
		iterPon.Qcalc[IDX1F(idx)] = acc; 
		
	}
}

__global__ void calcP_sha(const int k, sistema* d_sistema, barra* barra, ramo* ramo, float_type* d_Out) {
	extern __shared__ float_type s_data[];
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	int tid = threadIdx.x;
	if (idx <= d_sistema->nB) {
		
		float_type aux = 0;

		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; 

		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux) + _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux);
		
		aux *= barra->V[IDX1F(idx)];
		
		s_data[tid] = barra->V[IDX1F(k)] * aux;
	}
	else {
		s_data[tid] = 0.;
	}

	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			s_data[tid] += s_data[tid + s];
			
		}
		__syncthreads();
	}
	if (tid == 0) {
		d_Out[blockIdx.x] = s_data[0];
	}
}

__global__ void calcP(const int k, float_type* d_Pparcelas, sistema* d_sistema, barra* barra, ramo* ramo, int threads){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
	if(idx <= d_sistema->nB){
		float_type aux = 0;

		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; 

		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux) + _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux);

		aux *= barra->V[IDX1F(idx)];

		d_Pparcelas[idx-1] = barra->V[IDX1F(k)]*aux;
	}

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
			
		}
		__syncthreads();
	}
}

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

	calcP_sha<<<blocks, threadsPerBlock, sizeof(float_type)* threads>>> (k, d_sistema, d_barra, d_ramo, d_Out);

	checkCudaErrors(cudaMemcpy(h_Out, d_Out, sizeof(float_type) * blocks, cudaMemcpyDeviceToHost));

	float_type acc = 0;
	for (int i = 0; i < blocks; i++) {
		acc += h_Out[i];
	}

	checkCudaErrors(cudaMemcpy(iterPon.Pcalc + (k - 1), &(acc), sizeof(float_type), cudaMemcpyHostToDevice)); 

	if (d_Out == nullptr) { checkCudaErrors(cudaFree(d_Out)); }
	if (h_Out == nullptr) { free(d_Out); }

	if (global::verbose_mode) {
		float_type aux = 0;
		
		checkCudaErrors(cudaMemcpy(&aux, iterPon.Pcalc + (k - 1), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("P_%d\t= %f\n", k, aux);
	}
}

void d_dnCalcPf(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) { 
	float_type* d_Pparcelas = NULL;
	float_type aux = 0;

	checkCudaErrors(cudaMalloc(&d_Pparcelas, sizeof(float_type) * noThreadsPQ(sistPon.nB)));
	checkCudaErrors(cudaMemset(d_Pparcelas, 0, sizeof(float_type) * noThreadsPQ(sistPon.nB)));

	calcP <<<1, noThreadsPQ(sistPon.nB)>>> (k, d_Pparcelas, d_sistema, d_barra, d_ramo, noThreadsPQ(sistPon.nB));

	checkCudaErrors(cudaMemcpy(iterPon.Pcalc + (k - 1), &(d_Pparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); 
	if (d_Pparcelas) { checkCudaErrors(cudaFree(d_Pparcelas)); }

	if (global::verbose_mode) {
		
		checkCudaErrors(cudaMemcpy(&aux, iterPon.Pcalc + (k - 1), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("P_%d\t= %f\n", k, aux);
	}
}

__global__ void calcQ(const int k, float_type* d_Qparcelas, sistema* d_sistema, barra* d_barra, ramo* d_ramo, int threads){
	int idx = threadIdx.x + blockDim.x * blockIdx.x +1;
	if(idx <= d_sistema->nB){
		float_type aux = 0;

		aux = d_barra->theta[IDX1F(k)] - d_barra->theta[IDX1F(idx)]; 
		
		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux) - _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux);
		
		aux *= d_barra->V[IDX1F(idx)];

		d_Qparcelas[idx-1] = d_barra->V[IDX1F(k)]*aux;
	}

	idx--;

	__syncthreads();

	for (int s = threads/2; s > 0; s >>= 1){
		if (idx < s){
			d_Qparcelas[idx] += d_Qparcelas[idx + s];
			
		}
		__syncthreads();
	}

}

__global__ void calcQ_sha(const int k, sistema* d_sistema, barra* barra, ramo* ramo, float_type* d_Out) {
	extern __shared__ float_type s_data[];
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	int tid = threadIdx.x;
	if (idx <= d_sistema->nB) {
		
		float_type aux = 0;

		aux = barra->theta[IDX1F(k)] - barra->theta[IDX1F(idx)]; 

		aux = _cuReal(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * sin(aux) - _cuImag(d_sistema->Y[IDX2F(k, idx, d_sistema->nB)]) * cos(aux);
		
		aux *= barra->V[IDX1F(idx)];
		
		s_data[tid] = barra->V[IDX1F(k)] * aux;
	}
	else {
		s_data[tid] = 0.;
	}

	__syncthreads();

	for (int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			s_data[tid] += s_data[tid + s];
			
		}
		__syncthreads();
	}
	if (tid == 0) {
		d_Out[blockIdx.x] = s_data[0];
	}
}

void d_dnCalcQf_sha(const int k, sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, iterativo& iterPon) {
	const int threads = noThreadsPQ(sistPon.nB);
	const int threadsPerBlock = min_(512, threads);
	const int blocks = ceil((float_type)threads / (float_type)threadsPerBlock);

	float_type* d_Out = nullptr, * h_Out = nullptr;
	checkCudaErrors(cudaMalloc(&d_Out, sizeof(float_type) * blocks));
	h_Out = (float_type*)malloc(sizeof(float_type) * blocks);

	calcQ_sha <<<blocks, threadsPerBlock, sizeof(float_type)* threads >>> (k, d_sistema, d_barra, d_ramo, d_Out);

	checkCudaErrors(cudaMemcpy(h_Out, d_Out, sizeof(float_type) * blocks, cudaMemcpyDeviceToHost));

	float_type acc = 0;
	for (int i = 0; i < blocks; i++) {
		acc += h_Out[i];
	}

	checkCudaErrors(cudaMemcpy(iterPon.Qcalc + (k - 1), &(acc), sizeof(float_type), cudaMemcpyHostToDevice)); 

	if (d_Out == nullptr) { checkCudaErrors(cudaFree(d_Out)); }
	if (h_Out == nullptr) { free(d_Out); }

	if (global::verbose_mode) {
		float_type aux = 0;
		
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

	checkCudaErrors(cudaMemcpy(iterPon.Qcalc + (k - 1), &(d_Qparcelas[0]), sizeof(float_type), cudaMemcpyDeviceToDevice)); 
	if (d_Qparcelas) { checkCudaErrors(cudaFree(d_Qparcelas)); }

	if (global::verbose_mode) {
		checkCudaErrors(cudaMemcpy(&aux, (iterPon.Qcalc + (k - 1)), sizeof(float_type), cudaMemcpyDeviceToHost));
		printf("Q_%d\t= %f\n", k, aux);
	}
}

void d_dnCalculePQ(sistema* d_sistema, barra* d_barra, ramo* d_ramo, sistema& sistPon, barra& barraPon, iterativo& iterPon) {
	for (int i = 1; i <= sistPon.nB; i++) {
		d_dnCalcPf_sha(i, d_sistema, d_barra, d_ramo, sistPon, iterPon);
	}

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

	calcQeficiente_0based<<<(int) ceil(((float) sistPon.nB) / (float) (threadsPerBlock)), threadsPerBlock, 0, streams[1]>>>(d_sistema, d_ramo, sistPon, barraPon, iterPon);

	checkCudaErrors(cudaGetLastError());
}

void __global__ d_P_Eficiente(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
		ramoPon.Pdp[IDX1F(idx)] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
		ramoPon.Ppd[IDX1F(idx)] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x*cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y*sin(-aux));
	}
}

void __global__ d_Q_Eficiente(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
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

	d_P_Eficiente<<<dimGrid, dimBlock>>>(sistPon, barraPon, ramoPon);

	d_Q_Eficiente<<<dimGrid, dimBlock>>>(sistPon, barraPon, ramoPon);
}

void __global__ d_P_Eficiente_Sp(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
		ramoPon.Pdp[IDX1F(idx)] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(aux)  - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(k), IDX1F(m), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(aux));
		ramoPon.Ppd[IDX1F(idx)] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].x * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.spYval[coeffPos(IDX1F(m), IDX1F(k), sistPon.csrRowPtrY, sistPon.csrColIndY, sistPon.nnzY)].y * sin(-aux));
	}
}

void __global__ d_Q_Eficiente_Sp(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int idx = threadIdx.x + blockDim.x * blockIdx.x + 1;
	if (idx <= sistPon.nL){
		float_type aux = 0;
		const int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
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

	d_P_Eficiente_Sp<<<dimGrid, dimBlock, 0, streams[0]>>>(sistPon, barraPon, ramoPon);

	d_Q_Eficiente_Sp<<<dimGrid, dimBlock, 0, streams[1]>>>(sistPon, barraPon, ramoPon);
}

void __global__ d_P_dp(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); 
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
	ramoPon.Pdp[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
}

void __global__ d_Q_dp(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); 
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
	ramoPon.Qdp[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*sin(aux));
}

void __global__ d_P_pd(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); 
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
	ramoPon.Ppd[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*sin(aux));
}

void __global__ d_Q_pd(int k, int m, const sistema sistPon, const barra barraPon, const ramo ramoPon, const int i) {
	float_type aux = 0;
	aux = phif(k, m, sistPon, barraPon, ramoPon); 
	aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
	ramoPon.Qpd[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y*cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x*sin(aux));
}

void calcFluxf(sistema &h_sistema, barra &h_barra, ramo &h_ramo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nL; 
	dim3 dimBlock(3*deviceprop.warpSize, 1);
	dim3 dimGrid(
			(int) ceil(
					((float) tamanho) / (float) (3*deviceprop.warpSize)), 1);

	for (int i = 1; i <= sistPon.nL; i++) {
		d_P_dp<<<dimGrid, dimBlock>>>(h_ramo.de[IDX1F(i)], h_ramo.para[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}

	for (int i = 1; i <= sistPon.nL; i++) {
		d_Q_dp<<<dimGrid, dimBlock>>>(h_ramo.de[IDX1F(i)], h_ramo.para[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}

	for (int i = 1; i <= sistPon.nL; i++) {
		d_P_pd<<<dimGrid, dimBlock>>>(h_ramo.para[IDX1F(i)], h_ramo.de[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}

	for (int i = 1; i <= sistPon.nL; i++) {
		d_Q_pd<<<dimGrid, dimBlock>>>(h_ramo.para[IDX1F(i)], h_ramo.de[IDX1F(i)], sistPon, barraPon, ramoPon, IDX1F(i));
	}

	checkCudaErrors(cudaMemcpy(h_ramo.Ppd,   ramoPon.Ppd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Pdp,   ramoPon.Pdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qpd,   ramoPon.Qpd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qdp,   ramoPon.Qdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

}

void __global__ d_P_ef(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int i = (threadIdx.x + blockDim.x * blockIdx.x);
	int idx = i + 1;

	if (i < sistPon.nL) {
		int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];

		float_type aux = 0;
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
		ramoPon.Pdp[i] = (barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x * cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y * sin(aux));
		ramoPon.Ppd[i] = (barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y * sin(-aux));
	}

}

void __global__ d_Q_ef(const sistema sistPon, const barra barraPon, const ramo ramoPon) {
	int i = (threadIdx.x + blockDim.x * blockIdx.x);
	int idx = i + 1;

	if (i < sistPon.nL) {
		int k = ramoPon.de[IDX1F(idx)], m = ramoPon.para[IDX1F(idx)];

		float_type aux = 0;
		aux = phif(k, m, sistPon, barraPon, ramoPon); 
		aux += barraPon.theta[IDX1F(k)] - barraPon.theta[IDX1F(m)]; 
		ramoPon.Qdp[i] = (-barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(k)] * (sistPon.Y[IDX2F(k, m, sistPon.nB)].y + bshf(k, m, sistPon, ramoPon)) + barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].y * cos(aux) - barraPon.V[IDX1F(k)] * barraPon.V[IDX1F(m)] * sistPon.Y[IDX2F(k, m, sistPon.nB)].x * sin(aux));
		ramoPon.Qpd[i] = (-barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(m)] * (sistPon.Y[IDX2F(m, k, sistPon.nB)].y + bshf(m, k, sistPon, ramoPon)) + barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].y * cos(-aux) - barraPon.V[IDX1F(m)] * barraPon.V[IDX1F(k)] * sistPon.Y[IDX2F(m, k, sistPon.nB)].x * sin(-aux));
	}

}

void calcFluxf_ef(sistema &h_sistema, barra &h_barra, ramo &h_ramo, const sistema &sistPon, const barra &barraPon, const ramo &ramoPon, cudaDeviceProp deviceprop) {
	int tamanho = sistPon.nL;

	int threadsPerBlock = 128;
	int blocksPerGrid = (int)ceil((float)(tamanho) / (float)threadsPerBlock);

	d_P_ef<<< blocksPerGrid, threadsPerBlock >>>(sistPon, barraPon, ramoPon);
	d_Q_ef<<< blocksPerGrid, threadsPerBlock >>>(sistPon, barraPon, ramoPon);

	checkCudaErrors(cudaMemcpy(h_ramo.Ppd,   ramoPon.Ppd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Pdp,   ramoPon.Pdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qpd,   ramoPon.Qpd,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(h_ramo.Qdp,   ramoPon.Qdp,   sistPon.nL * sizeof(float_type),          cudaMemcpyDeviceToHost));

}