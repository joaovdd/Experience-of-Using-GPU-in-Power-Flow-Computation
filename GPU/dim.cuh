__global__ void d_showVec(float_type* show, int dim) {
	printf("ShowVec:\n");
	printf("[%f\n", show[0]);

	for (int i = 1; i < dim - 1; i++){
		printf(" %f\n", show[i]);
	}
	printf(" %f]\n", show[dim-1]);
}

__global__ void d_showVecI(int* show, int dim){
	printf("ShowVec:\n");
	printf("[%d\n", show[0]);

	for (int i = 1; i < dim - 1; i++){
		printf(" %d\n", show[i]);
	}
	printf(" %d]\n", show[dim-1]);
}

__global__ void d_showVecC(complex_type* show, int dim){
	printf("ShowVec:\n");
	printf("[%f + j * %f\n", show[0].x, show[0].y);

	for (int i = 1; i < dim - 1; i++){
		printf(" %f + j * %f\n", show[i].x, show[i].y);
	}
	printf(" %f + j * %f]\n", show[dim-1].x, show[dim-1].y);
}

void d_showVecf(float_type* show, int dim){
	cudaDeviceSynchronize();
	d_showVec<<<1,1>>>(show, dim);
	cudaDeviceSynchronize();
}

void d_showVecf(int* show, int dim){
	cudaDeviceSynchronize();
	d_showVecI<<<1,1>>>(show, dim);
	cudaDeviceSynchronize();
}

void inline d_showVecf(complex_type* show, int dim){
	cudaDeviceSynchronize();
	d_showVecC<<<1,1>>>(show, dim);
	cudaDeviceSynchronize();
}

__global__ void d_showMat(float_type* show, int dim){
	printf("\n[");

	for (int j = 0; j < dim; j++){
    	for (int i = 0; i < dim; i++){
    		printf("%5.1f ", show[j+i*dim]);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

void d_showMatf(float_type* show, int dim){
	d_showMat<<<1,1>>>(show, dim);
//	cudaDeviceSynchronize();
}

__global__ void d_showMat(const complex_type* show, const int dim){
	printf("\n[");

	for (int j = 0; j < dim; j++){
    	for (int i = 0; i < dim; i++){
    		printf("(%7.3f + j*%7.3f) ", show[j+i*dim].x ,show[j+i*dim].y);
	    }
		printf("\n ");
	}
	printf("]\n\n");
}

void d_showMatf(const complex_type* show, const int dim){
	d_showMat<<<1,1>>>(show, dim);
//	cudaDeviceSynchronize();
}
