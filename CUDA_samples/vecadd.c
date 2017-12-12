// Kernel definition, see also section 4.2.3 of Nvidia Cuda Programming Guide
__global__ void vecAdd(float* A, float* B, float* C) {
	//threadIdx.x is a build-in variable provided by CUDA runtime
	int i = threadIdx.x;
	A[i] = 0;
	B[i] = 0;
	C[i] = A[i] + B[i];
}

#include <stdio.h>
#define SIZE 10

int main() {
	int N=SIZE;
	float A[SIZE], B[SIZE], C[SIZE];
	float *devPtrA;
	float *devPtrB;
	float *devPtrC;
	int memsize = SIZE * sizeof(float);

	cudaMalloc((void**)&devPtrA, memsize);
	cudaMalloc((void**)&devPtrB, memsize);
	cudaMalloc((void**)&devPtrC, memsize);
	cudaMemcpy(devPtrA, A, memsize, cudaMemcpyHostToDevice);
	cudaMemcpy(devPtrB, B, memsize, cudaMemcpyHostToDevice);
	// __global__ funcitons are called: Func<<< Dg, Db, Ns >>>(paramter);
	vecAdd<<<1, N>>>(devPtrA, devPtrB, devPtrC);
	cudaMemcpy(C, devPtrC, memsize, cudaMemcpyDeviceToHost);

	for (int i = 0; i < SIZE; i++) printf("C[%d]=%f\n",i,C[i]);

	cudaFree(devPtrA);
	cudaFree(devPtrB);
	cudaFree(devPtrC);

}
