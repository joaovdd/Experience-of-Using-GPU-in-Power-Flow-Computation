#pragma once
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

inline int _ConvertSMVer2Cores(int major, int minor) {
	typedef struct {
		int SM;  
		
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] = {
		{0x30, 192},
		{0x32, 192},
		{0x35, 192},
		{0x37, 192},
		{0x50, 128},
		{0x52, 128},
		{0x53, 128},
		{0x60,  64},
		{0x61, 128},
		{0x62, 128},
		{0x70,  64},
		{0x72,  64},
		{0x75,  64},
		{0x80,  64},
		{0x86, 128},
		{-1, -1} };

	int index = 0;

	while (nGpuArchCoresPerSM[index].SM != -1) {
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
			return nGpuArchCoresPerSM[index].Cores;
		}

		index++;
	}

	printf(
		"MapSMtoCores for SM %d.%d is undefined."
		"  Default to use %d Cores/SM\n",
		major, minor, nGpuArchCoresPerSM[index - 1].Cores);
	return nGpuArchCoresPerSM[index - 1].Cores;
}

inline int gpuGetMaxGflopsDeviceId()
{
	int current_device = 0, sm_per_multiproc = 0;
	int max_compute_perf = 0, max_perf_device = 0;
	int device_count = 0, best_SM_arch = 0;
	cudaDeviceProp deviceProp;
	cudaGetDeviceCount(&device_count);

	while (current_device < device_count)
	{
		cudaGetDeviceProperties(&deviceProp, current_device);
		if (deviceProp.major > 0 && deviceProp.major < 9999)
		{
			best_SM_arch = MAX(best_SM_arch, deviceProp.major);
		}
		current_device++;
	}

	current_device = 0;
	while (current_device < device_count)
	{
		cudaGetDeviceProperties(&deviceProp, current_device);
		if (deviceProp.major == 9999 && deviceProp.minor == 9999)
		{
			sm_per_multiproc = 1;
		}
		else
		{
			sm_per_multiproc = _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
		}

		int compute_perf = deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;

		if (compute_perf > max_compute_perf)
		{
			
			if (best_SM_arch > 2)
			{
				
				if (deviceProp.major == best_SM_arch)
				{
					max_compute_perf = compute_perf;
					max_perf_device = current_device;
				}
			}
			else
			{
				max_compute_perf = compute_perf;
				max_perf_device = current_device;
			}
		}
		++current_device;
	}
	return max_perf_device;
}