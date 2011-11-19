#include "defines.h"
#include <cuda.h>

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
	blockIdx.x*gridDim.x+threadIdx.x,\
	blockIdx.y*gridDim.y+threadIdx.y,\
	blockIdx.z*gridDim.z+threadIdx.z,\
	__VA_ARGS__)

__device__ int device_is_active_point(int i, int j, int k, localN locN, int rank, parts_sizes parts);
