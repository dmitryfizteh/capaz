#include "defines.h"
#include <cuda.h>

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
	blockIdx.x*gridDim.x+threadIdx.x,\
	blockIdx.y*gridDim.y+threadIdx.y,\
	blockIdx.z*gridDim.z+threadIdx.z,\
	__VA_ARGS__)


