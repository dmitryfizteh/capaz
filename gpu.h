#include "defines.h"
#include <cuda.h>

void gpu_test_nan(double* HostArrayPtr, double* DevArrayPtr, int rank, int size, int localNx, consts def, char* file, int line);
void checkErrors(char *label);

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
	blockIdx.x*gridDim.x+threadIdx.x,\
	blockIdx.y*gridDim.y+threadIdx.y,\
	blockIdx.z*gridDim.z+threadIdx.z,\
	__VA_ARGS__)


