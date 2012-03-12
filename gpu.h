#ifndef GPU_H
#define GPU_H

#include "defines.h"
#include <cuda.h>

#define CUPRINTF(fmt, ...) printf("[%d, %d , %d]:\t" fmt, \
                                  blockIdx.x*gridDim.x+threadIdx.x,\
                                  blockIdx.y*gridDim.y+threadIdx.y,\
                                  blockIdx.z*gridDim.z+threadIdx.z,\
                                  __VA_ARGS__)

__device__ int device_is_active_point(int i, int j, int k);
__device__ int device_local_to_global(int local_index, char axis);
__device__ double cu_ro_eff_gdy(ptr_Arrays DevArraysPtr, int i, int j, int k);

#endif
