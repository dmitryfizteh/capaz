#include "../defines.h"
#include "../gpu.h"
#include "two-phase.h"

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
__device__ int device_is_active_point_local(int i, int localNx, int rank, int size)
{
	if((rank!=0 && i==0) || (rank!=size-1 && i==localNx-1))
		return 0;
	else
		return 1;
}

