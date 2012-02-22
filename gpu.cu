#include "gpu.h"
#include "shared_test.cu"

__constant__ consts gpu_def [1];

#ifdef TWO_PHASE
#include "Two-phase/two-phase.cu"
#endif
#ifdef THREE_PHASE
#include "Three-phase/three-phase.cu"
#endif
#ifdef B_L
#include "Backley-Leverett/b-l.cu"
#endif


// Преобразование локальных координат процессора к глобальным
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа
// (если 2 соседа с обеих сторон,то +2 точки).
// Глобальные границы хранятся как обычные точки (отсюда и условие на rank==0)
__device__ int device_local_to_global(int local_index, char axis, consts def)
{
	int global_index = local_index;
	switch (axis)
	{
	case 'x':
	{
		global_index += (*gpu_def).rankx * (*gpu_def).Nx / (*gpu_def).sizex + min((*gpu_def).rankx, (*gpu_def).Nx % (*gpu_def).sizex);
		break;
	}
	case 'y':
	{
		global_index += (*gpu_def).ranky * (*gpu_def).Ny / (*gpu_def).sizey + min((*gpu_def).ranky, (*gpu_def).Ny % (*gpu_def).sizey);
		break;
	}
	case 'z':
	{
		global_index += (*gpu_def).rankz * (*gpu_def).Nz / (*gpu_def).sizez + min((*gpu_def).rankz, (*gpu_def).Nz % (*gpu_def).sizez);
		break;
	}
	default:
	{
		;
	} //{printf("Error!");}
	}
	//some_test(global_index);
	return global_index;
}

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
__device__ int device_is_active_point(int i, int j, int k, consts def)
{
	if (((*gpu_def).rank % (*gpu_def).sizex != 0 && i == 0) || ((*gpu_def).rank % (*gpu_def).sizex != (*gpu_def).sizex - 1 && i == (*gpu_def).locNx - 1)
	    || (((*gpu_def).rank % (((*gpu_def).sizex) * ((*gpu_def).sizey))) / (*gpu_def).sizex != 0 && j == 0)	|| (((*gpu_def).rank % (((*gpu_def).sizex) * ((*gpu_def).sizey))) / (*gpu_def).sizex != (*gpu_def).sizey - 1 && j == (*gpu_def).locNy - 1)
	    || ((((*gpu_def).rank / ((*gpu_def).sizex) / ((*gpu_def).sizey) != 0 && k == 0) || ((*gpu_def).rank / ((*gpu_def).sizex) / ((*gpu_def).sizey) == (*gpu_def).sizez - 1 && k == (*gpu_def).locNz - 1)) && (*gpu_def).sizez > 1))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

// Функция вычисления "эффективной" плотности
__device__ double cu_ro_eff_gdy(ptr_Arrays DevArraysPtr, int i, int j, int k, consts def)
{
	int media = DevArraysPtr.media[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

#ifdef THREE_PHASE
	double ro_g_dy = (DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1. - DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])
	                  + DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
	                  + DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) * ((*gpu_def).m[media]) * ((*gpu_def).g_const) * ((*gpu_def).hy);
#else
	double ro_g_dy = (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
	                  + DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) * ((*gpu_def).m[media]) * ((*gpu_def).g_const) * ((*gpu_def).hy);
#endif
	return ro_g_dy;
}

// Расчет плотностей, давления NAPL P2 и Xi во всех точках сетки
void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_ro_Pn_Xi_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign Pn, Xi, ro", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Расчет давления воды P1 и насыщенности NAPL S2 во всех точках сетки
void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int w = 1; w <= def.newton_iterations; w++)
	{
		Newton_method_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	}
	checkErrors("assign Pw and Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Расчет скорости в каждой точке сетки
__global__ void assign_u_kernel(ptr_Arrays DevArraysPtr, consts def)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx)) && (j < ((*gpu_def).locNy)) && (k < ((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def) == 1))
	{
		//CUPRINTF("assign u\n");
		double Xi_w = DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double Xi_n = DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double P_w = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double P_n = DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
#ifdef THREE_PHASE
		double Xi_g = DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double P_g = DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
#endif
		if (((*gpu_def).Nx) > 2)
		{
			if (i == 0)
			{
				DevArraysPtr.ux_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (DevArraysPtr.P_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_w) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (DevArraysPtr.P_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_n) / ((*gpu_def).hx);
				//CUPRINTF("assign u=%e\n",DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]);
#ifdef THREE_PHASE
				DevArraysPtr.ux_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (DevArraysPtr.P_g[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_g) / ((*gpu_def).hx);
#endif
			}
			if (i == ((*gpu_def).locNx) - 1)
			{
				DevArraysPtr.ux_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (P_w - DevArraysPtr.P_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (P_n - DevArraysPtr.P_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hx);
#ifdef THREE_PHASE
				DevArraysPtr.ux_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (P_g - DevArraysPtr.P_g[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hx);
#endif
			}
			if ((i != 0) && (i != ((*gpu_def).locNx) - 1))
			{
				DevArraysPtr.ux_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * ((DevArraysPtr.P_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx)));
				DevArraysPtr.ux_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * ((DevArraysPtr.P_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx)));
#ifdef THREE_PHASE
				DevArraysPtr.ux_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * ((DevArraysPtr.P_g[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_g[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx)));
#endif
			}
		}
		else
		{
			DevArraysPtr.ux_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
			DevArraysPtr.ux_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#ifdef THREE_PHASE
			DevArraysPtr.ux_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#endif
		}

		if (((*gpu_def).Ny) > 2)
		{
			if (j == 0)
			{
				DevArraysPtr.uy_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (
				            (DevArraysPtr.P_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (P_w - DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (
				            (DevArraysPtr.P_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (P_n - DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;
#ifdef THREE_PHASE
				DevArraysPtr.uy_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (
				            (DevArraysPtr.P_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (P_g - DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;
#endif
			}

			if (j == ((*gpu_def).locNy) - 1)
			{
				DevArraysPtr.uy_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (
				            (P_w - DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)
				            - (DevArraysPtr.P_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (
				            (P_n - DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)
				            - (DevArraysPtr.P_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (*gpu_def).hy;
#ifdef THREE_PHASE
				DevArraysPtr.uy_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (
				            (P_g - DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)
				            - (DevArraysPtr.P_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (*gpu_def).hy;
#endif
			}
			if ((j != 0) && (j != ((*gpu_def).locNy) - 1))
			{
				DevArraysPtr.uy_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (
				            (DevArraysPtr.P_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (DevArraysPtr.P_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (2 * (*gpu_def).hy);

				DevArraysPtr.uy_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (
				            (DevArraysPtr.P_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (DevArraysPtr.P_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (2 * (*gpu_def).hy);
#ifdef THREE_PHASE
				DevArraysPtr.uy_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (
				            (DevArraysPtr.P_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j + 1) * (*gpu_def).hy)
				            - (DevArraysPtr.P_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] - DevArraysPtr.ro_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] * (*gpu_def).g_const * (j - 1) * (*gpu_def).hy)) / (2 * (*gpu_def).hy);
#endif
			}
		}
		else
		{
			DevArraysPtr.uy_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
			DevArraysPtr.uy_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#ifdef THREE_PHASE
			DevArraysPtr.uy_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#endif
		}

		if (((*gpu_def).Nz) > 2)
		{
			if (k == 0)
			{
				DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_w) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_n) / ((*gpu_def).hz);
#ifdef THREE_PHASE
				DevArraysPtr.uz_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P_g) / ((*gpu_def).hz);
#endif
			}
			if (k == ((*gpu_def).Nz) - 1)
			{
				DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * (P_w - DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * (P_n - DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hz);
#ifdef THREE_PHASE
				DevArraysPtr.uz_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * (P_g - DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / ((*gpu_def).hz);
#endif
			}
			if ((k != 0) && (k != ((*gpu_def).Nz) - 1))
			{
				DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_w * ((DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hz)));
				DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_n * ((DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hz)));
#ifdef THREE_PHASE
				DevArraysPtr.uz_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = Xi_g * ((DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hz)));
#endif
			}
		}
		else
		{
			DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
			DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#ifdef THREE_PHASE
			DevArraysPtr.uz_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0;
#endif
		}

		device_test_nan(DevArraysPtr.ux_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.ux_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uy_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uy_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
#ifdef THREE_PHASE
		device_test_nan(DevArraysPtr.ux_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uy_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.uz_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
#endif
	}
}

// Расчет скоростей во всех точках сетки
void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_u_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign u", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Расчет ro*S в каждой точке сетки методом направленных разностей
__global__ void assign_roS_kernel_nr(ptr_Arrays DevArraysPtr, double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx) - 1) && (j < (*gpu_def).locNy - 1) && (k < ((*gpu_def).locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != ((*gpu_def).locNz) - 1)) || (((*gpu_def).locNz) < 2)))
	{
		int media = DevArraysPtr.media[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

		double q_w = 0;
		double q_n = 0;

#ifdef B_L
		/*
		double F_bl=0;
		// В левом нижнем углу нагнетающая скважина
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			q_w=(*gpu_def).Q;
			q_n=0;
		}

		// В правом верхнем углу добывающая скважина
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			double S_e = (1. - DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]);
			double k_w = pow(S_e, (2. + 3. * ((*gpu_def).lambda[media])) / (*gpu_def).lambda[media]);
			double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));

			F_bl = (k_w/((*gpu_def).mu_w)) / (k_w/((*gpu_def).mu_w) + k_n/((*gpu_def).mu_n));
			q_w=-1 * (*gpu_def).Q * F_bl;
			q_n=-1 * (*gpu_def).Q * (1-F_bl);
		}
		*/
#endif

		double S2 = DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double roS2 = DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * S2;
		double P1 = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double P2 = DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
#ifdef THREE_PHASE
		double S1 = DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double roS1 = DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * S1;
		double roS3 = DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - S1 - S2);
		double P3 = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
#else
		double roS1 = DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - S2);
#endif
		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3;
#ifdef THREE_PHASE
		double b1, b2, b3;
#endif

		if (((*gpu_def).Nz) < 2)
		{
			f3 = 0;
			g3 = 0;
#ifdef THREE_PHASE
			b3 = 0;
#endif
		}
		else
		{
			z2 = -(DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P1) / (*gpu_def).hz;
			z1 = -(P1 - DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;

			f3 = (((z2 + abs(z2)) / 2.0 - (z1 - abs(z1)) / 2.0) * (-1) * DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
			      (z1 + abs(z1)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
			      (z2 - abs(z2)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;

			z2 = -(DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P2) / (*gpu_def).hz;
			z1 = -(P2 - DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;

			g3 = (((z2 + abs(z2)) / 2.0 - (z1 - abs(z1)) / 2.0) * (-1) * DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
			      (z1 + abs(z1)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
			      (z2 - abs(z2)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;
#ifdef THREE_PHASE
			z2 = -(DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P3) / (*gpu_def).hz;
			z1 = -(P3 - DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;

			b3 = (((z2 + abs(z2)) / 2.0 - (z1 - abs(z1)) / 2.0) * (-1) * DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
			      (z1 + abs(z1)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
			      (z2 - abs(z2)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hz;
#endif
		}

		x2 = -(DevArraysPtr.P_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P1) / (*gpu_def).hx;
		x1 = -(P1 - DevArraysPtr.P_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		y2 = -(DevArraysPtr.P_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P1) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j + 1) - DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j) / (*gpu_def).hy;
		y1 = -(P1 - DevArraysPtr.P_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j - DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j - 1)) / (*gpu_def).hy;

		f1 = (((x2 + abs(x2)) / 2.0 - (x1 - abs(x1)) / 2.0) * (-1) * DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (x1 + abs(x1)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (x2 - abs(x2)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		f2 = (((y2 + abs(y2)) / 2.0 - (y1 - abs(y1)) / 2.0) * (-1) * DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (y1 + abs(y1)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (y2 - abs(y2)) / 2.0 * (-1) * DevArraysPtr.Xi_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy;


		x2 = -(DevArraysPtr.P_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P2) / (*gpu_def).hx;
		x1 = -(P2 - DevArraysPtr.P_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		y2 = -(DevArraysPtr.P_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P2) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j + 1) - DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j) / (*gpu_def).hy;
		y1 = -(P2 - DevArraysPtr.P_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j - DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j - 1)) / (*gpu_def).hy;

		g1 = (((x2 + abs(x2)) / 2.0 - (x1 - abs(x1)) / 2.0) * (-1) * DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (x1 + abs(x1)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (x2 - abs(x2)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		g2 = (((y2 + abs(y2)) / 2.0 - (y1 - abs(y1)) / 2.0) * (-1) * DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (y1 + abs(y1)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (y2 - abs(y2)) / 2.0 * (-1) * DevArraysPtr.Xi_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy;

		DevArraysPtr.roS_w_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS1;
		DevArraysPtr.roS_n_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS2;
		DevArraysPtr.roS_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS1 - ((*gpu_def).dt / (*gpu_def).m[media]) * (q_w + f1 + f2 + f3);
		DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS2 - ((*gpu_def).dt / (*gpu_def).m[media]) * (q_n + g1 + g2 + g3);

		device_test_positive(DevArraysPtr.roS_w_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
#ifdef THREE_PHASE
		x2 = -(DevArraysPtr.P_g[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P3) / (*gpu_def).hx;
		x1 = -(P3 - DevArraysPtr.P_g[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		y2 = -(DevArraysPtr.P_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - P3) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j + 1) - DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j) / (*gpu_def).hy;
		y1 = -(P3 - DevArraysPtr.P_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * j - DevArraysPtr.ro_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (j - 1)) / (*gpu_def).hy;

		b1 = (((x2 + abs(x2)) / 2.0 - (x1 - abs(x1)) / 2.0) * (-1) * DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (x1 + abs(x1)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (x2 - abs(x2)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hx;

		b2 = (((y2 + abs(y2)) / 2.0 - (y1 - abs(y1)) / 2.0) * (-1) * DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] -
		      (y1 + abs(y1)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] +
		      (y2 - abs(y2)) / 2.0 * (-1) * DevArraysPtr.Xi_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ro_g[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (*gpu_def).hy;

		DevArraysPtr.roS_g_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS3;
		DevArraysPtr.roS_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = roS3 - ((*gpu_def).dt / (*gpu_def).m[media]) * (b1 + b2 + b3);

		device_test_positive(DevArraysPtr.roS_g_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
#endif
	}
}

// Расчет ro*S в каждой точке сетки
__global__ void assign_roS_kernel(ptr_Arrays DevArraysPtr, double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx) - 1) && (j < (*gpu_def).locNy - 1) && (k < ((*gpu_def).locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != ((*gpu_def).locNz))) || (((*gpu_def).locNz) < 2)))
	{
		int local = i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy);
		int media = DevArraysPtr.media[local];
		double S_n = DevArraysPtr.S_n[local];
		double roS_w = DevArraysPtr.ro_w[local] * (1 - S_n);
		double roS_n = DevArraysPtr.ro_n[local] * S_n;

		double divgrad1, divgrad2, Tx1, Ty1, Tx2, Ty2, Tz1, Tz2, A1 = 0, A2 = 0;

		if (((*gpu_def).Nz) < 2)
		{
			divgrad1 = 0;
			divgrad2 = 0;
			Tz1 = 0;
			Tz2 = 0;
		}
		else
		{
			divgrad1 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) * (DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1. - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[local] * (1. - S_n) + DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1. - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			divgrad2 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) * (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[local] * S_n + DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			Tz1 = (DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uz_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2. * ((*gpu_def).hz));
			Tz2 = (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uz_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2. * ((*gpu_def).hz));
		}

		divgrad1 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) *
		            ((DevArraysPtr.ro_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		             (DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		divgrad2 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) *
		            ((DevArraysPtr.ro_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * S_n + DevArraysPtr.ro_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		             (DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * S_n + DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.S_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		Tx1 = (DevArraysPtr.ro_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ux_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ux_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx));
		Ty1 = (DevArraysPtr.ro_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uy_w[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uy_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hy));
		Tx2 = (DevArraysPtr.ro_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ux_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.ux_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx));
		Ty2 = (DevArraysPtr.ro_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uy_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.ro_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * DevArraysPtr.uy_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) / (2 * ((*gpu_def).hy));

		double q_w = 0;
		double q_n = 0;

#ifdef B_L
		double F_bl = 0;
		// В левом нижнем углу нагнетающая скважина
		if (((i == 0) && (j == (*gpu_def).Ny - 2)) || ((i == 1) && (j == (*gpu_def).Ny - 1)))
		{
			q_w = (*gpu_def).Q;
			q_n = 0;
		}

		// В правом верхнем углу добывающая скважина
		if (((i == 0) && (j == (*gpu_def).Ny - 2)) || ((i == 1) && (j == (*gpu_def).Ny - 1)))
		{
			//int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			double S_e = (1. - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]);
			double k_w = pow(S_e, (2. + 3. * ((*gpu_def).lambda[media])) / (*gpu_def).lambda[media]);
			double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));

			F_bl = (k_w / ((*gpu_def).mu_w)) / (k_w / ((*gpu_def).mu_w) + k_n / ((*gpu_def).mu_n));
			q_w = -1 * (*gpu_def).Q * F_bl;
			q_n = -1 * (*gpu_def).Q * (1 - F_bl);
		}
#endif

		if ((t < 2 * ((*gpu_def).dt)) || TWO_LAYERS)
		{
			A1 = roS_w + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = roS_n + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1 + (2 * roS_w * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_w_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
			A2 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2 + (2 * roS_n * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_n_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
		}

		DevArraysPtr.roS_w_old[local] = roS_w;
		DevArraysPtr.roS_n_old[local] = roS_n;
		DevArraysPtr.roS_w[local] = A1;
		DevArraysPtr.roS_n[local] = A2;

		device_test_positive(DevArraysPtr.roS_w_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// Расчет ro*S во всех точках сетки
void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def)
{
#ifdef NR
	assign_roS_kernel_nr <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, t);
#else
	assign_roS_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, t);
#endif
	checkErrors("assign roS", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Граничные условия на S2
__global__ void Sn_boundary_kernel(ptr_Arrays DevArraysPtr, consts def)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx)) && (j < ((*gpu_def).locNy)) && (k < ((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def) == 1))
	{
		if ((i == 0) && ((*gpu_def).Nx > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			return;
		}

		if ((i == ((*gpu_def).locNx) - 1) && (((*gpu_def).Nx) > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			return;
		}

		if ((j == ((*gpu_def).locNy) - 1) && (((*gpu_def).locNy) > 2))
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

			if (i == 0)
			{
				DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			}
			if (i == (*gpu_def).Nx - 1)
			{
				DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			}
			return;
		}

		if ((j == 0) && (((*gpu_def).locNy) > 2))
		{
			int I = device_local_to_global(i, 'x', gpu_def[0]);
			if ((I >= ((*gpu_def).Nx) / 2 - ((*gpu_def).source)) && (I <= ((*gpu_def).Nx) / 2 + ((*gpu_def).source)) && (k >= ((*gpu_def).Nz) / 2 - ((*gpu_def).source)) && (k <= ((*gpu_def).Nz) / 2 + ((*gpu_def).source)))
			{
				DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).S_n_gr;
			}
			else
			{
				DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + (j + 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			}
			return;
		}

		if ((k == 0) && ((*gpu_def).Nz > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			return;
		}

		if ((k == (*gpu_def).Nz - 1) && ((*gpu_def).Nz > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			return;
		}

		device_test_S(DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// Граничные условия на P1
__global__ void Pw_boundary_kernel(ptr_Arrays DevArraysPtr, consts def)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx)) && (j < ((*gpu_def).locNy)) && (k < ((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def) == 1))
	{
		if ((i == 0) && ((*gpu_def).Nx > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			//return;
		}

		if ((i == ((*gpu_def).locNx) - 1) && ((*gpu_def).Nx > 2) && (j > 0) && (j < (*gpu_def).locNy - 1))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			//return;
		}

		if ((j == ((*gpu_def).locNy) - 1) && (((*gpu_def).locNy) > 2))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + (j - 1) * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * (*gpu_def).locNy] + cu_ro_eff_gdy(DevArraysPtr, i, j - 1, k, def);
			if (i == 0)
			{
				DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			}
			if (i == (*gpu_def).Nx - 1)
			{
				DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i - 1 + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			}
			//return;
		}

		if ((j == 0) && (((*gpu_def).locNy) > 2))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).P_atm;
			//return;
		}

		if ((k == 0) && (((*gpu_def).Nz) > 2) && (j > 0) && (j < ((*gpu_def).locNy) - 1))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k + 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			//return;
		}

		if ((k == ((*gpu_def).locNz) - 1) && (((*gpu_def).locNz) > 2) && (j > 0) && (j < ((*gpu_def).locNy) - 1))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + (k - 1) * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			//return;
		}

		device_test_positive(DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// Применение граничных условий
void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
#ifndef THREE_PHASE
	Sn_boundary_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	Pw_boundary_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign Pw", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
#else
	Border_S_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign S", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	Border_P_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr, def);
	checkErrors("assign Pw", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
#endif
}

// Функция загрузки данных в память хоста
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy(HostArrayPtr, DevArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy double data to device", __FILE__, __LINE__);
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, consts def)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(int), cudaMemcpyHostToDevice);
	checkErrors("copy int data to device", __FILE__, __LINE__);
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, consts def)
{
	cudaMalloc((void**) DevBuffer,  2 * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).P_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).P_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).S_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ro_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ro_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ux_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uy_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uz_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ux_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uy_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uz_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).Xi_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).Xi_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_n), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_w_old), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_n_old), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).media), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(int));

#ifdef THREE_PHASE
	cudaMalloc((void**) & ((*ArraysPtr).P_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).S_w), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ro_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).ux_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uy_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).uz_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).Xi_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_g), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).roS_g_old), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
#endif

#ifdef B_L
	cudaMalloc((void**) & ((*ArraysPtr).K), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
#endif

	checkErrors("memory allocation", __FILE__, __LINE__);
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	cudaFree(DevBuffer);
	cudaFree(DevArraysPtr.P_w);
	cudaFree(DevArraysPtr.P_n);
	cudaFree(DevArraysPtr.S_n);
	cudaFree(DevArraysPtr.ro_w);
	cudaFree(DevArraysPtr.ro_n);
	cudaFree(DevArraysPtr.ux_w);
	cudaFree(DevArraysPtr.uy_w);
	cudaFree(DevArraysPtr.uz_w);
	cudaFree(DevArraysPtr.ux_n);
	cudaFree(DevArraysPtr.uy_n);
	cudaFree(DevArraysPtr.uz_n);
	cudaFree(DevArraysPtr.Xi_w);
	cudaFree(DevArraysPtr.Xi_n);
	cudaFree(DevArraysPtr.roS_w);
	cudaFree(DevArraysPtr.roS_n);
	cudaFree(DevArraysPtr.roS_w_old);
	cudaFree(DevArraysPtr.roS_n_old);
	cudaFree(DevArraysPtr.media);

#ifdef THREE_PHASE
	cudaFree(DevArraysPtr.P_g);
	cudaFree(DevArraysPtr.S_w);
	cudaFree(DevArraysPtr.ro_g);
	cudaFree(DevArraysPtr.ux_g);
	cudaFree(DevArraysPtr.uy_g);
	cudaFree(DevArraysPtr.uz_g);
	cudaFree(DevArraysPtr.Xi_g);
	cudaFree(DevArraysPtr.roS_g);
	cudaFree(DevArraysPtr.roS_g_old);
#endif

#ifdef B_L
	cudaFree(DevArraysPtr.K);
#endif

	checkErrors("memory reliase", __FILE__, __LINE__);
}

// Инициализация ускорителя
// Расчет происходит на ускорителе, номер которого равен
// номеру запускающего процессора
void device_initialization(consts* def)
{
	// Было бы очень неплохо вместо GPU_PER_NODE использовать cudaGetDeviceCount
	//int deviceCount;
	//cudaGetDeviceCount ( &deviceCount );

	// Считаем, что ядер на узле не меньше, чем ускорителей
	int device = (*def).rank % GPU_PER_NODE;
	cudaSetDevice(device);

	// Количество запускаемых блоков
	// Если число точек сетки не кратно размеру блока,
	// то количество блоков будет на 1 больше.
	(*def).blocksX = (*def).locNx / BlockNX;
	if (((*def).locNx % BlockNX) != 0)
	{
		((*def).blocksX)++;
	}
	(*def).blocksY = ((*def).locNy) / BlockNY;
	if (((*def).locNy) % BlockNY != 0)
	{
		((*def).blocksY)++;
	}
	(*def).blocksZ = ((*def).locNz) / BlockNZ;
	if (((*def).locNz) % BlockNZ != 0)
	{
		((*def).blocksZ)++;
	}

	consts* deff = new consts[1];
	deff[0] = (*def);
	cudaMemcpyToSymbol(gpu_def, deff, sizeof(consts));
	checkErrors("constant memory copy", __FILE__, __LINE__);

	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, device);

	if (devProp.major < 2)
	{
		printf("\nError! Compute capability < 2, rank=%d\n", (*def).rank);
	}

	if (!(*def).rank)
	{
		//printf ( "Device %d\n", device );
		printf("Name : %s\n", devProp.name);
		printf("Compute capability : %d.%d\n", devProp.major, devProp.minor);
		printf("Total Global Memory : %ld\n", devProp.totalGlobalMem);
		printf("Shared memory per block: %d\n", devProp.sharedMemPerBlock);
		printf("Registers per block : %d\n", devProp.regsPerBlock);
		printf("Warp size : %d\n", devProp.warpSize);
		printf("Max threads per block : %d\n", devProp.maxThreadsPerBlock);
		printf("Total constant memory : %d\n", devProp.totalConstMem);
		printf("Number of multiprocessors: %d\n",  devProp.multiProcessorCount);
		printf("Kernel execution timeout: %s\n\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
		}
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
		}


		// Максимальный размер расчетной сетки для ускорителя
		// sizeof(ptr_Arrays)/4 - количество параметров в точке, т.к. 4 -размер одного указателя
		printf("\nTotal NAPL_Filtration grid size : %d\n\n", devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4));
	}

	// (def.locNx)+2 потому что 2NyNz на буфер обмена выделяется
	if (((*def).locNx + 2) * ((*def).locNy) * ((*def).locNz) > (devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4)))
	{
		printf("\nError! Not enough memory at GPU, rank=%d\n", (*def).rank);
	}
	fflush(stdout);

	// Инициализируем библиотеку cuPrintf для вывода текста на консоль
	// прямо из kernel
	cudaPrintfInit();
}

// Финализация ускорителя
void device__finalization(void)
{
	// Останавливаем библиотеку cuPrintf для вывода текста на консоль
	// прямо из kernel
	cudaPrintfEnd();
}

__global__ void load_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < (*gpu_def).locNy && k < ((*gpu_def).locNz))
	{
		DevBuffer[j + ((*gpu_def).locNy)*k] = DevArrayPtr[1 + ((*gpu_def).locNx) * j + ((*gpu_def).locNx) * ((*gpu_def).locNy) * k];
		DevBuffer[j + ((*gpu_def).locNy)*k + ((*gpu_def).locNy) * ((*gpu_def).locNz)] = DevArrayPtr[((*gpu_def).locNx) - 2 + ((*gpu_def).locNx) * j + ((*gpu_def).locNx) * ((*gpu_def).locNy) * k];

		device_test_nan(DevBuffer[j + ((*gpu_def).locNy)*k], __FILE__, __LINE__);
		device_test_nan(DevBuffer[j + ((*gpu_def).locNy)*k + ((*gpu_def).locNy) * ((*gpu_def).locNz)], __FILE__, __LINE__);
	}
}


void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, 2 * (def.Ny) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void save_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < ((*gpu_def).locNy) && k < ((*gpu_def).locNz))
	{
		if ((*gpu_def).rank != 0)
		{
			DevArrayPtr[((*gpu_def).locNx)*j + ((*gpu_def).locNx) * ((*gpu_def).locNy)*k] = DevBuffer[j + ((*gpu_def).locNy) * k];
			device_test_nan(DevArrayPtr[((*gpu_def).locNx)*j + ((*gpu_def).locNx) * ((*gpu_def).locNy)*k], __FILE__, __LINE__);
		}
		if ((*gpu_def).rank != (*gpu_def).sizex - 1)
		{
			DevArrayPtr[((*gpu_def).locNx) - 1 + ((*gpu_def).locNx)*j + ((*gpu_def).locNx) * ((*gpu_def).locNy)*k] = DevBuffer[j + ((*gpu_def).locNy) * k + ((*gpu_def).locNy) * ((*gpu_def).locNz)];
			device_test_nan(DevArrayPtr[((*gpu_def).locNx) - 1 + ((*gpu_def).locNx)*j + ((*gpu_def).locNx) * ((*gpu_def).locNy)*k], __FILE__, __LINE__);
		}
	}
}

void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, 2 * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

