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
__device__ int device_i_to_I(int i, int rank, parts_sizes parts, consts def)
{
	int I;
	if (rank <= ((*gpu_def).Nx)%parts.x)
	{
		if(rank==0)
			I=i;
		else
			I=(((*gpu_def).Nx)/parts.x+1)*rank+i-1;
	}
	else
		I=(((*gpu_def).Nx)/parts.x+1)*rank-(rank-((*gpu_def).Nx)%parts.x)+i-1;
	return I;
}

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
__device__ int device_is_active_point(int i, int j, int k, localN locN, int rank, parts_sizes parts)
{
	if((rank % parts.x != 0 && i == 0) || (rank % parts.x != parts.x - 1 && i == locN.x - 1)
		|| ((rank % ((parts.x) * (parts.y))) / parts.x != 0 && j == 0)	|| ((rank % ((parts.x) * (parts.y))) / parts.x != parts.y - 1 && j == locN.y - 1)
		|| (((rank / (parts.x) / (parts.y) != 0 && k == 0) || (rank / (parts.x) / (parts.y) == parts.z - 1 && k == locN.z - 1)) && parts.z > 1))
		return 0;
	else
		return 1;
}

// Расчет плотностей, давления NAPL P2 и Xi во всех точках сетки
void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	assign_ro_Pn_Xi_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr,locN,rank,parts); 
	checkErrors("assign Pn, Xi, ro", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Расчет давления воды P1 и насыщенности NAPL S2 во всех точках сетки
void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	for (int w=1;w<=def.newton_iterations;w++)
	{
		Newton_method_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, locN); 
		checkErrors("assign Pw and Sn", __FILE__, __LINE__);
		cudaPrintfDisplay(stdout, true);
	}
}

// Расчет скорости в каждой точке сетки
__global__ void assign_u_kernel(ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)) && (j<(locN.y)) && (k<(locN.z)) && (device_is_active_point(i, j, k, locN, rank, parts)==1))
	{
		//CUPRINTF("assign u\n");
		double Xi_w = DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double Xi_n = DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double P_w = DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double P_n = DevArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		if (((*gpu_def).Nx)>2)
		{
			if (i == 0)
			{
				DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (DevArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - P_w) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (DevArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - P_n) / ((*gpu_def).hx);
				//CUPRINTF("assign u=%e\n",DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]);
			}
			if (i == (locN.x) - 1)
			{
				DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (P_w - DevArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (P_n - DevArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / ((*gpu_def).hx);
			}
			if ((i != 0) && (i != (locN.x) - 1))
			{
				DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * ((DevArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hx)));
				DevArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * ((DevArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hx)));
			}
		}
		else
		{
			DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			DevArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		}
	
		if (((*gpu_def).Ny)>2)
		{
			if (j == 0)
			{
				DevArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (
					(DevArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
				    - (P_w - DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (
					(DevArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (P_n - DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;			}

			if (j == (locN.y) - 1)
			{
				DevArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (
					(P_w - DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * j * (*gpu_def).hy)
					- (DevArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (
					(P_n - DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * j * (*gpu_def).hy)
					- (DevArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (*gpu_def).hy;
			}
			if ((j != 0) && (j != (locN.y) - 1))
			{
				DevArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (
					(DevArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (DevArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (2*(*gpu_def).hy);

				DevArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (
					(DevArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (DevArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] - DevArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (2*(*gpu_def).hy);
			}
		}
		else
		{
			DevArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			DevArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		}

		if (((*gpu_def).Nz)>2)
		{
			if (k == 0)
			{
				DevArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (DevArraysPtr.P_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - P_w) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (DevArraysPtr.P_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - P_n) / ((*gpu_def).hz);
			}
			if (k == ((*gpu_def).Nz) - 1)
			{
				DevArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * (P_w - DevArraysPtr.P_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * (P_n - DevArraysPtr.P_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / ((*gpu_def).hz);
			}
			if ((k != 0) && (k != ((*gpu_def).Nz) - 1))
			{
				DevArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_w * ((DevArraysPtr.P_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - DevArraysPtr.P_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hz)));
				DevArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = Xi_n * ((DevArraysPtr.P_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - DevArraysPtr.P_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hz)));
			}
		}
		else
		{
			DevArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			DevArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		}

	device_test_nan(DevArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}

// Расчет скоростей во всех точках сетки
void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def)
{
	assign_u_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, locN, rank, parts); 
	checkErrors("assign u", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// Расчет ro*S в каждой точке сетки методом направленных разностей
__global__ void assign_roS_kernel_nr(ptr_Arrays DevArraysPtr, localN locN, double t)
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;
	
	if ((i<(locN.x)-1) && (j<locN.y-1) && (k<(locN.z)) && (i!=0) && (j!=0) && (((k!=0) && (k!=(locN.z)-1)) || ((locN.z)<2)))
	{
		int media = DevArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		double S2 = DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double roS1 = DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - S2);
		double roS2 = DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * S2;
 		double P1 = DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double P2 = DevArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3;

		if (((*gpu_def).Nz)<2)
		{
			f3=0;
			g3=0;
		}
		else
		{
			z2 = -(DevArraysPtr.P_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - P1)/(*gpu_def).hz;
			z1 = -(P1 - DevArraysPtr.P_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])/(*gpu_def).hz;

			f3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                      (z1 + abs(z1))/2.0*(-1)* DevArraysPtr.Xi_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] +
                      (z2 - abs(z2))/2.0*(-1)* DevArraysPtr.Xi_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)])/(*gpu_def).hz;

			z2 = -(DevArraysPtr.P_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - P2)/(*gpu_def).hz;
			z1 = -(P2 - DevArraysPtr.P_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])/(*gpu_def).hz;

			g3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                      (z1 + abs(z1))/2.0*(-1)* DevArraysPtr.Xi_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] +
                      (z2 - abs(z2))/2.0*(-1)* DevArraysPtr.Xi_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)])/(*gpu_def).hz;
		}

		x2 = -(DevArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - P1)/(*gpu_def).hx;
        x1 = -(P1 - DevArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hx;

        y2 = -(DevArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - P1)/(*gpu_def).hy + DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).g_const;
        y1 = -(P1 - DevArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hy + DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).g_const;

        f1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (x1 + abs(x1))/2.0*(-1)* DevArraysPtr.Xi_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] +
                (x2 - abs(x2))/2.0*(-1)* DevArraysPtr.Xi_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hx;

        f2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (y1 + abs(y1))/2.0*(-1)* DevArraysPtr.Xi_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] +
                (y2 - abs(y2))/2.0*(-1)* DevArraysPtr.Xi_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hy;


        x2 = -(DevArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - P2)/(*gpu_def).hx;
        x1 = -(P2 - DevArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hx;

        y2 = -(DevArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - P2)/(*gpu_def).hy + DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).g_const;
        y1 = -(P2 - DevArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hy + DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).g_const;

        g1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (x1 + abs(x1))/2.0*(-1)* DevArraysPtr.Xi_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] +
                (x2 - abs(x2))/2.0*(-1)* DevArraysPtr.Xi_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hx;

        g2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (y1 + abs(y1))/2.0*(-1)* DevArraysPtr.Xi_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] +
                (y2 - abs(y2))/2.0*(-1)* DevArraysPtr.Xi_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])/(*gpu_def).hy;

		DevArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = roS1;
		DevArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = roS2;
		DevArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = roS1 - ((*gpu_def).dt/(*gpu_def).m[media])*(f1 + f2 + f3);
		DevArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = roS2 - ((*gpu_def).dt/(*gpu_def).m[media])*(g1 + g2 + g3);

		device_test_positive(DevArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}

// Расчет ro*S в каждой точке сетки
__global__ void assign_roS_kernel(ptr_Arrays DevArraysPtr, localN locN, double t) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)-1) && (j<locN.y-1) && (k<(locN.z)) && (i!=0) && (j!=0) && (((k!=0) && (k!=(locN.z))) || ((locN.z)<2)))
	{
		int local=i+j*(locN.x)+k*(locN.x)*(locN.y);
		int media = DevArraysPtr.media[local];
		double S_n = DevArraysPtr.S_n[local];
		double roS_w = DevArraysPtr.ro_w[local] * (1 - S_n);
		double roS_n = DevArraysPtr.ro_n[local] * S_n;

		double divgrad1, divgrad2, Tx1, Ty1, Tx2, Ty2, Tz1, Tz2, A1=0, A2=0;

		if (((*gpu_def).Nz)<2)
		{
			divgrad1=0;
			divgrad2=0;
			Tz1=0;
			Tz2=0;
		}
		else
		{
			divgrad1 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) * (DevArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * (1. - DevArraysPtr.S_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)]) - 2 * DevArraysPtr.ro_w[local] * (1. - S_n) + DevArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * (1. - DevArraysPtr.S_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			divgrad2 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) * (DevArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * DevArraysPtr.S_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - 2 * DevArraysPtr.ro_n[local] * S_n + DevArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * (DevArraysPtr.S_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			Tz1 = (DevArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * DevArraysPtr.uz_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - DevArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * DevArraysPtr.uz_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2. * ((*gpu_def).hz));
			Tz2 = (DevArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * DevArraysPtr.uz_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - DevArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * DevArraysPtr.uz_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2. * ((*gpu_def).hz));
		}

		divgrad1 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) *
		((DevArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - DevArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)]) - 2 * DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - S_n) + DevArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - DevArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		(DevArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (1 - DevArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]) - 2 * DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - S_n) + DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (1 - DevArraysPtr.S_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		divgrad2 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) *
		((DevArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - 2 * DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * S_n + DevArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		(DevArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - 2 * DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * S_n + DevArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.S_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		Tx1 = (DevArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ux_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ux_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hx));
		Ty1 = (DevArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.uy_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.uy_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hy));
		Tx2 = (DevArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ux_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.ux_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hx));
		Ty2 = (DevArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.uy_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - DevArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * DevArraysPtr.uy_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (2 * ((*gpu_def).hy));

		double q_w=0;
		double q_n=0;

#ifdef B_L
		// В левом нижнем углу нагнетающая скважина
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			q_w=0;
			//q_n=def.Q;
		}

		// В правом верхнем углу добывающая скважина
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			//q_w=def.Q * F_bl(S);
			//q_n=def.Q * (1-F_bl(S));
		}
#endif

		if ((t < 2 * ((*gpu_def).dt)) || TWO_LAYERS)
		{
			A1 = roS_w + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = roS_n + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1 + (2 * roS_w * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
			A2 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2 + (2 * roS_n * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
		}

		DevArraysPtr.roS_w_old[local] = roS_w;
		DevArraysPtr.roS_n_old[local] = roS_n;
		DevArraysPtr.roS_w[local] = A1;
		DevArraysPtr.roS_n[local] = A2;

		device_test_positive(DevArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}

// Расчет ro*S во всех точках сетки
void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	#ifdef NR
		assign_roS_kernel_nr<<<dim3(blocksX,blocksY*blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, locN, t); 
	#else
		assign_roS_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, locN, t);
	#endif
		checkErrors("assign roS", __FILE__, __LINE__);
		cudaPrintfDisplay(stdout, true);
}

// Граничные условия на S2
__global__ void Sn_boundary_kernel(ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)) && (j<(locN.y)) && (k<(locN.z)) && (device_is_active_point(i, j, k, locN, rank, parts)==1))
	{
		if ((i == 0) && ((*gpu_def).Nx>2) && (j>0) && (j<locN.y - 1))
		{
		   DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)];
		   return;
		}

		if ((i == (locN.x) - 1) && (((*gpu_def).Nx)>2) && (j>0) && (j<locN.y - 1))
		{
			DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)];
			 return;
		}

		if ((j == (locN.y) - 1) && ((locN.y)>2))
		{
			DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)];

			if (i==0) 
				DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)];
			if (i==(*gpu_def).Nx - 1)
				DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)];
			return;
		}
		
		if ((j==0) && ((locN.y)>2))
		{
			int I=device_i_to_I(i, rank, parts, gpu_def[0]);
			if ((I>=((*gpu_def).Nx)/2-((*gpu_def).source)) && (I<=((*gpu_def).Nx)/2+((*gpu_def).source)) && (k>=((*gpu_def).Nz)/2-((*gpu_def).source)) && (k<=((*gpu_def).Nz)/2+((*gpu_def).source)))
				DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = (*gpu_def).S_n_gr;
			else
				DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)];
			return;
		}

		if ((k == 0) && ((*gpu_def).Nz > 2) && (j>0) && (j<locN.y - 1))
		{
			DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)];
			return;
		}
		
		if ((k == (*gpu_def).Nz - 1) && ((*gpu_def).Nz > 2) && (j>0) && (j<locN.y - 1))
		{
			DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.S_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)];
			return;
		}

		device_test_S(DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}

// Граничные условия на P1
__global__ void Pw_boundary_kernel(ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)) && (j<(locN.y)) && (k<(locN.z)) && (device_is_active_point(i, j, k, locN, rank, parts)==1))
	{
		if ((i == 0) && ((*gpu_def).Nx > 2) && (j>0) && (j<locN.y - 1))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)]; 
			//return;
		}

		if ((i == (locN.x) - 1) && ((*gpu_def).Nx>2) && (j>0) && (j<locN.y - 1))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)];
			//return;
		}

		if ((j == (locN.y) - 1) && ((locN.y)>2))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] + DevArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*locN.y] * (*gpu_def).g_const * (*gpu_def).hy;
			if (i==0) 
				DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)];
			if (i==(*gpu_def).Nx - 1)
				DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)];
			//return;
		}

		if ((j==0) && ((locN.y)>2))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = (*gpu_def).P_atm;
			//return;
		}

		if ((k == 0) && (((*gpu_def).Nz)>2) && (j>0) && (j<(locN.y) - 1))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)]; 
			//return;
		}

		if ((k == (locN.z) - 1) && ((locN.z)>2) && (j>0) && (j<(locN.y) - 1))
		{
			DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = DevArraysPtr.P_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)];
			//return;
		}

		device_test_positive(DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}

// Применение граничных условий
void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def)
{
	Sn_boundary_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr,locN, rank,parts); 
	checkErrors("assign Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	Pw_boundary_kernel<<<dim3(blocksX,blocksY,blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr,locN, rank,parts); 
	checkErrors("assign Pw", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}	

// Функция загрузки данных в память хоста
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def)
{
	cudaMemcpy( HostArrayPtr, DevArrayPtr, (locN.x)*(locN.y)*(locN.z)*sizeof(double), cudaMemcpyDeviceToHost );
	checkErrors("copy data to host", __FILE__, __LINE__);
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def)
{
	cudaMemcpy( DevArrayPtr, HostArrayPtr, (locN.x)*(locN.y)*(locN.z)*sizeof(double), cudaMemcpyHostToDevice );
	checkErrors("copy double data to device", __FILE__, __LINE__);
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, localN locN, consts def)
{
	cudaMemcpy( DevArrayPtr, HostArrayPtr, (locN.x)*(locN.y)*(locN.z)*sizeof(int), cudaMemcpyHostToDevice );
	checkErrors("copy int data to device", __FILE__, __LINE__);
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, localN locN, consts def)
{
	cudaMalloc((void**) DevBuffer,  2 * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).P_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).P_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).S_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ro_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ro_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ux_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uy_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uz_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ux_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uy_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uz_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).Xi_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).Xi_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_w),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_n),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_w_old),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_n_old),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).media),  (locN.x) * (locN.y) * (locN.z) * sizeof(int));

#ifdef B_L
	cudaMalloc((void**) &((*ArraysPtr).K),  (locN.x) * (locN.y) * (locN.z) * sizeof(double));
#endif

	checkErrors("memory allocation", __FILE__, __LINE__);
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	cudaFree (DevBuffer);
	cudaFree (DevArraysPtr.P_w);
	cudaFree (DevArraysPtr.P_n);
	cudaFree (DevArraysPtr.S_n);
	cudaFree (DevArraysPtr.ro_w);
	cudaFree (DevArraysPtr.ro_n);
	cudaFree (DevArraysPtr.ux_w);
	cudaFree (DevArraysPtr.uy_w);
	cudaFree (DevArraysPtr.uz_w);
	cudaFree (DevArraysPtr.ux_n);
	cudaFree (DevArraysPtr.uy_n);
	cudaFree (DevArraysPtr.uz_n);
	cudaFree (DevArraysPtr.Xi_w);
	cudaFree (DevArraysPtr.Xi_n);
	cudaFree (DevArraysPtr.roS_w);
	cudaFree (DevArraysPtr.roS_n);
	cudaFree (DevArraysPtr.roS_w_old);
	cudaFree (DevArraysPtr.roS_n_old);
	cudaFree (DevArraysPtr.media);

#ifdef B_L
	cudaFree (DevArraysPtr.K);
#endif

	checkErrors("memory reliase", __FILE__, __LINE__);
}

// Инициализация ускорителя
// Расчет происходит на ускорителе, номер которого равен
// номеру запускающего процессора
void device_initialization(int rank, int* blocksX, int* blocksY, int* blocksZ, localN locN, consts def)
{
	// Было бы очень неплохо вместо GPU_PER_NODE использовать cudaGetDeviceCount
	//int deviceCount;
	//cudaGetDeviceCount ( &deviceCount );

	// Считаем, что ядер на узле не меньше, чем ускорителей
	int device=rank % GPU_PER_NODE;
	cudaSetDevice(device);

	// Количество запускаемых блоков
	// Если число точек сетки не кратно размеру блока,
	// то количество блоков будет на 1 больше.
	*blocksX=(locN.x)/BlockNX;
	if ((locN.x)%BlockNX!=0)
		(*blocksX)++;
	*blocksY=(locN.y)/BlockNY;
	if ((locN.y)%BlockNY!=0)
		(*blocksY)++;
	*blocksZ=(locN.z)/BlockNZ;
	if ((locN.z)%BlockNZ!=0)
		(*blocksZ)++;

	consts* deff=new consts[1];
	deff[0]=def;
	cudaMemcpyToSymbol ( gpu_def, deff, sizeof ( consts ));
	checkErrors("constant memory copy", __FILE__, __LINE__);

	cudaDeviceProp devProp;
    cudaGetDeviceProperties ( &devProp, device );
        
	if (devProp.major < 2)
		printf ("\nError! Compute capability < 2, rank=%d\n",rank);
        
	if (!rank)
	{
		//printf ( "Device %d\n", device );
		printf ( "Name : %s\n", devProp.name );
		printf ( "Compute capability : %d.%d\n", devProp.major, devProp.minor );
		printf ( "Total Global Memory : %ld\n", devProp.totalGlobalMem );
		printf ( "Shared memory per block: %d\n", devProp.sharedMemPerBlock );
		printf ( "Registers per block : %d\n", devProp.regsPerBlock );
		printf ( "Warp size : %d\n", devProp.warpSize );
		printf ( "Max threads per block : %d\n", devProp.maxThreadsPerBlock );
		printf ( "Total constant memory : %d\n", devProp.totalConstMem );
		printf ( "Number of multiprocessors: %d\n",  devProp.multiProcessorCount);
		printf ( "Kernel execution timeout: %s\n\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
		for (int i = 0; i < 3; ++i)
			printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
		for (int i = 0; i < 3; ++i)
			printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);


		// Максимальный размер расчетной сетки для ускорителя
		// sizeof(ptr_Arrays)/4 - количество параметров в точке, т.к. 4 -размер одного указателя
		printf ( "\nTotal NAPL_Filtration grid size : %d\n\n", devProp.totalGlobalMem/(sizeof(ptr_Arrays)*sizeof(double)/4) );
	}

		// (locN.x)+2 потому что 2NyNz на буфер обмена выделяется
	if ( (locN.x+2)*(locN.y)*(locN.z) > (devProp.totalGlobalMem/(sizeof(ptr_Arrays)*sizeof(double)/4)))
		printf ("\nError! Not enough memory at GPU, rank=%d\n",rank);
	fflush( stdout);

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

__global__ void load_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer, localN locN)
{
	int j=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.y+blockIdx.y*blockDim.y;

	if (j<locN.y && k<(locN.z))
	{
		DevBuffer[j+(locN.y)*k]=DevArrayPtr[1+(locN.x)*j+(locN.x)*(locN.y)*k];
		DevBuffer[j+(locN.y)*k+(locN.y)*(locN.z)]=DevArrayPtr[(locN.x)-2+(locN.x)*j+(locN.x)*(locN.y)*k];

		device_test_nan(DevBuffer[j+(locN.y)*k], __FILE__, __LINE__);
		device_test_nan(DevBuffer[j+(locN.y)*k+(locN.y)*(locN.z)], __FILE__, __LINE__);
	}
}


void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	load_exchange_data_kernel<<<dim3(blocksY,blocksZ), dim3(BlockNY,BlockNZ)>>>(DevArrayPtr, DevBuffer, locN); 
	checkErrors("load_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy( HostBuffer, DevBuffer, 2*(def.Ny)*(def.Nz)*sizeof(double), cudaMemcpyDeviceToHost );
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void save_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer, localN locN, int rank, parts_sizes parts)
{
	int j=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.y+blockIdx.y*blockDim.y;

	if (j<(locN.y) && k<(locN.z))
	{
		if (rank!=0)
		{
			DevArrayPtr[(locN.x)*j+(locN.x)*(locN.y)*k]=DevBuffer[j+(locN.y)*k];
			device_test_nan(DevArrayPtr[(locN.x)*j+(locN.x)*(locN.y)*k], __FILE__, __LINE__);
		}
		if (rank!=parts.x-1)
		{
			DevArrayPtr[(locN.x)-1+(locN.x)*j+(locN.x)*(locN.y)*k]=DevBuffer[j+(locN.y)*k+(locN.y)*(locN.z)];
			device_test_nan(DevArrayPtr[(locN.x)-1+(locN.x)*j+(locN.x)*(locN.y)*k], __FILE__, __LINE__);
		}
	}
}

void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	cudaMemcpy( DevBuffer, HostBuffer, 2*(locN.y)*(locN.z)*sizeof(double), cudaMemcpyHostToDevice );
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_kernel<<<dim3(blocksY,blocksZ), dim3(BlockNY,BlockNZ)>>>(DevArrayPtr, DevBuffer, locN, rank, parts); 
	checkErrors("save_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}