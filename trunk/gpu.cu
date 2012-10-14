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

// �������� �� ����� �������������� ���������
__device__ int device_is_injection_well(int i, int j, int k)
{
#ifdef B_L
	if (((i == 1) && (j == 1)) || ((i == 0) && (j == 0)) || ((i == 1) && (j == 0)) || ((i == 0) && (j == 1)))
#endif
#ifdef THREE_PHASE
	if ((j == 1) && (i > 0) && (i < (gpu_def->locNx) / 3 - 1) && (((gpu_def->locNz) < 2) || (k > 0) && (k < (gpu_def->locNz) / 3 - 1)))
#endif
#ifndef TWO_PHASE
			return 1;
		else
#endif
			return 0;
}

// �������� �� ����� ���������� ���������
__device__ int device_is_output_well(int i, int j, int k)
{
#ifdef B_L
	if (((i == gpu_def->Nx - 2) && (j == gpu_def->Ny - 2)) || ((i == gpu_def->Nx - 1) && (j == gpu_def->Ny - 1)) || ((i == gpu_def->Nx - 1) && (j == gpu_def->Ny - 2)) || ((i == gpu_def->Nx - 2) && (j == gpu_def->Ny - 1)))
#endif
#ifdef THREE_PHASE
	if ((j == 1) && (i > 0) && (i < (gpu_def->locNx) / 3 - 1) && (((gpu_def->locNz) < 2) || (k > 0) && (k < (gpu_def->locNz) / 3 - 1)))
#endif
#ifndef TWO_PHASE
			return 1;
		else
#endif
			return 0;
}

// ������������� �������� ���������/���������� ��������� q_i �� ���������
__device__ void device_wells_q(ptr_Arrays DevArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g)
{
#ifdef B_L
	// �������������� ��������
	if (device_is_injection_well(i, j, k))
	{
		*q_w = gpu_def->Q;
		*q_n = 0.;
		*q_g = 0.;
	}

	// ���������� ��������
	if (device_is_output_well(i, j, k))
	{
		*q_g = 0;
		double k_w=0., k_n=0.;
		device_assing_k(&k_w, &k_n, 1. - DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);

		double F_bl = (k_w / gpu_def->mu_w) / (k_w / gpu_def->mu_w + k_n / gpu_def->mu_n);
		*q_w = -1. * gpu_def->Q * F_bl;
		*q_n = -1. * gpu_def->Q * (1. - F_bl);
	}
#endif

#ifdef THREE_PHASE

	double q = 0.;
	*q_w = 0.0;
	*q_g = 0.0;
	*q_n = 0.0;

/*	if (device_is_injection_well(i, j, k))
	{
		*q_w = 0.01;
		*q_g = 0.005;
		*q_n = 0.02;
	}
	if (device_is_output_well(i, j, k))
	{
		q = 0.035;

		*q_w = -q * DevArraysPtr.S_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)];
		*q_g = -q * (1 - DevArraysPtr.S_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);
		*q_n = -q * DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)];
	}*/
#endif
}

__global__ void assign_ro_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		DevArraysPtr.ro_w[local] = (gpu_def->ro0_w) * (1. + (gpu_def->beta_w) * (DevArraysPtr.P_w[local] - (gpu_def->P_atm)));
		DevArraysPtr.ro_n[local] = (gpu_def->ro0_n) * (1. + (gpu_def->beta_n) * (DevArraysPtr.P_n[local] - (gpu_def->P_atm)));
	#ifdef THREE_PHASE
		DevArraysPtr.ro_g[local] = (gpu_def->ro0_g) * DevArraysPtr.P_g[local] / (gpu_def->P_atm);
		device_test_positive(DevArraysPtr.ro_g[local], __FILE__, __LINE__);
	#endif
		device_test_positive(DevArraysPtr.ro_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_n[local], __FILE__, __LINE__);
	}
}

// ���������� ���������� �����, ����� ������� ����� ����������� �������� �� ������� (i1, j1, k1)
__device__ int device_set_boundary_basic_coordinate(int i, int j, int k)
{
	int i1, j1, k1;

	i1 = i;
	j1 = j;
	k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (gpu_def->locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (gpu_def->locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((gpu_def->locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (gpu_def->locNz) - 1) && ((gpu_def->locNz) > 2))
	{
		k1 --;
	}

	return (i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy));
}

// ������ ����������� ��������
__device__ double central_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *(ptr-1) )/ (2. * (gpu_def->hx));	
		}
	case 'y':
		{
			return (*(ptr + gpu_def->locNx) - *(ptr - gpu_def->locNx) )/ (2. * (gpu_def->hy));
		}
	case 'z':
		{
			return (*(ptr + gpu_def->locNx * (gpu_def->locNy)) - *(ptr - gpu_def->locNx * (gpu_def->locNy)) )/ (2. * (gpu_def->hz));
		}
	default:
		{
			device_print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// ������ ������������ ��������
__device__ double directed_difference (double x1, double x2, double* Xi, double* ro, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-1)) * (*(ro-1)) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+1)) * (*(ro+1))) / gpu_def->hx;
		}
	case 'y':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi - gpu_def->locNx)) * (*(ro - gpu_def->locNx)) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi + gpu_def->locNx)) * (*(ro + gpu_def->locNx))) / gpu_def->hy;
		}
	case 'z':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
				(x1 + fabs(x1)) / 2. * (-1.) * (*(Xi - gpu_def->locNx * (gpu_def->locNy))) * (*(ro - gpu_def->locNx * (gpu_def->locNy))) +
				(x2 - fabs(x2)) / 2. * (-1.) * (*(Xi + gpu_def->locNx * (gpu_def->locNy))) * (*(ro + gpu_def->locNx * (gpu_def->locNy)))) / gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// ������ ����� ��������
__device__ double left_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*ptr - *(ptr-1) )/ gpu_def->hx;	
		}
	case 'y':
		{
			return (*ptr - *(ptr-gpu_def->locNx) )/ gpu_def->hy;
		}
	case 'z':
		{
			return (*ptr - *(ptr - gpu_def->locNx * (gpu_def->locNy)) )/ gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [left_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// ������ ������ ��������
__device__ double right_difference (double* ptr, char axis)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *ptr )/ gpu_def->hx;	
		}
	case 'y':
		{
			return (*(ptr + gpu_def->locNx) - *ptr )/ gpu_def->hy;
		}
	case 'z':
		{
			return (*(ptr + gpu_def->locNx * (gpu_def->locNy)) - *ptr )/ gpu_def->hz;
		}
	default:
		{
			device_print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// �������������� ��������� ��������� ���������� � ����������
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������
// (���� 2 ������ � ����� ������,�� +2 �����).
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� rank==0)
__device__ int device_local_to_global(int local_index, char axis)
{
	int global_index = local_index;
	switch (axis)
	{
		case 'x':
		{
			global_index += gpu_def->rankx * (gpu_def->Nx) / gpu_def->sizex + min(gpu_def->rankx, gpu_def->Nx % gpu_def->sizex);
			break;
		}
		case 'y':
		{
			global_index += gpu_def->ranky * (gpu_def->Ny) / gpu_def->sizey + min(gpu_def->ranky, gpu_def->Ny % gpu_def->sizey);
			break;
		}
		case 'z':
		{
			global_index += gpu_def->rankz * (gpu_def->Nz) / gpu_def->sizez + min(gpu_def->rankz, gpu_def->Nz % gpu_def->sizez);
			break;
		}
		default:
		{
			//CUPRINTF("Error!");
		}
	}
	//some_test(global_index);
	return global_index;
}

// �������� �� ����� �������� (�.�. �� ��������������� ������ ��� ������ �� ��������)
__device__ int device_is_active_point(int i, int j, int k)
{
	if ((gpu_def->rank % gpu_def->sizex != 0 && i == 0) || (gpu_def->rank % gpu_def->sizex != gpu_def->sizex - 1 && i == gpu_def->locNx - 1)
	    || ((gpu_def->rank % ((gpu_def->sizex) * (gpu_def->sizey))) / gpu_def->sizex != 0 && j == 0)	|| ((gpu_def->rank % ((gpu_def->sizex) * (gpu_def->sizey))) / gpu_def->sizex != gpu_def->sizey - 1 && j == gpu_def->locNy - 1)
	    || (((gpu_def->rank / (gpu_def->sizex) / (gpu_def->sizey) != 0 && k == 0) || (gpu_def->rank / (gpu_def->sizex) / (gpu_def->sizey) == gpu_def->sizez - 1 && k == gpu_def->locNz - 1)) && gpu_def->sizez > 1))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

// ������� ���������� "�����������" ���������
__device__ double device_ro_eff_gdy(ptr_Arrays DevArraysPtr, int local)
{
#ifdef THREE_PHASE
	double ro_g_dy = (DevArraysPtr.ro_g[local] * (1. - DevArraysPtr.S_w[local] - DevArraysPtr.S_n[local]) + DevArraysPtr.ro_w[local] * DevArraysPtr.S_w[local]
	                  + DevArraysPtr.ro_n[local] * DevArraysPtr.S_n[local]) * (DevArraysPtr.m[local]) * (gpu_def->g_const) * (gpu_def->hy);
#else
	double ro_g_dy = (DevArraysPtr.ro_n[local] * DevArraysPtr.S_n[local] + DevArraysPtr.ro_w[local] * (1 - DevArraysPtr.S_n[local]))
					  * (DevArraysPtr.m[local]) * (gpu_def->g_const) * (gpu_def->hy);
#endif
	return ro_g_dy;
}

// ������ ����������, �������� NAPL P2 � Xi �� ���� ������ �����
void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_P_Xi_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign P, Xi", __FILE__, __LINE__);
	assign_ro_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign ro", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ �������� ���� P1 � ������������ NAPL S2 �� ���� ������ �����
void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	Newton_method_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign Pw and Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ �������� � ������ ����� �����
__global__ void assign_u_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int local=i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		if ((gpu_def->Nx) > 2)
		{
			if (i == 0)
			{
				DevArraysPtr.ux_w[local] = DevArraysPtr.Xi_w[local] * right_difference(DevArraysPtr.P_w+local, 'x');
				DevArraysPtr.ux_n[local] = DevArraysPtr.Xi_n[local] * right_difference(DevArraysPtr.P_n+local, 'x');
#ifdef THREE_PHASE
				DevArraysPtr.ux_g[local] = DevArraysPtr.Xi_g[local] * right_difference(DevArraysPtr.P_g+local, 'x');
#endif
			}
			else
			{
				if (i == (gpu_def->locNx) - 1)
				{
					DevArraysPtr.ux_w[local] = DevArraysPtr.Xi_w[local] * left_difference(DevArraysPtr.P_w+local, 'x');
					DevArraysPtr.ux_n[local] = DevArraysPtr.Xi_n[local] * left_difference(DevArraysPtr.P_n+local, 'x');
#ifdef THREE_PHASE
					DevArraysPtr.ux_g[local] = DevArraysPtr.Xi_g[local] * left_difference(DevArraysPtr.P_g+local, 'x');
#endif
				}
				else
				{
					DevArraysPtr.ux_w[local] = DevArraysPtr.Xi_w[local] * central_difference (DevArraysPtr.P_w+local, 'x');
					DevArraysPtr.ux_n[local] = DevArraysPtr.Xi_n[local] * central_difference (DevArraysPtr.P_n+local, 'x');
#ifdef THREE_PHASE
					DevArraysPtr.ux_g[local] = DevArraysPtr.Xi_g[local] * central_difference (DevArraysPtr.P_g+local, 'x');
#endif
				}
			}
		}
		else
		{
			DevArraysPtr.ux_w[local] = 0.;
			DevArraysPtr.ux_n[local] = 0.;
#ifdef THREE_PHASE
			DevArraysPtr.ux_g[local] = 0.;
#endif
		}

		if ((gpu_def->Ny) > 2)
		{
			if (j == 0)
			{
				DevArraysPtr.uy_w[local] = DevArraysPtr.Xi_w[local] * (right_difference (DevArraysPtr.P_w+local, 'y') - DevArraysPtr.ro_w[local] * (gpu_def->g_const));
				DevArraysPtr.uy_n[local] = DevArraysPtr.Xi_n[local] * (right_difference (DevArraysPtr.P_n+local, 'y') - DevArraysPtr.ro_n[local] * (gpu_def->g_const));
#ifdef THREE_PHASE
				DevArraysPtr.uy_g[local] = DevArraysPtr.Xi_g[local] * (right_difference (DevArraysPtr.P_g+local, 'y') - DevArraysPtr.ro_g[local] * (gpu_def->g_const));
#endif
			}
			else
			{
				if (j == (gpu_def->locNy) - 1)
				{
					DevArraysPtr.uy_w[local] = DevArraysPtr.Xi_w[local] * (left_difference (DevArraysPtr.P_w+local, 'y') - DevArraysPtr.ro_w[local] * (gpu_def->g_const));
					DevArraysPtr.uy_n[local] = DevArraysPtr.Xi_n[local] * (left_difference (DevArraysPtr.P_n+local, 'y') - DevArraysPtr.ro_n[local] * (gpu_def->g_const));
#ifdef THREE_PHASE
					DevArraysPtr.uy_g[local] = DevArraysPtr.Xi_g[local] * (left_difference (DevArraysPtr.P_g+local, 'y') - DevArraysPtr.ro_g[local] * (gpu_def->g_const));
#endif
				}
				else
				{
					DevArraysPtr.uy_w[local] = DevArraysPtr.Xi_w[local] * (central_difference (DevArraysPtr.P_w+local, 'y')	- DevArraysPtr.ro_w[local] * (gpu_def->g_const));
					DevArraysPtr.uy_n[local] = DevArraysPtr.Xi_n[local] * (central_difference (DevArraysPtr.P_n+local, 'y')	- DevArraysPtr.ro_n[local] * (gpu_def->g_const));
#ifdef THREE_PHASE
					DevArraysPtr.uy_g[local] = DevArraysPtr.Xi_g[local] * (central_difference (DevArraysPtr.P_g+local, 'y')	- DevArraysPtr.ro_g[local] * (gpu_def->g_const));
#endif
				}
			}
		}
		else
		{
			DevArraysPtr.uy_w[local] = 0.;
			DevArraysPtr.uy_n[local] = 0.;
#ifdef THREE_PHASE
			DevArraysPtr.uy_g[local] = 0.;
#endif
		}

		if ((gpu_def->Nz) > 2)
		{
			if (k == 0)
			{
				DevArraysPtr.uz_w[local] = DevArraysPtr.Xi_w[local] * right_difference (DevArraysPtr.P_w+local, 'z');
				DevArraysPtr.uz_n[local] = DevArraysPtr.Xi_n[local] * right_difference (DevArraysPtr.P_n+local, 'z');
#ifdef THREE_PHASE
				DevArraysPtr.uz_g[local] = DevArraysPtr.Xi_g[local] * right_difference (DevArraysPtr.P_g+local, 'z');
#endif
			}
			else
			{
				if (k == (gpu_def->locNz) - 1)
				{
					DevArraysPtr.uz_w[local] = DevArraysPtr.Xi_w[local] * left_difference (DevArraysPtr.P_w+local, 'z');
					DevArraysPtr.uz_n[local] = DevArraysPtr.Xi_n[local] * left_difference (DevArraysPtr.P_n+local, 'z');
#ifdef THREE_PHASE
					DevArraysPtr.uz_g[local] = DevArraysPtr.Xi_g[local] * left_difference (DevArraysPtr.P_g+local, 'z');
#endif
				}
				else
				{
					DevArraysPtr.uz_w[local] = DevArraysPtr.Xi_w[local] * central_difference (DevArraysPtr.P_w+local, 'z');
					DevArraysPtr.uz_n[local] = DevArraysPtr.Xi_n[local] * central_difference (DevArraysPtr.P_n+local, 'z');
#ifdef THREE_PHASE
					DevArraysPtr.uz_g[local] = DevArraysPtr.Xi_g[local] * central_difference (DevArraysPtr.P_g+local, 'z');
#endif
				}
			}
		}
		else
		{
			DevArraysPtr.uz_w[local] = 0.;
			DevArraysPtr.uz_n[local] = 0.;
#ifdef THREE_PHASE
			DevArraysPtr.uz_g[local] = 0.;
#endif
		}

		device_test_u(DevArraysPtr.ux_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.ux_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uy_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uy_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uz_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uz_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
#ifdef THREE_PHASE
		device_test_u(DevArraysPtr.ux_g[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uy_g[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_u(DevArraysPtr.uz_g[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
#endif
	}
}

// ������ ��������� �� ���� ������ �����
void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_u_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign u", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ ro*S � ������ ����� ����� ������� ������������ ���������
__global__ void assign_roS_kernel_nr(ptr_Arrays DevArraysPtr, double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx) - 1) && (j < gpu_def->locNy - 1) && (k < (gpu_def->locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2)))
	{
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		if(! DevArraysPtr.m[local])
			return;

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// �������� q �� ���������
		device_wells_q(DevArraysPtr, i, j, k, &q_w, &q_n, &q_g);

#ifdef THREE_PHASE
		DevArraysPtr.roS_w[local] = DevArraysPtr.ro_w[local] * DevArraysPtr.S_w[local];
		DevArraysPtr.roS_g[local] = DevArraysPtr.ro_g[local]
		* (1. - DevArraysPtr.S_w[local] - DevArraysPtr.S_n[local]);

		double fx_g, fy_g, fz_g, A3 = 0.;
#else
		DevArraysPtr.roS_w[local] = DevArraysPtr.ro_w[local] * (1. - DevArraysPtr.S_n[local]);
#endif
		DevArraysPtr.roS_n[local] = DevArraysPtr.ro_n[local] * DevArraysPtr.S_n[local];
		//double Pw = DevArraysPtr.P_w[local];
		//double Pn = DevArraysPtr.P_n[local];

		double x1, x2, y1, y2, z1, z2, fx_w, fy_w, fz_w, fx_n, fy_n, fz_n, A1 = 0., A2 = 0.;

		if ((gpu_def->Nz) < 2)
		{
			fz_w = 0.;
			fz_n = 0.;
#ifdef THREE_PHASE
			fz_g = 0.;
#endif
		}
		else
		{
			z2 = -1. * right_difference (DevArraysPtr.P_w+local, 'z');
			z1 = -1. * left_difference (DevArraysPtr.P_w+local, 'z');
			fz_w = directed_difference (z1, z2, DevArraysPtr.Xi_w+local, DevArraysPtr.ro_w+local, 'z');

			z2 = -1. * right_difference (DevArraysPtr.P_n+local, 'z'); 
			z1 = -1. * left_difference (DevArraysPtr.P_n+local, 'z'); 
			fz_n = directed_difference (z1, z2, DevArraysPtr.Xi_n+local, DevArraysPtr.ro_n+local, 'z');
#ifdef THREE_PHASE
			z2 = -1. * right_difference (DevArraysPtr.P_g+local, 'z'); 
			z1 = -1. * left_difference (DevArraysPtr.P_g+local, 'z'); 
			fz_g = directed_difference (z1, z2, DevArraysPtr.Xi_g+local, DevArraysPtr.ro_g+local, 'z');
#endif
		}

		x2 = -1. * right_difference (DevArraysPtr.P_w+local, 'x'); //-(DevArraysPtr.P_w[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - Pw) / gpu_def->hx;
		x1 = -1. * left_difference (DevArraysPtr.P_w+local, 'x'); //-(Pw - DevArraysPtr.P_w[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / gpu_def->hx;
		y2 = -1. * right_difference (DevArraysPtr.P_w+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_w[local]);
		y1 = -1. * left_difference (DevArraysPtr.P_w+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_w[local]);
		fx_w = directed_difference (x1, x2, DevArraysPtr.Xi_w+local, DevArraysPtr.ro_w+local, 'x');
		fy_w = directed_difference (y1, y2, DevArraysPtr.Xi_w+local, DevArraysPtr.ro_w+local, 'y');

		x2 = -1. * right_difference (DevArraysPtr.P_n+local, 'x'); 
		x1 = -1. * left_difference (DevArraysPtr.P_n+local, 'x'); 
		y2 = -1. * right_difference (DevArraysPtr.P_n+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_n[local]);
		y1 = -1. * left_difference (DevArraysPtr.P_n+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_n[local]);
		fx_n = directed_difference (x1, x2, DevArraysPtr.Xi_n+local, DevArraysPtr.ro_n+local, 'x');
		fy_n = directed_difference (y1, y2, DevArraysPtr.Xi_n+local, DevArraysPtr.ro_n+local, 'y');

		A1 = DevArraysPtr.roS_w[local] - (gpu_def->dt / DevArraysPtr.m[local]) * (-q_w + fx_w + fy_w + fz_w);
		A2 = DevArraysPtr.roS_n[local] - (gpu_def->dt / DevArraysPtr.m[local]) * (-q_n + fx_n + fy_n + fz_n);

		DevArraysPtr.roS_w_old[local] = DevArraysPtr.roS_w[local];
		DevArraysPtr.roS_n_old[local] = DevArraysPtr.roS_n[local];
		DevArraysPtr.roS_w[local] = A1;
		DevArraysPtr.roS_n[local] = A2;

		device_test_positive(DevArraysPtr.roS_w[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[local], __FILE__, __LINE__);

#ifdef THREE_PHASE
		x2 = -1. * right_difference (DevArraysPtr.P_g+local, 'x'); 
		x1 = -1. * left_difference (DevArraysPtr.P_g+local, 'x'); 
		y2 = -1. * right_difference (DevArraysPtr.P_g+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_g[local]);
		y1 = -1. * left_difference (DevArraysPtr.P_g+local, 'y') + gpu_def->g_const * (DevArraysPtr.ro_g[local]);
		fx_g = directed_difference (x1, x2, DevArraysPtr.Xi_g+local, DevArraysPtr.ro_g+local, 'x');
		fy_g = directed_difference (y1, y2, DevArraysPtr.Xi_g+local, DevArraysPtr.ro_g+local, 'y');

		A3 = DevArraysPtr.roS_g[local] - (gpu_def->dt / DevArraysPtr.m[local]) * (q_g + fx_g + fy_g + fz_g);

		DevArraysPtr.roS_g_old[local] = DevArraysPtr.roS_g[local];
		DevArraysPtr.roS_g[local] = A3;

		device_test_positive(DevArraysPtr.roS_g[local], __FILE__, __LINE__);
#endif
	}
}

// ������ ro*S � ������ ����� �����
__global__ void assign_roS_kernel(ptr_Arrays DevArraysPtr, double t)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx) - 1) && (j < gpu_def->locNy - 1) && (k < (gpu_def->locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != (gpu_def->locNz))) || ((gpu_def->locNz) < 2)))
	{
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		double S_n = DevArraysPtr.S_n[local];
		double roS_w = DevArraysPtr.ro_w[local] * (1 - S_n);
		double roS_n = DevArraysPtr.ro_n[local] * S_n;

		double divgrad1, divgrad2, Tx1, Ty1, Tx2, Ty2, Tz1, Tz2, A1 = 0, A2 = 0;

		if ((gpu_def->Nz) < 2)
		{
			divgrad1 = 0;
			divgrad2 = 0;
			Tz1 = 0;
			Tz2 = 0;
		}
		else
		{
			divgrad1 = (DevArraysPtr.m[local] * (gpu_def->l) * (gpu_def->c_w) / 2.) * (DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] * (1. - DevArraysPtr.S_n[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)]) - 2 * DevArraysPtr.ro_w[local] * (1. - S_n) + DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)] * (1. - DevArraysPtr.S_n[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)])) / ((gpu_def->hz) * (gpu_def->hz));
			divgrad2 = (DevArraysPtr.m[local] * (gpu_def->l) * (gpu_def->c_n) / 2.) * (DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.S_n[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] - 2 * DevArraysPtr.ro_n[local] * S_n + DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)] * (DevArraysPtr.S_n[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)])) / ((gpu_def->hz) * (gpu_def->hz));
			Tz1 = (DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uz_w[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uz_w[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)]) / (2. * (gpu_def->hz));
			Tz2 = (DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uz_n[i + j * (gpu_def->locNx) + (k + 1) * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uz_n[i + j * (gpu_def->locNx) + (k - 1) * (gpu_def->locNx) * (gpu_def->locNy)]) / (2. * (gpu_def->hz));
		}

		divgrad1 += (DevArraysPtr.m[local] * (gpu_def->l) * (gpu_def->c_w) / 2.) *
		            ((DevArraysPtr.ro_w[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - DevArraysPtr.S_n[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) - 2 * DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - DevArraysPtr.S_n[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)])) / ((gpu_def->hx) * (gpu_def->hx)) +
		             (DevArraysPtr.ro_w[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - DevArraysPtr.S_n[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) - 2 * DevArraysPtr.ro_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * (1 - DevArraysPtr.S_n[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)])) / ((gpu_def->hy) * (gpu_def->hy)));

		divgrad2 += (DevArraysPtr.m[local] * (gpu_def->l) * (gpu_def->c_n) / 2.) *
		            ((DevArraysPtr.ro_n[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.S_n[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - 2 * DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * S_n + DevArraysPtr.ro_n[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.S_n[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / ((gpu_def->hx) * (gpu_def->hx)) +
		             (DevArraysPtr.ro_n[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.S_n[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - 2 * DevArraysPtr.ro_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * S_n + DevArraysPtr.ro_n[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.S_n[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / ((gpu_def->hy) * (gpu_def->hy)));

		Tx1 = (DevArraysPtr.ro_w[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.ux_w[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_w[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.ux_w[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / (2 * (gpu_def->hx));
		Ty1 = (DevArraysPtr.ro_w[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_w[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_w[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_w[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / (2 * (gpu_def->hy));
		Tx2 = (DevArraysPtr.ro_n[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.ux_n[i + 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_n[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.ux_n[i - 1 + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / (2 * (gpu_def->hx));
		Ty2 = (DevArraysPtr.ro_n[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_n[i + (j + 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - DevArraysPtr.ro_n[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.uy_n[i + (j - 1) * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]) / (2 * (gpu_def->hy));

		double q_w = 0;
		double q_n = 0;

#ifdef B_L
		// � ������ ���������� ��������� ����������� ��������
		if ((i == gpu_def->Nx / 2) && (j == gpu_def->Ny - 3) && (k == gpu_def->Nz / 2))
		{
			q_w = gpu_def->Q;
			q_n = 0;
		}

		// � ������ ���������� ��������� ���������� ��������
		if ((i == gpu_def->Nx - 3) && (j == 3) && (k == gpu_def->Nz - 3))
		{
			double k_w=0., k_n=0.;
			device_assing_k(&k_w, &k_n, 1. - DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)]);
			double F_bl = 0;
			F_bl = (k_w / gpu_def->mu_w) / (k_w / gpu_def->mu_w + k_n / gpu_def->mu_n);
			q_w = -1 * (gpu_def->Q) * F_bl;
			q_n = -1 * (gpu_def->Q) * (1 - F_bl);
		}
#endif

		if ((t < 2 * (gpu_def->dt)) || TWO_LAYERS)
		{
			A1 = roS_w + ((gpu_def->dt) / DevArraysPtr.m[local]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = roS_n + ((gpu_def->dt) / DevArraysPtr.m[local]) * (q_w + divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2 * (gpu_def->dt) * (gpu_def->dt)) / (DevArraysPtr.m[local] * ((gpu_def->dt) + 2 * (gpu_def->tau))) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1 + (2 * roS_w * (DevArraysPtr.m[local]) * (gpu_def->tau)) / ((gpu_def->dt) * (gpu_def->dt)) + DevArraysPtr.roS_w_old[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.m[local] * ((gpu_def->dt) - 2 * (gpu_def->tau)) / (2 * (gpu_def->dt) * (gpu_def->dt)));
			A2 = (2 * (gpu_def->dt) * (gpu_def->dt)) / (DevArraysPtr.m[local] * ((gpu_def->dt) + 2 * (gpu_def->tau))) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2 + (2 * roS_n * (DevArraysPtr.m[local]) * (gpu_def->tau)) / ((gpu_def->dt) * (gpu_def->dt)) + DevArraysPtr.roS_n_old[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] * DevArraysPtr.m[local] * ((gpu_def->dt) - 2 * (gpu_def->tau)) / (2 * (gpu_def->dt) * (gpu_def->dt)));
		}

		DevArraysPtr.roS_w_old[local] = roS_w;
		DevArraysPtr.roS_n_old[local] = roS_n;
		DevArraysPtr.roS_w[local] = A1;
		DevArraysPtr.roS_n[local] = A2;

		device_test_positive(DevArraysPtr.roS_w_old[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
	}
}

// ������ ro*S �� ���� ������ �����
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

// ���������� ��������� �������
void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	Border_S_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign S", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	Border_P_kernel <<< dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX, BlockNY, BlockNZ)>>>(DevArraysPtr);
	checkErrors("assign Pw", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������� �������� ������ � ������ �����
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy(HostArrayPtr, DevArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
}

// ������� �������� ������ ���� double � ������ ����������
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy double data to device", __FILE__, __LINE__);
}

// ������� �������� ������ ���� int � ������ ����������
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, consts def)
{
	cudaMemcpy(DevArrayPtr, HostArrayPtr, (def.locNx) * (def.locNy) * (def.locNz)*sizeof(int), cudaMemcpyHostToDevice);
	checkErrors("copy int data to device", __FILE__, __LINE__);
}

// ��������� ������ ���������� ��� ������ ����� ��������� �������
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, consts def)
{
	int buffer_size = 0;

	if(def.sizex > 1)
		buffer_size = (def.locNy) * (def.locNz);
	if(def.sizey > 1 && (def.locNx) * (def.locNz) > buffer_size)
		buffer_size = (def.locNx) * (def.locNz);
	if(def.sizez > 1 && (def.locNx) * (def.locNy) > buffer_size)
		buffer_size = (def.locNx) * (def.locNy);

	if(buffer_size)
		cudaMalloc((void**) DevBuffer,  buffer_size * sizeof(double));
	
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
	cudaMalloc((void**) & ((*ArraysPtr).m), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
	cudaMalloc((void**) & ((*ArraysPtr).K), (def.locNx) * (def.locNy) * (def.locNz) * sizeof(double));
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

	checkErrors("memory allocation", __FILE__, __LINE__);
}

// ����������� ������ ���������� �� ��� ������� ����� ��������� �������
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
	cudaFree(DevArraysPtr.m);
	cudaFree(DevArraysPtr.K);
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

	checkErrors("memory release", __FILE__, __LINE__);
}

// ������������� ����������
// ������ ���������� �� ����������, ����� �������� �����
// ������ ������������ ����������
void device_initialization(consts* def)
{
	// ���� �� ����� ������� ������ GPU_PER_NODE ������������ cudaGetDeviceCount
	//int deviceCount;
	//cudaGetDeviceCount ( &deviceCount );

	// �������, ��� ���� �� ���� �� ������, ��� �����������
	int device = (*def).rank % GPU_PER_NODE;
	cudaSetDevice(device);

	// ���������� ����������� ������
	// ���� ����� ����� ����� �� ������ ������� �����,
	// �� ���������� ������ ����� �� 1 ������.
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
		//printf("Kernel execution timeout: %s\n\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
		}
		for (int i = 0; i < 3; ++i)
		{
			printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
		}


		// ������������ ������ ��������� ����� ��� ����������
		// sizeof(ptr_Arrays)/4 - ���������� ���������� � �����, �.�. 4 -������ ������ ���������
		printf("\nTotal NAPL_Filtration grid size : %d\n\n", devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4));
	}

	// (def.locNx)+2 ������ ��� 2NyNz �� ����� ������ ����������
	if (((*def).locNx + 2) * ((*def).locNy) * ((*def).locNz) > (devProp.totalGlobalMem / (sizeof(ptr_Arrays)*sizeof(double) / 4)))
	{
		printf("\nError! Not enough memory at GPU, rank=%d\n", (*def).rank);
	}
	fflush(stdout);

	// �������������� ���������� cuPrintf ��� ������ ������ �� �������
	// ����� �� kernel
	cudaPrintfInit();
}

// ����������� ����������
void device_finalization(void)
{
	// ������������� ���������� cuPrintf ��� ������ ������ �� �������
	// ����� �� kernel
	cudaPrintfEnd();
}

__global__ void load_exchange_data_part_xl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevBuffer[j + (gpu_def->locNy)*k] = DevArrayPtr[1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[j + (gpu_def->locNy)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_xr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevBuffer[j + (gpu_def->locNy)*k] = DevArrayPtr[(gpu_def->locNx) - 2 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[j + (gpu_def->locNy)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_yl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevBuffer[i + (gpu_def->locNx)*k] = DevArrayPtr[i + (gpu_def->locNx) + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_yr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevBuffer[i + (gpu_def->locNx)*k] = DevArrayPtr[i + (gpu_def->locNx) * (gpu_def->locNy - 2) + (gpu_def->locNx) * (gpu_def->locNy) * k];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*k], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_zl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevBuffer[i + (gpu_def->locNx)*j] = DevArrayPtr[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy)];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*j], __FILE__, __LINE__);
	}
}

__global__ void load_exchange_data_part_zr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevBuffer[i + (gpu_def->locNx)*j] = DevArrayPtr[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 2)];
		device_test_nan(DevBuffer[i + (gpu_def->locNx)*j], __FILE__, __LINE__);
	}
}

void load_exchange_data_part_xl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_xl_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_xl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Ny) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void load_exchange_data_part_xr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_xr_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_xr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Ny) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void load_exchange_data_part_yl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_yl_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_yl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Nx) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void load_exchange_data_part_yr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_yr_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_yr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Nx) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void load_exchange_data_part_zl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_zl_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_zl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Nx) * (def.Ny) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void load_exchange_data_part_zr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_part_zr_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArrayPtr, DevBuffer);
	checkErrors("load_exchange_data_part_zr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy(HostBuffer, DevBuffer, (def.Nx) * (def.Ny) * sizeof(double), cudaMemcpyDeviceToHost);
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void save_exchange_data_part_xl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevArrayPtr[(gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[j + (gpu_def->locNy)*k];
		device_test_nan(DevArrayPtr[(gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_xr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (j < gpu_def->locNy && k < (gpu_def->locNz))
	{
		DevArrayPtr[(gpu_def->locNx) - 1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[j + (gpu_def->locNy)*k];
		device_test_nan(DevArrayPtr[(gpu_def->locNx) - 1 + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_yl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevArrayPtr[i + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[i + (gpu_def->locNx)*k];
		device_test_nan(DevArrayPtr[i + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_yr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && k < (gpu_def->locNz))
	{
		DevArrayPtr[i + (gpu_def->locNx) * (gpu_def->locNy - 1) + (gpu_def->locNx) * (gpu_def->locNy) * k] = DevBuffer[i + (gpu_def->locNx)*k];
		device_test_nan(DevArrayPtr[i + (gpu_def->locNx) * (gpu_def->locNy - 1) + (gpu_def->locNx) * (gpu_def->locNy) * k], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_zl_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevArrayPtr[i + (gpu_def->locNx) * j] = DevBuffer[i + (gpu_def->locNx)*j];
		device_test_nan(DevArrayPtr[i + (gpu_def->locNx) * j], __FILE__, __LINE__);
	}
}

__global__ void save_exchange_data_part_zr_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;

	if (i < gpu_def->locNx && j < (gpu_def->locNy))
	{
		DevArrayPtr[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 1)] = DevBuffer[i + (gpu_def->locNx)*j];
		device_test_nan(DevArrayPtr[i + (gpu_def->locNx) * j + (gpu_def->locNx) * (gpu_def->locNy) * (gpu_def->locNz - 1)], __FILE__, __LINE__);
	}
}

void save_exchange_data_part_xl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_xl_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_xl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void save_exchange_data_part_xr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNy) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_xr_kernel <<< dim3(def.blocksY, def.blocksZ), dim3(BlockNY, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_xr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void save_exchange_data_part_yl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNx) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_yl_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_yl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void save_exchange_data_part_yr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNx) * (def.locNz)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_yr_kernel <<< dim3(def.blocksX, def.blocksZ), dim3(BlockNX, BlockNZ)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_yr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void save_exchange_data_part_zl(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNx) * (def.locNy)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_zl_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_zl", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

void save_exchange_data_part_zr(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy(DevBuffer, HostBuffer, (def.locNx) * (def.locNy)*sizeof(double), cudaMemcpyHostToDevice);
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_part_zr_kernel <<< dim3(def.blocksX, def.blocksY), dim3(BlockNX, BlockNY)>>>(DevArrayPtr, DevBuffer);
	checkErrors("save_exchange_data_part_zr", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}
