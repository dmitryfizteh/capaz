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


// �������������� ��������� ��������� ���������� � ����������
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������ 
// (���� 2 ������ � ����� ������,�� +2 �����). 
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� rank==0)
__device__ int device_i_to_I(int i, consts def)
{
	int I;
	if ((*gpu_def).rank <= ((*gpu_def).Nx)%(*gpu_def).sizex)
	{
		if((*gpu_def).rank==0)
			I=i;
		else
			I=(((*gpu_def).Nx)/(*gpu_def).sizex+1)*(*gpu_def).rank+i-1;
	}
	else
		I=(((*gpu_def).Nx)/(*gpu_def).sizex+1)*(*gpu_def).rank-((*gpu_def).rank-((*gpu_def).Nx)%(*gpu_def).sizex)+i-1;
	return I;
}

// �������� �� ����� �������� (�.�. �� ��������������� ������ ��� ������ �� ��������)
__device__ int device_is_active_point(int i, int j, int k, consts def)
{
	if(((*gpu_def).rank % (*gpu_def).sizex != 0 && i == 0) || ((*gpu_def).rank % (*gpu_def).sizex != (*gpu_def).sizex - 1 && i == (*gpu_def).locNx - 1)
		|| (((*gpu_def).rank % (((*gpu_def).sizex) * ((*gpu_def).sizey))) / (*gpu_def).sizex != 0 && j == 0)	|| (((*gpu_def).rank % (((*gpu_def).sizex) * ((*gpu_def).sizey))) / (*gpu_def).sizex != (*gpu_def).sizey - 1 && j == (*gpu_def).locNy - 1)
		|| ((((*gpu_def).rank / ((*gpu_def).sizex) / ((*gpu_def).sizey) != 0 && k == 0) || ((*gpu_def).rank / ((*gpu_def).sizex) / ((*gpu_def).sizey) == (*gpu_def).sizez - 1 && k == (*gpu_def).locNz - 1)) && (*gpu_def).sizez > 1))
		return 0;
	else
		return 1;
}

// ������� ���������� "�����������" ���������
__device__ double cu_ro_eff_gdy(ptr_Arrays DevArraysPtr, int i, int j, int k, consts def)
{
	int media = DevArraysPtr.media[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

#ifdef THREE_PHASE
	double ro_g_dy = (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) *((*gpu_def).locNy)] * (1. - DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.S_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) 
		+ DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) *((*gpu_def).locNy)] * DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
	+ DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) *((*gpu_def).locNy)] * DevArraysPtr.S_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) * ((*gpu_def).m[media]) * ((*gpu_def).g_const) * ((*gpu_def).hy);
#else
	double ro_g_dy = (DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) *((*gpu_def).locNy)] * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] 
	+ DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) *((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) * ((*gpu_def).m[media]) * ((*gpu_def).g_const) * ((*gpu_def).hy);
#endif
	return ro_g_dy;
}

// ������ ����������, �������� NAPL P2 � Xi �� ���� ������ �����
void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_ro_Pn_Xi_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, def); 
	checkErrors("assign Pn, Xi, ro", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ �������� ���� P1 � ������������ NAPL S2 �� ���� ������ �����
void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int w=1;w<=def.newton_iterations;w++)
	{
		Newton_method_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr); 
	}
	checkErrors("assign Pw and Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ �������� � ������ ����� �����
__global__ void assign_u_kernel(ptr_Arrays DevArraysPtr, consts def) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)) && (j<((*gpu_def).locNy)) && (k<((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def)==1))
	{
		//CUPRINTF("assign u\n");
		double Xi_w = DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double Xi_n = DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P_w = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P_n = DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		if (((*gpu_def).Nx)>2)
		{
			if (i == 0)
			{
				DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (DevArraysPtr.P_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P_w) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (DevArraysPtr.P_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P_n) / ((*gpu_def).hx);
				//CUPRINTF("assign u=%e\n",DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]);
			}
			if (i == ((*gpu_def).locNx) - 1)
			{
				DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (P_w - DevArraysPtr.P_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / ((*gpu_def).hx);
				DevArraysPtr.ux_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (P_n - DevArraysPtr.P_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / ((*gpu_def).hx);
			}
			if ((i != 0) && (i != ((*gpu_def).locNx) - 1))
			{
				DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * ((DevArraysPtr.P_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.P_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx)));
				DevArraysPtr.ux_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * ((DevArraysPtr.P_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.P_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx)));
			}
		}
		else
		{
			DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
			DevArraysPtr.ux_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
		}
	
		if (((*gpu_def).Ny)>2)
		{
			if (j == 0)
			{
				DevArraysPtr.uy_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (
					(DevArraysPtr.P_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
				    - (P_w - DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (
					(DevArraysPtr.P_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (P_n - DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)) / (*gpu_def).hy;			}

			if (j == ((*gpu_def).locNy) - 1)
			{
				DevArraysPtr.uy_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (
					(P_w - DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)
					- (DevArraysPtr.P_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (*gpu_def).hy;

				DevArraysPtr.uy_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (
					(P_n - DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * j * (*gpu_def).hy)
					- (DevArraysPtr.P_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (*gpu_def).hy;
			}
			if ((j != 0) && (j != ((*gpu_def).locNy) - 1))
			{
				DevArraysPtr.uy_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (
					(DevArraysPtr.P_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (DevArraysPtr.P_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (2*(*gpu_def).hy);

				DevArraysPtr.uy_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (
					(DevArraysPtr.P_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j+1) * (*gpu_def).hy)
					- (DevArraysPtr.P_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] - DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] * (*gpu_def).g_const * (j-1) * (*gpu_def).hy)) / (2*(*gpu_def).hy);
			}
		}
		else
		{
			DevArraysPtr.uy_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
			DevArraysPtr.uy_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
		}

		if (((*gpu_def).Nz)>2)
		{
			if (k == 0)
			{
				DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - P_w) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - P_n) / ((*gpu_def).hz);
			}
			if (k == ((*gpu_def).Nz) - 1)
			{
				DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * (P_w - DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / ((*gpu_def).hz);
				DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * (P_n - DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / ((*gpu_def).hz);
			}
			if ((k != 0) && (k != ((*gpu_def).Nz) - 1))
			{
				DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_w * ((DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hz)));
				DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = Xi_n * ((DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hz)));
			}
		}
		else
		{
			DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
			DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0;
		}

	device_test_nan(DevArraysPtr.ux_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.ux_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uy_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uy_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	device_test_nan(DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ������ ��������� �� ���� ������ �����
void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	assign_u_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, def); 
	checkErrors("assign u", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

// ������ ro*S � ������ ����� ����� ������� ������������ ���������
__global__ void assign_roS_kernel_nr(ptr_Arrays DevArraysPtr, double t)
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;
	
	if ((i<((*gpu_def).locNx)-1) && (j<(*gpu_def).locNy-1) && (k<((*gpu_def).locNz)) && (i!=0) && (j!=0) && (((k!=0) && (k!=((*gpu_def).locNz)-1)) || (((*gpu_def).locNz)<2)))
	{
		int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		double q_w=0;
		double q_n=0;

#ifdef B_L
		/*
		double F_bl=0;
		// � ����� ������ ���� ����������� ��������
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			q_w=(*gpu_def).Q;
			q_n=0;
		}

		// � ������ ������� ���� ���������� ��������
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

		double S2 = DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double roS1 = DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - S2);
		double roS2 = DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * S2;
 		double P1 = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P2 = DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3;

		if (((*gpu_def).Nz)<2)
		{
			f3=0;
			g3=0;
		}
		else
		{
			z2 = -(DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - P1)/(*gpu_def).hz;
			z1 = -(P1 - DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hz;

			f3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                      (z1 + abs(z1))/2.0*(-1)* DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                      (z2 - abs(z2))/2.0*(-1)* DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hz;

			z2 = -(DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - P2)/(*gpu_def).hz;
			z1 = -(P2 - DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hz;

			g3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                      (z1 + abs(z1))/2.0*(-1)* DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                      (z2 - abs(z2))/2.0*(-1)* DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hz;
		}

		x2 = -(DevArraysPtr.P_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P1)/(*gpu_def).hx;
        x1 = -(P1 - DevArraysPtr.P_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hx;

        y2 = -(DevArraysPtr.P_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P1)/(*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (j+1) - DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * j);
        y1 = -(P1 - DevArraysPtr.P_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * j - DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (j-1));

        f1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                (x1 + abs(x1))/2.0*(-1)* DevArraysPtr.Xi_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                (x2 - abs(x2))/2.0*(-1)* DevArraysPtr.Xi_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hx;

        f2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                (y1 + abs(y1))/2.0*(-1)* DevArraysPtr.Xi_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                (y2 - abs(y2))/2.0*(-1)* DevArraysPtr.Xi_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hy;


        x2 = -(DevArraysPtr.P_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P2)/(*gpu_def).hx;
        x1 = -(P2 - DevArraysPtr.P_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hx;

        y2 = -(DevArraysPtr.P_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - P2)/(*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (j+1) - DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * j);
        y1 = -(P2 - DevArraysPtr.P_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hy + (*gpu_def).g_const * (DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * j - DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (j-1));

        g1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                (x1 + abs(x1))/2.0*(-1)* DevArraysPtr.Xi_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                (x2 - abs(x2))/2.0*(-1)* DevArraysPtr.Xi_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hx;

        g2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] -
                (y1 + abs(y1))/2.0*(-1)* DevArraysPtr.Xi_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] +
                (y2 - abs(y2))/2.0*(-1)* DevArraysPtr.Xi_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])/(*gpu_def).hy;

		DevArraysPtr.roS_w_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = roS1;
		DevArraysPtr.roS_n_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = roS2;
		DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = roS1 - ((*gpu_def).dt/(*gpu_def).m[media])*(q_w + f1 + f2 + f3);
		DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = roS2 - ((*gpu_def).dt/(*gpu_def).m[media])*(q_n + g1 + g2 + g3);

		device_test_positive(DevArraysPtr.roS_w_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ������ ro*S � ������ ����� �����
__global__ void assign_roS_kernel(ptr_Arrays DevArraysPtr, double t) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)-1) && (j<(*gpu_def).locNy-1) && (k<((*gpu_def).locNz)) && (i!=0) && (j!=0) && (((k!=0) && (k!=((*gpu_def).locNz))) || (((*gpu_def).locNz)<2)))
	{
		int local=i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy);
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
			divgrad1 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) * (DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1. - DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[local] * (1. - S_n) + DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1. - DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			divgrad2 = ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) * (DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[local] * S_n + DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * (DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)])) / (((*gpu_def).hz) * ((*gpu_def).hz));
			Tz1 = (DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uz_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2. * ((*gpu_def).hz));
			Tz2 = (DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uz_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2. * ((*gpu_def).hz));
		}

		divgrad1 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_w) / 2.) *
		((DevArraysPtr.ro_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		(DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) - 2 * DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - S_n) + DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (1 - DevArraysPtr.S_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)])) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		divgrad2 += ((*gpu_def).m[media] * ((*gpu_def).l) * ((*gpu_def).c_n) / 2.) *
		((DevArraysPtr.ro_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.S_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * S_n + DevArraysPtr.ro_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.S_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (((*gpu_def).hx) * ((*gpu_def).hx)) +
		(DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.S_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - 2 * DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * S_n + DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.S_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (((*gpu_def).hy) * ((*gpu_def).hy)));

		Tx1 = (DevArraysPtr.ro_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ux_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ux_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx));
		Ty1 = (DevArraysPtr.ro_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uy_w[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uy_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hy));
		Tx2 = (DevArraysPtr.ro_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ux_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.ux_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hx));
		Ty2 = (DevArraysPtr.ro_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uy_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - DevArraysPtr.ro_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * DevArraysPtr.uy_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) / (2 * ((*gpu_def).hy));

		double q_w=0;
		double q_n=0;

#ifdef B_L
		double F_bl=0;
		// � ����� ������ ���� ����������� ��������
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			q_w=(*gpu_def).Q;
			q_n=0;
		}

		// � ������ ������� ���� ���������� ��������
		if (((i==0) && (j==(*gpu_def).Ny-2)) || ((i==1) && (j==(*gpu_def).Ny-1)))
		{
			//int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			double S_e = (1. - DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]);
			double k_w = pow(S_e, (2. + 3. * ((*gpu_def).lambda[media])) / (*gpu_def).lambda[media]);
			double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));

			F_bl = (k_w/((*gpu_def).mu_w)) / (k_w/((*gpu_def).mu_w) + k_n/((*gpu_def).mu_n));
			q_w=-1 * (*gpu_def).Q * F_bl;
			q_n=-1 * (*gpu_def).Q * (1-F_bl);
		}
#endif

		if ((t < 2 * ((*gpu_def).dt)) || TWO_LAYERS)
		{
			A1 = roS_w + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = roS_n + (((*gpu_def).dt) / (*gpu_def).m[media]) * (q_w + divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1 + (2 * roS_w * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_w_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
			A2 = (2 * ((*gpu_def).dt) * ((*gpu_def).dt)) / ((*gpu_def).m[media] * (((*gpu_def).dt) + 2 * ((*gpu_def).tau))) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2 + (2 * roS_n * (*gpu_def).m[media] * ((*gpu_def).tau)) / (((*gpu_def).dt) * ((*gpu_def).dt)) + DevArraysPtr.roS_n_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] * (*gpu_def).m[media] * (((*gpu_def).dt) - 2 * ((*gpu_def).tau)) / (2 * ((*gpu_def).dt) * ((*gpu_def).dt)));
		}

		DevArraysPtr.roS_w_old[local] = roS_w;
		DevArraysPtr.roS_n_old[local] = roS_n;
		DevArraysPtr.roS_w[local] = A1;
		DevArraysPtr.roS_n[local] = A2;

		device_test_positive(DevArraysPtr.roS_w_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n_old[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ������ ro*S �� ���� ������ �����
void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def)
{
	#ifdef NR
		assign_roS_kernel_nr<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, t); 
	#else
		assign_roS_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, t);
	#endif
		checkErrors("assign roS", __FILE__, __LINE__);
		cudaPrintfDisplay(stdout, true);
}

// ��������� ������� �� S2
__global__ void Sn_boundary_kernel(ptr_Arrays DevArraysPtr, consts def) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)) && (j<((*gpu_def).locNy)) && (k<((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def)==1))
	{
		if ((i == 0) && ((*gpu_def).Nx>2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
		   DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		   return;
		}

		if ((i == ((*gpu_def).locNx) - 1) && (((*gpu_def).Nx)>2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			 return;
		}

		if ((j == ((*gpu_def).locNy) - 1) && (((*gpu_def).locNy)>2))
		{
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

			if (i==0) 
				DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			if (i==(*gpu_def).Nx - 1)
				DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			return;
		}
		
		if ((j==0) && (((*gpu_def).locNy)>2))
		{
			int I=device_i_to_I(i, gpu_def[0]);
			if ((I>=((*gpu_def).Nx)/2-((*gpu_def).source)) && (I<=((*gpu_def).Nx)/2+((*gpu_def).source)) && (k>=((*gpu_def).Nz)/2-((*gpu_def).source)) && (k<=((*gpu_def).Nz)/2+((*gpu_def).source)))
				DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).S_n_gr;
			else
				DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+(j+1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			return;
		}

		if ((k == 0) && ((*gpu_def).Nz > 2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)];
			return;
		}
		
		if ((k == (*gpu_def).Nz - 1) && ((*gpu_def).Nz > 2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)];
			return;
		}

		device_test_S(DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ��������� ������� �� P1
__global__ void Pw_boundary_kernel(ptr_Arrays DevArraysPtr, consts def) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)) && (j<((*gpu_def).locNy)) && (k<((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def)==1))
	{
		if ((i == 0) && ((*gpu_def).Nx > 2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]; 
			//return;
		}

		if ((i == ((*gpu_def).locNx) - 1) && ((*gpu_def).Nx>2) && (j>0) && (j<(*gpu_def).locNy - 1))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			//return;
		}

		if ((j == ((*gpu_def).locNy) - 1) && (((*gpu_def).locNy)>2))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+(j-1)*((*gpu_def).locNx)+k*((*gpu_def).locNx)*(*gpu_def).locNy] + cu_ro_eff_gdy(DevArraysPtr, i, j-1, k, def);
			if (i==0) 
				DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			if (i==(*gpu_def).Nx - 1)
				DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i-1+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
			//return;
		}

		if ((j==0) && (((*gpu_def).locNy)>2))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).P_atm;
			//return;
		}

		if ((k == 0) && (((*gpu_def).Nz)>2) && (j>0) && (j<((*gpu_def).locNy) - 1))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k+1)*((*gpu_def).locNx)*((*gpu_def).locNy)]; 
			//return;
		}

		if ((k == ((*gpu_def).locNz) - 1) && (((*gpu_def).locNz)>2) && (j>0) && (j<((*gpu_def).locNy) - 1))
		{
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+(k-1)*((*gpu_def).locNx)*((*gpu_def).locNy)];
			//return;
		}

		device_test_positive(DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ���������� ��������� �������
void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	Sn_boundary_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, def); 
	checkErrors("assign Sn", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	Pw_boundary_kernel<<<dim3(def.blocksX, def.blocksY, def.blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr, def); 
	checkErrors("assign Pw", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}	

// ������� �������� ������ � ������ �����
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy( HostArrayPtr, DevArrayPtr, ((*gpu_def).locNx)*((*gpu_def).locNy)*((*gpu_def).locNz)*sizeof(double), cudaMemcpyDeviceToHost );
	checkErrors("copy data to host", __FILE__, __LINE__);
}

// ������� �������� ������ ���� double � ������ ����������
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, consts def)
{
	cudaMemcpy( DevArrayPtr, HostArrayPtr, ((*gpu_def).locNx)*((*gpu_def).locNy)*((*gpu_def).locNz)*sizeof(double), cudaMemcpyHostToDevice );
	checkErrors("copy double data to device", __FILE__, __LINE__);
}

// ������� �������� ������ ���� int � ������ ����������
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, consts def)
{
	cudaMemcpy( DevArrayPtr, HostArrayPtr, ((*gpu_def).locNx)*((*gpu_def).locNy)*((*gpu_def).locNz)*sizeof(int), cudaMemcpyHostToDevice );
	checkErrors("copy int data to device", __FILE__, __LINE__);
}

// ��������� ������ ���������� ��� ������ ����� ��������� �������
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, consts def)
{
	cudaMalloc((void**) DevBuffer,  2 * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).P_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).P_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).S_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ro_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ro_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ux_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uy_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uz_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).ux_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uy_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).uz_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).Xi_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).Xi_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_w),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_n),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_w_old),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).roS_n_old),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
	cudaMalloc((void**) &((*ArraysPtr).media),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(int));

#ifdef B_L
	cudaMalloc((void**) &((*ArraysPtr).K),  ((*gpu_def).locNx) * ((*gpu_def).locNy) * ((*gpu_def).locNz) * sizeof(double));
#endif

	checkErrors("memory allocation", __FILE__, __LINE__);
}

// ����������� ������ ���������� �� ��� ������� ����� ��������� �������
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
	(*def).blocksX=((*gpu_def).locNx)/BlockNX;
	if (((*gpu_def).locNx)%BlockNX!=0)
		((*def).blocksX)++;
	(*def).blocksY=((*gpu_def).locNy)/BlockNY;
	if (((*gpu_def).locNy)%BlockNY!=0)
		((*def).blocksY)++;
	(*def).blocksZ=((*gpu_def).locNz)/BlockNZ;
	if (((*gpu_def).locNz)%BlockNZ!=0)
		((*def).blocksZ)++;

	consts* deff=new consts[1];
	deff[0] = (*def);
	cudaMemcpyToSymbol ( gpu_def, deff, sizeof ( consts ));
	checkErrors("constant memory copy", __FILE__, __LINE__);

	cudaDeviceProp devProp;
    cudaGetDeviceProperties ( &devProp, device );
        
	if (devProp.major < 2)
		printf ("\nError! Compute capability < 2, rank=%d\n",(*def).rank);
        
	if (!(*def).rank)
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


		// ������������ ������ ��������� ����� ��� ����������
		// sizeof(ptr_Arrays)/4 - ���������� ���������� � �����, �.�. 4 -������ ������ ���������
		printf ( "\nTotal NAPL_Filtration grid size : %d\n\n", devProp.totalGlobalMem/(sizeof(ptr_Arrays)*sizeof(double)/4) );
	}

		// ((*gpu_def).locNx)+2 ������ ��� 2NyNz �� ����� ������ ����������
	if ( ((*gpu_def).locNx+2)*((*gpu_def).locNy)*((*gpu_def).locNz) > (devProp.totalGlobalMem/(sizeof(ptr_Arrays)*sizeof(double)/4)))
		printf ("\nError! Not enough memory at GPU, rank=%d\n",(*def).rank);
	fflush( stdout);

	// �������������� ���������� cuPrintf ��� ������ ������ �� �������
	// ����� �� kernel
	cudaPrintfInit();
}

// ����������� ����������
void device__finalization(void)
{
	// ������������� ���������� cuPrintf ��� ������ ������ �� �������
	// ����� �� kernel
	cudaPrintfEnd();
}

__global__ void load_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.y+blockIdx.y*blockDim.y;

	if (j<(*gpu_def).locNy && k<((*gpu_def).locNz))
	{
		DevBuffer[j+((*gpu_def).locNy)*k]=DevArrayPtr[1+((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k];
		DevBuffer[j+((*gpu_def).locNy)*k+((*gpu_def).locNy)*((*gpu_def).locNz)]=DevArrayPtr[((*gpu_def).locNx)-2+((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k];

		device_test_nan(DevBuffer[j+((*gpu_def).locNy)*k], __FILE__, __LINE__);
		device_test_nan(DevBuffer[j+((*gpu_def).locNy)*k+((*gpu_def).locNy)*((*gpu_def).locNz)], __FILE__, __LINE__);
	}
}


void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	load_exchange_data_kernel<<<dim3(def.blocksY, def.blocksZ), dim3(BlockNY,BlockNZ)>>>(DevArrayPtr, DevBuffer); 
	checkErrors("load_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	cudaMemcpy( HostBuffer, DevBuffer, 2 * (def.Ny) * (def.Nz) * sizeof(double), cudaMemcpyDeviceToHost );
	checkErrors("copy data to host", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}

__global__ void save_exchange_data_kernel(double* DevArrayPtr, double* DevBuffer)
{
	int j=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.y+blockIdx.y*blockDim.y;

	if (j<((*gpu_def).locNy) && k<((*gpu_def).locNz))
	{
		if ((*gpu_def).rank!=0)
		{
			DevArrayPtr[((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k]=DevBuffer[j+((*gpu_def).locNy)*k];
			device_test_nan(DevArrayPtr[((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k], __FILE__, __LINE__);
		}
		if ((*gpu_def).rank!=(*gpu_def).sizex-1)
		{
			DevArrayPtr[((*gpu_def).locNx)-1+((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k]=DevBuffer[j+((*gpu_def).locNy)*k+((*gpu_def).locNy)*((*gpu_def).locNz)];
			device_test_nan(DevArrayPtr[((*gpu_def).locNx)-1+((*gpu_def).locNx)*j+((*gpu_def).locNx)*((*gpu_def).locNy)*k], __FILE__, __LINE__);
		}
	}
}

void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	cudaMemcpy( DevBuffer, HostBuffer, 2*((*gpu_def).locNy)*((*gpu_def).locNz)*sizeof(double), cudaMemcpyHostToDevice );
	checkErrors("copy data to device", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);

	save_exchange_data_kernel<<<dim3(def.blocksY, def.blocksZ), dim3(BlockNY,BlockNZ)>>>(DevArrayPtr, DevBuffer); 
	checkErrors("save_exchange_data", __FILE__, __LINE__);
	cudaPrintfDisplay(stdout, true);
}