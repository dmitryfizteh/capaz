#include "../defines.h"
#include "../gpu.h"
#include "b-l.h"

// явл€етс€ ли точка активной (т.е. не предназначенной только дл€ обмена на границах)
__device__ int device_is_active_point_local(int i, int localNx, int rank, int size)
{
	if((rank!=0 && i==0) || (rank!=size-1 && i==localNx-1))
		return 0;
	else
		return 1;
}

// –асчет плотностей, давлени€ NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr, int localNx, int rank, int size) 
{
	int blockIdxz=blockIdx.y / BlockNY;
	int blockIdxy=blockIdx.y % BlockNY;
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.z+blockIdxz*blockDim.z;
	int j=threadIdx.y+blockIdxy*blockDim.y;

	if ((i<localNx) && (j<((*gpu_def).Ny)) && (k<((*gpu_def).Nz)) && (device_is_active_point_local(i, localNx, rank, size)==1))
	{
		int media = DevArraysPtr.media[i+j*localNx+k*localNx*((*gpu_def).Ny)];
		double S_n = DevArraysPtr.S_n[i+j*localNx+k*localNx*((*gpu_def).Ny)];
		double P_w = DevArraysPtr.P_w[i+j*localNx+k*localNx*((*gpu_def).Ny)];

		double S_e = (1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]);
		double k_w = pow(S_e, (2. + 3. * (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]);
		double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));
		double P_k = (*gpu_def).P_d[media] * pow((1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]), -1 / (*gpu_def).lambda[media]);

		DevArraysPtr.P_n[i+j*localNx+k*localNx*((*gpu_def).Ny)] = P_w + P_k;
		DevArraysPtr.Xi_w[i+j*localNx+k*localNx*((*gpu_def).Ny)] = -1 * (*gpu_def).K[media] * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i+j*localNx+k*localNx*((*gpu_def).Ny)] = -1 * (*gpu_def).K[media] * k_n / (*gpu_def).mu_n;
		DevArraysPtr.ro_w[i+j*localNx+k*localNx*((*gpu_def).Ny)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i+j*localNx+k*localNx*((*gpu_def).Ny)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm));
	}
}

// –асчет плотностей, давлени€ NAPL P2 и Xi во всех точках сетки
void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ)
{
	assign_ro_Pn_Xi_kernel<<<dim3(blocksX,blocksY*blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr,localNx,rank,size); 
	checkErrors("assign Pn, Xi, ro");
}

// ћетод Ќьютона дл€ каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr, int localNx) 
{
	int blockIdxz=blockIdx.y / BlockNY;
	int blockIdxy=blockIdx.y % BlockNY;
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int k=threadIdx.z+blockIdxz*blockDim.z;
	int j=threadIdx.y+blockIdxy*blockDim.y;

	if ((i<localNx-1) && (j<(*gpu_def).Ny-1) && (k<(*gpu_def).Nz) && (i!=0) && (j!=0) && (((k!=0) && (k!=(*gpu_def).Nz-1)) || ((*gpu_def).Nz<2)))
	{
		int media = DevArraysPtr.media[i+j*localNx+k*localNx*((*gpu_def).Ny)];
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;
		double S_n=DevArraysPtr.S_n[i+j*localNx+k*localNx*((*gpu_def).Ny)];
		double P_w=DevArraysPtr.P_w[i+j*localNx+k*localNx*((*gpu_def).Ny)];

		S_e = (1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]);
		P_k = (*gpu_def).P_d[media] * pow(S_e, -1 / (*gpu_def).lambda[media]);
		AAA = pow(S_e, ((-1 / (*gpu_def).lambda[media]) - 1));
		F1 = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm)) * (1 - S_n) - DevArraysPtr.roS_w[i+j*localNx+k*localNx*((*gpu_def).Ny)];
		F2 = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm)) * S_n - DevArraysPtr.roS_n[i+j*localNx+k*localNx*((*gpu_def).Ny)];

		PkS = AAA * (*gpu_def).P_d[media] / ((*gpu_def).lambda[media] * (1 - (*gpu_def).S_wr[media]));
		F1P = (*gpu_def).ro0_w * ((*gpu_def).beta_w) * (1 - S_n);
		F2P = (*gpu_def).ro0_n * ((*gpu_def).beta_n) * S_n;
		F1S = (-1) * (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		F2S = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm + (S_n * PkS)));

		det = F1P * F2S - F1S * F2P;

		DevArraysPtr.P_w[i+j*localNx+k*localNx*((*gpu_def).Ny)] = P_w - (1 / det) * (F2S * F1 - F1S * F2);
		DevArraysPtr.S_n[i+j*localNx+k*localNx*((*gpu_def).Ny)] = S_n - (1 / det) * (F1P * F2 - F2P * F1);
	}
}

// –асчет давлени€ воды P1 и насыщенности NAPL S2 во всех точках сетки
void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ)
{
	Newton_method_kernel<<<dim3(blocksX,blocksY*blocksZ), dim3(BlockNX,BlockNY,BlockNZ)>>>(DevArraysPtr,localNx); 
	checkErrors("assign Pw and Sn");
}
