//#include "../defines.h"
#include "../gpu.h"

// –асчет плотностей, давлени€ NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)) && (j<(locN.y)) && (k<(locN.z)) && (device_is_active_point(i, j, k, locN, rank, parts)==1))
	{
		int media = DevArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double S_n = DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double P_w = DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		double S_e = (1.- S_n - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]);
		double k_w = pow(S_e, (2. + 3. * (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]);
		double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));
		
		DevArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = P_w;
		DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = -1 * (DevArraysPtr.K[i+j*(locN.x)+k*(locN.x)*(locN.y)]) * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = -1 * (DevArraysPtr.K[i+j*(locN.x)+k*(locN.x)*(locN.y)]) * k_n / (*gpu_def).mu_n;
		DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w - (*gpu_def).P_atm));

		device_test_positive(DevArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}


// ћетод Ќьютона дл€ каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr, localN locN) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<(locN.x)-1) && (j<locN.y-1) && (k<(locN.z)) && (i!=0) && (j!=0) && (((k!=0) && (k!=(locN.z)-1)) || ((locN.z)<2)))
	{
		double F1, F2, F1P, F2P, F1S, F2S, det;
		double S_n=DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double P_w=DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		//S_e = (1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]);

		//AAA = pow(S_e, ((-1 / (*gpu_def).lambda[media]) - 1));
		F1 = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm)) * (1 - S_n) - DevArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		F2 = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w - (*gpu_def).P_atm)) * S_n - DevArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		F1P = (*gpu_def).ro0_w * ((*gpu_def).beta_w) * (1 - S_n);
		F2P = (*gpu_def).ro0_n * ((*gpu_def).beta_n) * S_n;
		F1S = (-1) * (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		F2S = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w - (*gpu_def).P_atm));

		det = F1P * F2S - F1S * F2P;

		DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = P_w - (1 / det) * (F2S * F1 - F1S * F2);
		DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = S_n - (1 / det) * (F1P * F2 - F2P * F1);

		device_test_positive(DevArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	}
}