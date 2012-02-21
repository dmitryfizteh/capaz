//#include "../defines.h"
#include "../gpu.h"
//#include "two-phase.h"

// –асчет плотностей, давлени€ NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr, consts def) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)) && (j<((*gpu_def).locNy)) && (k<((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def)==1))
	{
		int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double S_n = DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P_w = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		double S_e = (1.- S_n - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]);
		double k_w = pow(S_e, (2. + 3. * (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]);
		double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + (*gpu_def).lambda[media]) / (*gpu_def).lambda[media]));
		double P_k = (*gpu_def).P_d[media] * pow((1. - S_n - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media]), -1. / (*gpu_def).lambda[media]);

		DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = P_w + P_k;
		DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = -1 * (*gpu_def).K[media] * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = -1 * (*gpu_def).K[media] * k_n / (*gpu_def).mu_n;
		DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm));

		device_test_positive(DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}


// ћетод Ќьютона дл€ каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)-1) && (j<(*gpu_def).locNy-1) && (k<((*gpu_def).locNz)) && (i!=0) && (j!=0) && (((k!=0) && (k!=((*gpu_def).locNz)-1)) || (((*gpu_def).locNz)<2)))
	{
		int media = DevArraysPtr.media[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;
		double S_n=DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P_w=DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		S_e = (1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]);
		P_k = (*gpu_def).P_d[media] * pow(S_e, -1 / (*gpu_def).lambda[media]);
		AAA = pow(S_e, ((-1 / (*gpu_def).lambda[media]) - 1));
		F1 = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm)) * (1 - S_n) - DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		F2 = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm)) * S_n - DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		PkS = AAA * (*gpu_def).P_d[media] / ((*gpu_def).lambda[media] * (1 - (*gpu_def).S_wr[media]));
		F1P = (*gpu_def).ro0_w * ((*gpu_def).beta_w) * (1 - S_n);
		F2P = (*gpu_def).ro0_n * ((*gpu_def).beta_n) * S_n;
		F1S = (-1) * (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		F2S = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w + P_k - (*gpu_def).P_atm + (S_n * PkS)));

		det = F1P * F2S - F1S * F2P;

		DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = P_w - (1 / det) * (F2S * F1 - F1S * F2);
		DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = S_n - (1 / det) * (F1P * F2 - F2P * F1);

		device_test_positive(DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// «апись начальных условий дл€ каждой точки сетки (независимо от остальных точек)
__global__ void data_initialization(ptr_Arrays DefArraysPtr, long int* t, consts def)
{
	*t = 0;

	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if(device_is_active_point(i, j, k, def))
	{
		// ѕреобразование локальных координат процессора к глобальным
		int I = device_local_to_global(i, 'x', def);

		// ≈сли точка на верхней границе, не далее ((*gpu_def).source) точек от центра,
		// то в ней начальна€ насыщенность. »наче, нулева€
		if ((j==0) && (I>=((*gpu_def).Nx)/2-((*gpu_def).source)) && (I<=((*gpu_def).Nx)/2+((*gpu_def).source)) && (k>=((*gpu_def).Nz)/2-((*gpu_def).source)) && (k<=((*gpu_def).Nz)/2+((*gpu_def).source)))
			DefArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=(*gpu_def).S_n_gr;
		else
			DefArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=0;

		if(j == 0)
			DefArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=(*gpu_def).P_atm;
		else
			DefArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=DefArraysPtr.P_w[i+(j-1)*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy] + cu_ro_eff_gdy(DefArraysPtr, i, j-1, k, def);
							
		DefArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=0;

		DefArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1. + ((*gpu_def).beta_w) * (DefArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));

		///!!!! Ќе учитываютс€ капилл€рные силы! »ли надо считать перед этим шагом P_n
		DefArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * (DefArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));

		/*
		if ((DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=((*gpu_def).NX)/2.*((*gpu_def).h1)) && (DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=4.*((*gpu_def).NX)/5.*((*gpu_def).h1)))
			if ((DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=2./5.*(*gpu_def).locNy*((*gpu_def).h2)) && (DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=(-1.)*DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]/4.+2./5.*(*gpu_def).locNy*((*gpu_def).h2)))
				DefArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=1;

		if ((DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=((*gpu_def).NX)/5.*((*gpu_def).h1)) && (DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=2.*((*gpu_def).NX)/5.*((*gpu_def).h1)))
			if ((DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=4./5.*(*gpu_def).locNy*((*gpu_def).h2)) && (DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=3./5.*(*gpu_def).locNy*((*gpu_def).h2)))
				DefArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=1;
				*/
					
		/*
		if ((DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=2.*((*gpu_def).NX)/5.*((*gpu_def).h1)) && (DefArraysPtr.x[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=3.*((*gpu_def).NX)/5.*((*gpu_def).h1)))
			if ((DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]>=1./10.*(*gpu_def).locNy*((*gpu_def).h2)) && (DefArraysPtr.y[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]<=3./10.*(*gpu_def).locNy*((*gpu_def).h2)))
				DefArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=1;
		*/

		device_test_nan(DefArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DefArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DefArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
	}
}
