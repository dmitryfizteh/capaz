//#include "../defines.h"
#include "../gpu.h"

// –асчет плотностей, давлени€ NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr, consts def) 
{
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;

	if ((i<((*gpu_def).locNx)) && (j<((*gpu_def).locNy)) && (k<((*gpu_def).locNz)) && (device_is_active_point(i, j, k, def)==1))
	{
		double k_w, k_n;
		double S = 1 - DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		double S_wc=0.2;
//		double S_wi=0.2;
		double S_or=0.2;
		double S_e=(S-S_wc)/(1-S_wc-S_or);

		k_w=S_e*S_e;
		k_n=(1-S_e)*(1-S_e);

		//krw(Sor) = kro(Swc) = 1.0

		if (S<S_wc)
		{
			k_w=0;
			k_n=1;
		}

		if (S>(1-S_or))
		{
			k_w=1;
			k_n=0;
		}
		
		DevArraysPtr.P_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		DevArraysPtr.Xi_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = -1 * (DevArraysPtr.K[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = -1 * (DevArraysPtr.K[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]) * k_n / (*gpu_def).mu_n;
		DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));

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
		double A1 = DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double A2 = DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double a = (*gpu_def).beta_w * ((*gpu_def).beta_n);
		double b = (*gpu_def).beta_w + (*gpu_def).beta_n - A2 * ((*gpu_def).beta_w) / ((*gpu_def).ro0_n) - A1 * ((*gpu_def).beta_n) / ((*gpu_def).ro0_w);
		double c = 1 - A2 / (*gpu_def).ro0_n  - A1 / (*gpu_def).ro0_w;
		double D = b*b - 4*a*c;
		double P1 = (*gpu_def).P_atm + (-1 * b + sqrt(D)) / (2 * a);
		double P2 = (*gpu_def).P_atm + (-1 * b - sqrt(D)) / (2 * a);

		if (P1 < 0)
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = P2;
		else
			DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = P1;

		DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] / ((*gpu_def).ro0_n * (1 + (*gpu_def).beta_n * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)]-(*gpu_def).P_atm)));

		if ( DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] < 0.2)
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0.2;

		if ( DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] > 0.8)
			DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = 0.8;

		/*
		double F1, F2, F1P, F2P, F1S, F2S, det;
		double S_n=DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		double P_w=DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		//S_e = (1 - S_n - (*gpu_def).S_wr[media]) / (1 - (*gpu_def).S_wr[media]);

		//AAA = pow(S_e, ((-1 / (*gpu_def).lambda[media]) - 1));
		F1 = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm)) * (1 - S_n) - DevArraysPtr.roS_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];
		F2 = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w - (*gpu_def).P_atm)) * S_n - DevArraysPtr.roS_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)];

		F1P = (*gpu_def).ro0_w * ((*gpu_def).beta_w) * (1 - S_n);
		F2P = (*gpu_def).ro0_n * ((*gpu_def).beta_n) * S_n;
		F1S = (-1) * (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (P_w - (*gpu_def).P_atm));
		F2S = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (P_w - (*gpu_def).P_atm));

		det = F1P * F2S - F1S * F2P;

		DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = P_w - (1 / det) * (F2S * F1 - F1S * F2);
		DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = S_n - (1 / det) * (F1P * F2 - F2P * F1);
		*/

		device_test_positive(DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// ƒавление на нагнетающей скважине
__device__ double device_Injection_well_P(ptr_Arrays DevArraysPtr, int i, int j, int k, consts def)
{
	// 10000psi in Pa
	return def.P_atm+10000000;//68947572.9;
}

// ƒавление на добывающей скважине
__device__ double device_Production_well_P(ptr_Arrays DevArraysPtr, int i, int j, int k, consts def)
{
	// 4000psi in Pa
	return def.P_atm;//27579029.16;
}

// ѕрисвоение начальных условий дл€ каждой точки сетки (независимо от остальных точек)
__global__ void data_initialization(ptr_Arrays DevArraysPtr, long int* t, consts def)
{
	*t = 0;
	int i=threadIdx.x+blockIdx.x*blockDim.x;
	int j=threadIdx.y+blockIdx.y*blockDim.y;
	int k=threadIdx.z+blockIdx.z*blockDim.z;
	if(device_is_active_point(i, j, k, def))
	{
		// ѕреобразование локальных координат процессора к глобальным
		int I = device_local_to_global(i, 'x', def);

		DevArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=0;
		DevArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=0.7;
		//DevArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy] =0.3 + 0.3 * j / (*gpu_def).Ny;

		double ro_g_dy = ((*gpu_def).ro0_n * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] 
		+ (*gpu_def).ro0_w * (1 - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) * ((*gpu_def).m[(DevArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy])]) * ((*gpu_def).g_const) * ((*gpu_def).hy);

		// 6000 pound per square inch = 41 368 543.8 ѕаскал€
		if(j == 0)
			DevArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=(*gpu_def).P_atm + 5000000;// 50368543.8;
		else
			DevArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy]=DevArraysPtr.P_w[i+(j-1)*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy] + ro_g_dy;

						
		// ¬ центре резервуара находитс€ нагнетающа€ скважина
		if ((i==(*gpu_def).Nx/2) && (j==(*gpu_def).Ny-3) && (k==(*gpu_def).Nz/2))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = device_Injection_well_P(DevArraysPtr, i, j, k, def);
			//DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.5;
		}

		// ¬ центре резервуара находитс€ добывающа€ скважина
		if ((i==(*gpu_def).Nx-3) && (j==3) && (k==(*gpu_def).Nz-3))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = device_Production_well_P(DevArraysPtr, i, j, k, def);
			//DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.5;
		}
						

		DevArraysPtr.ro_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1. + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i+j*((*gpu_def).locNx)+k*((*gpu_def).locNx)*((*gpu_def).locNy)] - (*gpu_def).P_atm));

		device_test_nan(DevArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.P_w[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.media[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy], __FILE__, __LINE__);
	}
}

