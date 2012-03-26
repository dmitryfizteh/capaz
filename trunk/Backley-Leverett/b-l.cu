//#include "../defines.h"
#include "../gpu.h"
#include "b-l.h"

//****************************
//TODO: copy to another place
//****************************

// Давление на нагнетающей скважине
double Injection_well_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	// 10000psi in Pa
	return INJECTION_WELL_Pw; //68947572.9;
}

// Давление на добывающей скважине
double Production_well_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	// 4000psi in Pa
	return OUTPUT_WELL_Pw;//27579029.16;
}

// Присвоение начальных условий
void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					// Преобразование локальных координат процессора к глобальным
					int I = local_to_global(i, 'x', def);

					//HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy] = 0;
					HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = BACKGROUND_Sn;
					//HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy] =0.3 + 0.3 * j / def.Ny;

					double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
					                  + def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[i + j * def.locNx + k * def.locNx * def.locNy]) * (def.g_const) * (def.hy);

					// 6000 pound per square inch = 41 368 543.8 Паскаля
					if (j == 0)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = BACKGROUND_Pw; //+ 5000000;    // 50368543.8;
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_g_dy;
					}

					
					// В центре резервуара находится нагнетающая скважина
					if ((i == def.Nx / 2) && (j == def.Ny / 2) && (k == def.Nz / 2))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Injection_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}

					// В центре резервуара находится добывающая скважина
					if ((i == def.Nx - 3) && (j == def.Ny / 2) && (k == def.Nz - 3))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Production_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}
					

					HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					test_nan(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}

//****************************
//TODO: copy to another place
//****************************


// Расчет плотностей, давления NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx)) && (j < ((*gpu_def).locNy)) && (k < ((*gpu_def).locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		double k_w, k_n;
		double S = 1 - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

		double S_wc = 0.2;
		//		double S_wi=0.2;
		double S_or = 0.2;
		double S_e = (S - S_wc) / (1 - S_wc - S_or);

		k_w = S_e * S_e;
		k_n = (1 - S_e) * (1 - S_e);

		//krw(Sor) = kro(Swc) = 1.0

		if (S < S_wc)
		{
			k_w = 0;
			k_n = 1;
		}

		if (S > (1 - S_or))
		{
			k_w = 1;
			k_n = 0;
		}

		DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = -1 * (DevArraysPtr.K[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = -1 * (DevArraysPtr.K[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]) * k_n / (*gpu_def).mu_n;
		DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));

		device_test_positive(DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}


// Метод Ньютона для каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx) - 1) && (j < (*gpu_def).locNy - 1) && (k < ((*gpu_def).locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != ((*gpu_def).locNz) - 1)) || (((*gpu_def).locNz) < 2)))
	{
		double A1 = DevArraysPtr.roS_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double A2 = DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double a = (*gpu_def).beta_w * ((*gpu_def).beta_n);
		double b = (*gpu_def).beta_w + (*gpu_def).beta_n - A2 * ((*gpu_def).beta_w) / ((*gpu_def).ro0_n) - A1 * ((*gpu_def).beta_n) / ((*gpu_def).ro0_w);
		double c = 1 - A2 / (*gpu_def).ro0_n  - A1 / (*gpu_def).ro0_w;
		double D = b * b - 4 * a * c;
		double P1 = (*gpu_def).P_atm + (-1 * b + sqrt(D)) / (2 * a);
		double P2 = (*gpu_def).P_atm + (-1 * b - sqrt(D)) / (2 * a);

		if (P1 < 0)
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = P2;
		}
		else
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = P1;
		}

		DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] / ((*gpu_def).ro0_n * (1 + (*gpu_def).beta_n * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm)));

		if (DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] < 0.2)
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.2;
		}

		if (DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] > 0.8)
		{
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.8;
		}

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

		device_test_positive(DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

// Давление на нагнетательной скважине
__device__ double device_Injection_well_P(ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	// 10000psi in Pa
	return INJECTION_WELL_Pw;//(*gpu_def).P_atm + 10000000; //68947572.9;
}

// Давление на добывающей скважине
__device__ double device_Production_well_P(ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	// 4000psi in Pa
	return OUTPUT_WELL_Pw;//(*gpu_def).P_atm;//27579029.16;
}

// Присвоение начальных условий для каждой точки сетки (независимо от остальных точек)
__global__ void data_initialization_kernel(ptr_Arrays DevArraysPtr, long int* t)
{
	*t = 0;
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;
	if (device_is_active_point(i, j, k))
	{
		// Преобразование локальных координат процессора к глобальным
		int I = device_local_to_global(i, 'x');

		//DevArraysPtr.m[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = 0;
		DevArraysPtr.S_n[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = 0.7;
		//DevArraysPtr.S_n[i+j*(*gpu_def).locNx+k*(*gpu_def).locNx*(*gpu_def).locNy] =0.3 + 0.3 * j / (*gpu_def).Ny;

		double ro_g_dy = ((*gpu_def).ro0_n * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
		                  + (*gpu_def).ro0_w * (1 - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)])) * ((*gpu_def).m[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy]) * ((*gpu_def).g_const) * ((*gpu_def).hy);

		// 6000 pound per square inch = 41 368 543.8 Паскаля
		if (j == 0)
		{
			DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).P_atm; //+ 5000000;    // 50368543.8;
		}
		else
		{
			DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = DevArraysPtr.P_w[i + (j - 1) * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] + ro_g_dy;
		}

		/*
		// В центре резервуара находится нагнетательная скважина
		if ((i == (*gpu_def).Nx / 2) && (j == (*gpu_def).Ny / 2) && (k == (*gpu_def).Nz / 2))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = device_Injection_well_P(DevArraysPtr, i, j, k);
			//DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.5;
		}

		// В центре резервуара находится добывающая скважина
		if ((i == (*gpu_def).Nx - 3) && (j == (*gpu_def).Ny / 2) && (k == (*gpu_def).Nz - 3))
		{
			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = device_Production_well_P(DevArraysPtr, i, j, k);
			//DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = 0.5;
		}
		*/

		DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1. + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));

		device_test_nan(DevArraysPtr.S_n[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.m[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
	}
}

