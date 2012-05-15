//#include "../defines.h"
#include "../gpu.h"
#include "b-l.h"

//****************************
//TODO: copy to another place
//****************************

// Присвоение начальных условий
void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy]=def.porosity[0];
					HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.Background_Sn;

					double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
					                  + def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy]) * (def.g_const) * (def.hy);

					if (j == 0)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_g_dy;
					}

					/*
					// нагнетательная скважина
					if (is_injection_well(i, j, k, def))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Injection_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}

					// добывающая скважина
					if (is_output_well(i, j, k, def))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Production_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}
					*/

					HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					test_S(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_positive(HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}
//****************************
//TODO: copy to another place
//****************************

// Расчет относительных проницаемостей в точке
__device__ void device_assing_k(double* k_w, double* k_n, double S_w)
{
	/*
	// SPE-постановка
	double S_wc = 0.2;
	double S_or = 0.2;
	double S_e = (S_w - S_wc) / (1. - S_wc - S_or);

	*k_w = S_e * S_e;
	*k_n = (1. - S_e) * (1. - S_e);

	if (S_w < S_wc)
	{
		*k_w = 0.;
		*k_n = 1.;
	}

	if (S_w > (1 - S_or))
	{
		*k_w = 1.;
		*k_n = 0.;
	}
	*/
	
	// постановка ИПМ
	double S_sv = 0.1;
	double S_zv = 0.8;
	double S_1 = 0.70324;

	if ((S_sv<=S_w) && (S_w<=S_zv))
		*k_n=((S_zv-S_w)/(S_zv-S_sv))*((S_zv-S_w)/(S_zv-S_sv));
	else
		if ((0<=S_w) && (S_w<=S_sv))
			*k_n=1.;
		else //S_zv<S<=1
			*k_n=0.;

	if ((S_sv<=S_w) && (S_w<=S_1))
		*k_w=((S_w-S_sv)/(S_zv-S_sv))*((S_w-S_sv)/(S_zv-S_sv));
	else
		if ((0<=S_w) && (S_w<=S_sv))
			*k_w=0.;
		else 
			if ((S_1<=S_w) && (S_w<=S_zv))
				*k_w=0.8*pow((S_w-S_sv)/(S_zv-S_sv), 0.5);
			else//S_zv<S<=1
				*k_w=1.;

	device_test_S(*k_n, __FILE__, __LINE__);
	device_test_S(*k_w, __FILE__, __LINE__);
}

// Расчет плотностей, давления NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_P_Xi_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		double k_w=0., k_n=0.;
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		device_assing_k(&k_w, &k_n, 1. - DevArraysPtr.S_n[local]);

		DevArraysPtr.P_n[local] = DevArraysPtr.P_w[local];
		DevArraysPtr.Xi_w[local] = -1 * (DevArraysPtr.K[local]) * k_w / gpu_def->mu_w;
		DevArraysPtr.Xi_n[local] = -1 * (DevArraysPtr.K[local]) * k_n / gpu_def->mu_n;

		device_test_positive(DevArraysPtr.P_n[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[local], __FILE__, __LINE__);
	}
}


// Метод Ньютона для каждой точки сетки (независимо от остальных точек)
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx) - 1) && (j < gpu_def->locNy - 1) && (k < (gpu_def->locNz)) && (i != 0) && (j != 0) && (((k != 0) && (k != (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2)))
	{
		double A1 = DevArraysPtr.roS_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)];
		double A2 = DevArraysPtr.roS_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)];
		double a = gpu_def->beta_w * (gpu_def->beta_n);
		double b = gpu_def->beta_w + gpu_def->beta_n - A2 * (gpu_def->beta_w) / (gpu_def->ro0_n) - A1 * (gpu_def->beta_n) / (gpu_def->ro0_w);
		double c = 1. - A2 / gpu_def->ro0_n  - A1 / gpu_def->ro0_w;
		double D = b * b - 4. * a * c;
		double P1 = gpu_def->P_atm + (-1. * b + sqrt(D)) / (2. * a);
		double P2 = gpu_def->P_atm + (-1. * b - sqrt(D)) / (2. * a);

		if (P1 < 0.)
		{
			DevArraysPtr.P_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = P2;
		}
		else
		{
			DevArraysPtr.P_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = P1;
		}

		DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = DevArraysPtr.roS_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] / (gpu_def->ro0_n * (1 + gpu_def->beta_n * (DevArraysPtr.P_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - gpu_def->P_atm)));

		device_test_positive(DevArraysPtr.P_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
	}
}

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
__global__ void Border_S_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if (((i == 0) || (i == (gpu_def->locNx) - 1) || (j == 0) || (j == (gpu_def->locNy) - 1) ||
		(((k == 0) || (k == (gpu_def->locNz) - 1)) && ((gpu_def->locNz) >= 2))) && (device_is_active_point(i, j, k) == 1))
	{
		int i1 = i, j1 = j, k1 = k;

		device_set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1);

		DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = DevArraysPtr.S_n[i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy)];

		/*
		// В центре резервуара находится нагнетательная скважина
		if (device_is_injection_well(i, j, k, def))
		{
			DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = INJECTION_WELL_Sn;
		}

		// В центре резервуара находится добывающая скважина
		if (device_is_output_well(i, j, k, def))
		{
			DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] = OUTPUT_WELL_Sn;
		}
		*/

		device_test_S(DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)], __FILE__, __LINE__);
	}
}

__global__ void Border_P_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if (((i == 0) || (i == (gpu_def->locNx) - 1) || (j == 0) || (j == (gpu_def->locNy) - 1) ||
		(((k == 0) || (k == (gpu_def->locNz) - 1)) && ((gpu_def->locNz) >= 2))) && (device_is_active_point(i, j, k) == 1))
	{
		int i1 = i, j1 = j, k1 = k;
		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		device_set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1);

		if ((j != 0) && (j != (gpu_def->locNy) - 1))
		{
			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy)];
		}
		//else if(j == 0)
		//	DevArraysPtr.P_w[local] = gpu_def->P_atm;
		else
		{
			double ro_g_dy = (gpu_def->ro0_n * DevArraysPtr.S_n[local]
							  + gpu_def->ro0_w * (1 - DevArraysPtr.S_n[local])) * (DevArraysPtr.m[ local]) * (gpu_def->g_const) * (gpu_def->hy);

			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy)] + ro_g_dy;//DevArraysPtr.ro_w[i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy)] * (gpu_def->g_const) * (gpu_def->hy);
		}

	
		// В центре резервуара находится нагнетающая скважина
		if (device_is_injection_well(i, j, k))
		//if (((i == 0) && (j == 0)) || ((i == 1) && (j == 0)) || ((i == 0) && (j == 1)))
		{
			DevArraysPtr.P_w[local] = gpu_def->InjWell_Pw;
		}

		// В центре резервуара находится добывающая скважина
		if (device_is_output_well(i, j, k))
		//if (((i == gpu_def->Nx - 1) && (j == gpu_def->Ny - 1)) || ((i == gpu_def->Nx - 1) && (j == gpu_def->Ny - 2)) || ((i == gpu_def->Nx - 2) && (j == gpu_def->Ny - 1)))
		{
			DevArraysPtr.P_w[local] = gpu_def->OutWell_Pw;
		}
	

		device_test_positive(DevArraysPtr.P_w[local], __FILE__, __LINE__);
	}
}