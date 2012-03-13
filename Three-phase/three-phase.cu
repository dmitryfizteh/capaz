#include "../gpu.h"

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

					int media = HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy] = 0;
	/*				int j1 = def.locNy / 2;

					if (j < j1)
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_gr + (def.S_n_init - def.S_n_gr) * j / j1;
					}
					else
					*/
					if (j == 0)
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_gr;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_init;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_init;
					}

					if (j == 0)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_eff_gdy(HostArraysPtr, i, j - 1, k, def);
					}

					HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					///!!!! Не учитываются каппилярные силы! Или надо считать перед этим шагом P_w, P_g
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
					HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					test_nan(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}

// Расчет плотностей, давления NAPL P2 и Xi в каждой точке сетки (независимо от остальных точек)
__global__ void assign_ro_Pn_Xi_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < ((*gpu_def).locNx)) && (j < ((*gpu_def).locNy)) && (k < ((*gpu_def).locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int media = DevArraysPtr.media[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double k_w, k_g, k_n, P_k_nw, P_k_gn;
		double A = (*gpu_def).lambda[media];
		double S_w_e = (DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_gr[media]);
		double S_n_e = (DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).S_nr[media]) / (1. - (*gpu_def).S_nr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_nr[media]);
		double S_g_e = 1. - S_w_e - S_n_e;

		if (S_w_e <= (*gpu_def).S_w_range[0])
		{
			S_w_e = (*gpu_def).S_w_range[0];
			k_w = 0.;
		}
		else
		{
			k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)), (A - 1.) / A), 2.);
		}

		if (S_g_e <= (*gpu_def).S_g_range[0])
		{
			S_g_e = (*gpu_def).S_g_range[0];
			k_g = 0.;
		}
		else
		{
			k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);
		}

		if (S_n_e <= 0.)
		{
			k_n = 0.;
		}
		else
		{
			double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);
			double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A / (A - 1.)), (A - 1.) / A), 2.);
			k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e);
		}

		if (S_w_e <= (*gpu_def).S_w_range[1])
		{
			P_k_nw = ((*gpu_def).aw[0]) * S_w_e + ((*gpu_def).bw[0]);
		}
		else if (S_w_e >= (*gpu_def).S_w_range[2])
		{
			P_k_nw = ((*gpu_def).aw[1]) * S_w_e + ((*gpu_def).bw[1]);
		}
		else
		{
			P_k_nw = (*gpu_def).P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A);
		}

		if (S_g_e <= (*gpu_def).S_g_range[1])
		{
			P_k_gn = ((*gpu_def).ag[0]) * S_g_e + ((*gpu_def).bg[0]);
		}
		else if (S_g_e >= (*gpu_def).S_g_range[2])
		{
			P_k_gn = ((*gpu_def).ag[1]) * S_g_e + ((*gpu_def).bg[1]);
		}
		else
		{
			P_k_gn = (*gpu_def).P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A);
		}

		DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw;
		DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_gn;

		DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (-1.) * ((*gpu_def).K[media]) * k_w / (*gpu_def).mu_w;
		DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (-1.) * ((*gpu_def).K[media]) * k_n / (*gpu_def).mu_n;
		DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (-1.) * ((*gpu_def).K[media]) * k_g / (*gpu_def).mu_g;

		DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1 + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1 + ((*gpu_def).beta_n) * (DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_g * DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] / (*gpu_def).P_atm;

		device_test_positive(DevArraysPtr.P_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.P_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);

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
		int media = DevArraysPtr.media[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		double S_w_e, S_g_e, S_n_e, P_k_nw, P_k_gn, A, Sg, F1, F2, F3;
		double PkSw, PkSn, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sn, F2Sn, F3Sn, det;
		double a[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		for (int w = 1; w <= (*gpu_def).newton_iterations; w++)
		{
			S_w_e = (DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).S_wr[media]) / (1. - (*gpu_def).S_wr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_gr[media]);
			S_n_e = (DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).S_nr[media]) / (1. - (*gpu_def).S_nr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_nr[media]);
			S_g_e = 1. - S_w_e - S_n_e;
			A = (*gpu_def).lambda[media];                                                                                                                                                                                                                                                                  /*2*/

			// По краям интервала [0, 1] функции капиллярных давлений гладко заменяем линейными, производные меняются соответственно.
			// Описание можно посмотреть в файле mathcad.
			if (S_w_e <= (*gpu_def).S_w_range[1])
			{
				P_k_nw = ((*gpu_def).aw[0]) * S_w_e + ((*gpu_def).bw[0]);
				PkSw = ((*gpu_def).aw[0]);
			}
			else if (S_w_e >= (*gpu_def).S_w_range[2])
			{
				P_k_nw = ((*gpu_def).aw[1]) * S_w_e + ((*gpu_def).bw[1]);
				PkSw = ((*gpu_def).aw[1]);
			}
			else
			{
				P_k_nw = (*gpu_def).P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A);
				PkSw = (*gpu_def).P_d_nw[media] * pow(pow(S_w_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(S_w_e, (A / (1. - A) - 1.)) / (1. - A)
				       / (1. - (*gpu_def).S_wr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_gr[media]);
			}

			if (S_g_e <= (*gpu_def).S_g_range[1])
			{
				P_k_gn = ((*gpu_def).ag[0]) * S_g_e + ((*gpu_def).bg[0]);
				PkSn = (-1) * ((*gpu_def).ag[0]);
			}
			else if (S_g_e >= (*gpu_def).S_g_range[2])
			{
				P_k_gn = ((*gpu_def).ag[1]) * S_g_e + ((*gpu_def).bg[1]);
				PkSn = (-1) * ((*gpu_def).ag[1]);
			}
			else
			{
				P_k_gn = (*gpu_def).P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A);
				PkSn = (*gpu_def).P_d_gn[media] * pow(pow(1. - S_g_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S_g_e, A / (1. - A) - 1.) / (1. - A)
				       / (1. - (*gpu_def).S_wr[media] - (*gpu_def).S_nr[media] - (*gpu_def).S_gr[media]);
			}

			Sg = 1. - DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

			F1 = (*gpu_def).ro0_w * (1. + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm))
			     * DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.roS_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			F2 = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw - (*gpu_def).P_atm))
			     * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - DevArraysPtr.roS_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			F3 = (*gpu_def).ro0_g * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw + P_k_gn) / (*gpu_def).P_atm
			     * Sg - DevArraysPtr.roS_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];

			F1P = (*gpu_def).ro0_w * (*gpu_def).beta_w * DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			F2P = (*gpu_def).ro0_n * (*gpu_def).beta_n * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			F3P = (*gpu_def).ro0_g * Sg / (*gpu_def).P_atm;

			F1Sw = (*gpu_def).ro0_w * (1 + (*gpu_def).beta_w * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
			F2Sw = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * PkSw) * DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
			F2Sn = (*gpu_def).ro0_n * (1. + (*gpu_def).beta_n * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw - (*gpu_def).P_atm));
			F1Sn = 0;
			F3Sn = (-1) * (*gpu_def).ro0_g * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw + P_k_gn - Sg * PkSn) / (*gpu_def).P_atm;
			F3Sw = (-1) * (*gpu_def).ro0_g * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + P_k_nw + P_k_gn - Sg * (PkSn + PkSw)) / (*gpu_def).P_atm;

			// Вычисление дополнительных миноров матрицы частных производных
			a[0] = F2Sw * F3Sn - F3Sw * F2Sn;
			a[3] = F3Sw * F1Sn - F1Sw * F3Sn;
			a[6] = F1Sw * F2Sn - F2Sw * F1Sn;
			a[1] = F3P * F2Sn - F2P * F3Sn;
			a[4] = F1P * F3Sn - F3P * F1Sn;
			a[7] = F2P * F1Sn - F1P * F2Sn;
			a[2] = F2P * F3Sw - F3P * F2Sw;
			a[5] = F3P * F1Sw - F1P * F3Sw;
			a[8] = F1P * F2Sw - F2P * F1Sw;

			det = F1P * a[0] + F2P * a[3] + F3P * a[6];

			DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
			        - (1. / det) * (a[0] * F1 + a[3] * F2 + a[6] * F3);
			DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
			        - (1. / det) * (a[1] * F1 + a[4] * F2 + a[7] * F3);
			DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)]
			        - (1. / det) * (a[2] * F1 + a[5] * F2 + a[8] * F3);
		}


		device_test_positive(DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)], __FILE__, __LINE__);
	}
}

//Задание граничных условий отдельно для (Sw,Sg),Pn

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
__global__ void Border_S_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == ((*gpu_def).locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == ((*gpu_def).locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && (((*gpu_def).locNz) > 2))
	{
		k1 ++;
	}
	if ((k == ((*gpu_def).locNz) - 1) && (((*gpu_def).locNz) > 2))
	{
		k1 --;
	}

	if ((j != 0) || (((*gpu_def).source) <= 0))
	{
		DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_w[i1 + j1 * ((*gpu_def).locNx) + k1 * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
		DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.S_n[i1 + j1 * ((*gpu_def).locNx) + k1 * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
	}

	if ((j == 0) && (((*gpu_def).source) > 0))
	{
		DevArraysPtr.S_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).S_w_gr;
		DevArraysPtr.S_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).S_n_gr;
	}
}

__global__ void Border_P_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == ((*gpu_def).locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == ((*gpu_def).locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && (((*gpu_def).locNz) > 2))
	{
		k1 ++;
	}
	if ((k == ((*gpu_def).locNz) - 1) && (((*gpu_def).locNz) > 2))
	{
		k1 --;
	}

	if ((j != 0) && (j != ((*gpu_def).locNy) - 1))
	{
		DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i1 + j1 * ((*gpu_def).locNx) + k1 * ((*gpu_def).locNx) * ((*gpu_def).locNy)];
	}
	else if (j == 0)
	{
		DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).P_atm;
	}
	else
	{
		DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = DevArraysPtr.P_w[i1 + j1 * ((*gpu_def).locNx) + k1 * ((*gpu_def).locNx) * ((*gpu_def).locNy)] + cu_ro_eff_gdy(DevArraysPtr, i1, j1, k1);
	}
}

// Запись начальных условий для каждой точки сетки (независимо от остальных точек)
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

		int media = DevArraysPtr.media[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = 0;
		int j1 = (*gpu_def).locNy / 2;

		if (j < j1)
		{
			DevArraysPtr.S_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).S_w_gr + ((*gpu_def).S_w_init - (*gpu_def).S_w_gr) * j / j1;
			DevArraysPtr.S_n[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).S_n_gr + ((*gpu_def).S_n_init - (*gpu_def).S_n_gr) * j / j1;
		}
		else
		{
			DevArraysPtr.S_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).S_w_init;
			DevArraysPtr.S_n[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).S_n_init;
		}

		if (j == 0)
		{
			DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = (*gpu_def).P_atm;
		}
		else
		{
			DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] = DevArraysPtr.P_w[i + (j - 1) * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy] + cu_ro_eff_gdy(DevArraysPtr, i, j - 1, k);
		}

		DevArraysPtr.ro_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_w * (1. + ((*gpu_def).beta_w) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));

		///!!!! Не учитываются каппилярные силы! Или надо считать перед этим шагом P_w, P_g
		DevArraysPtr.ro_n[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_n * (1. + ((*gpu_def).beta_n) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));
		DevArraysPtr.ro_g[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] = (*gpu_def).ro0_g * (1. + ((*gpu_def).beta_g) * (DevArraysPtr.P_w[i + j * ((*gpu_def).locNx) + k * ((*gpu_def).locNx) * ((*gpu_def).locNy)] - (*gpu_def).P_atm));

		device_test_nan(DevArraysPtr.S_n[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.P_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.media[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.S_w[i + j * (*gpu_def).locNx + k * (*gpu_def).locNx * (*gpu_def).locNy], __FILE__, __LINE__);
	}
}


