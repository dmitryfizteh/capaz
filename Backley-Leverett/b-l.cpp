#include "../defines.h"
#include "b-l.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	double k_w, k_n;
	double S = 1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

	double S_wc = 0.2;
	double S_wi = 0.2;
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

	HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = -1. * (HostArraysPtr.K[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = -1. * (HostArraysPtr.K[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) * k_n / def.mu_n;

	//	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		double A1 = HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double A2 = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double a = def.beta_w * (def.beta_n);
		double b = def.beta_w + def.beta_n - A2 * (def.beta_w) / (def.ro0_n) - A1 * (def.beta_n) / (def.ro0_w);
		double c = 1 - A2 / def.ro0_n  - A1 / def.ro0_w;
		double D = b * b - 4 * a * c;
		double P1 = def.P_atm + (-1 * b + sqrt(D)) / (2 * a);
		double P2 = def.P_atm + (-1 * b - sqrt(D)) / (2 * a);

		if (P1 < 0)
		{
			HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = P2;
		}
		else
		{
			HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = P1;
		}

		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] / (def.ro0_n * (1 + def.beta_n * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm)));
		/*
		if ( HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] < 0.2)
			HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = 0.2;

		if ( HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] > 0.8)
			HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = 0.8;
		*/

		test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	}

}

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (def.locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (def.locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		k1 --;
	}

	HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];


	// В центре резервуара находится нагнетательная скважина
	if ((i == def.Nx / 2) && (j == def.Ny - 3) && (k == def.Nz / 2))
	{
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = INJECTION_WELL_Sn;
	}

	// В центре резервуара находится добывающая скважина
	if ((i == def.Nx - 3) && (j == 3) && (k == def.Nz - 3))
	{
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = OUTPUT_WELL_Sn;
	}


	test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

// Давление на нагнетательной скважине
double Injection_well_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	// 10000psi in Pa
	return INJECTION_WELL_Pw;//def.P_atm + 100000; //68947572.9;
}

// Давление на добывающей скважине
double Production_well_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	// 4000psi in Pa
	return OUTPUT_WELL_Pw;//def.P_atm;//27579029.16;
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (def.locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (def.locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		k1 --;
	}


	if ((j != 0) && (j != (def.locNy) - 1))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}
	//else if(j == 0)
	//	HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	else
	{
		double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		                  + def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[(HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy])]) * (def.g_const) * (def.hy);

		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] + ro_g_dy;//HostArraysPtr.ro_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] * (def.g_const) * (def.hy);
	}


	// В центре резервуара находится нагнетающая скважина
	if ((i == def.Nx / 2) && (j == def.Ny - 3) && (k == def.Nz / 2))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Injection_well_P(HostArraysPtr, i, j, k, def);
	}

	// В центре резервуара находится добывающая скважина
	if ((i == def.Nx - 3) && (j == 3) && (k == def.Nz - 3))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Production_well_P(HostArraysPtr, i, j, k, def);
	}


	test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
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

					HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy] = 0;
					HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = BACKGROUND_Sn;
					//HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy] =0.3 + 0.3 * j / def.Ny;

					double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
					                  + def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[(HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy])]) * (def.g_const) * (def.hy);

					// 6000 pound per square inch = 41 368 543.8 Паскаля
					if (j == 0)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = BACKGROUND_Pw;//def.P_atm + 50000;    // 50368543.8;
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_g_dy;
					}

					
					// В центре резервуара находится нагнетающая скважина
					if ((i == def.Nx / 2) && (j == def.Ny - 3) && (k == def.Nz / 2))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Injection_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}

					// В центре резервуара находится добывающая скважина
					if ((i == def.Nx - 3) && (j == 3) && (k == def.Nz - 3))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = Production_well_P(HostArraysPtr, i, j, k, def);
						//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0.5;
					}
					

					HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					test_nan(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}

