#include "../defines.h"
#include "two-phase.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int media = 0;
	double S_e = (1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media]);
	//if (S_e<0)
	//	S_e=0;
	double k_w = pow(S_e, (2. + 3. * (def.lambda[media])) / def.lambda[media]);
	double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));
	double P_k = def.P_d[media] * pow((1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media]), -1. / def.lambda[media]);

	HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k;
	HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = -1. * (def.K[media]) * k_n / def.mu_n;

	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int media = 0;
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_e = (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1 - def.S_wr[media]);
			P_k = def.P_d[media] * pow(S_e, (-1.) / def.lambda[media]);
			AAA = pow(S_e, (((-1.) / def.lambda[media]) - 1.));
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) - HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k - def.P_atm)) * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			PkS = AAA * (def.P_d[media]) / (def.lambda[media] * (1 - def.S_wr[media]));
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k - def.P_atm + (HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * PkS)));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - (1 / det) * (F1P * F2 - F2P * F1);
		}

		test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);

	}
}

// «адание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1, def);

	if ((j != 0) || ((def.source) <= 0))
	{
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}

	if ((j == 0) && ((def.source) > 0))
	{
		int I = local_to_global(i, 'x', def);
		if ((I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
		{
			HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.S_n_gr;
		}
		else
			//HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = 0;
		{
			HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
		}
	}

	test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1, def);

	if ((j != 0) && (j != (def.locNy) - 1))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}
	else if (j == 0)
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	}
	else
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] + ro_eff_gdy(HostArraysPtr, i1, j1, k1, def);
	}

	test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					// ѕреобразование локальных координат процессора к глобальным
					int I = local_to_global(i, 'x', def);
					HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy]=def.porosity[0];

					// ≈сли точка на верхней границе, не далее (def.source) точек от центра,
					// то в ней начальна€ насыщенность. »наче, нулева€
					if ((j == 0) && (I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
					{
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = 0;
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

					///!!!! Ќе учитываютс€ капилл€рные силы! »ли надо считать перед этим шагом P_n
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					/*
					if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=(def.NX)/2.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=4.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=2./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=(-1.)*HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]/4.+2./5.*def.locNy*(def.h2)))
							HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;

					if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=2.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=4./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=3./5.*def.locNy*(def.h2)))
							HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;
							*/

					/*
					if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=2.*(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=3.*(def.NX)/5.*(def.h1)))
						if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=1./10.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=3./10.*def.locNy*(def.h2)))
							HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;
					*/

					test_nan(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}


