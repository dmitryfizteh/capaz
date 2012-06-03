#include "../defines.h"
#include "b-l.h"

// Расчет относительных проницаемостей в точке
void assing_k(double* k_w, double* k_n, double S_w)
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

	test_S(*k_n, __FILE__, __LINE__);
	test_S(*k_w, __FILE__, __LINE__);
}

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	double k_w=0., k_n=0.;

	assing_k(&k_w, &k_n, 1. - HostArraysPtr.S_n[local]);

	HostArraysPtr.P_n[local] = HostArraysPtr.P_w[local];
	HostArraysPtr.Xi_w[local] = -1. * (HostArraysPtr.K[local]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[local] = -1. * (HostArraysPtr.K[local]) * k_n / def.mu_n;

	//	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[local], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[local], __FILE__, __LINE__);
}

void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int media = 0;
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		double F1, F2, F1P, F2P, F1S, F2S, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm)) * (1. - HostArraysPtr.S_n[local]) - HostArraysPtr.roS_w[local];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[local] - def.P_atm)) * HostArraysPtr.S_n[local] - HostArraysPtr.roS_n[local];
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[local]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[local];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[local]  - def.P_atm));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local] - (1. / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local] - (1. / det) * (F1P * F2 - F2P * F1);
		}

		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);

	}
}

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	//if ((i == 0) || (i == (def.locNx) - 1) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	//{
		int local1=set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		HostArraysPtr.S_n[local] = HostArraysPtr.S_n[local1];

		/*
		// В центре резервуара находится нагнетательная скважина
		if (is_injection_well(i, j, k, def))
		{
			HostArraysPtr.S_n[local] = INJECTION_WELL_Sn;
		}

		// В центре резервуара находится добывающая скважина
		if (is_output_well(i, j, k, def))
		{
			HostArraysPtr.S_n[local] = OUTPUT_WELL_Sn;
		}
		*/

		test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
	//}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	//if ((i == 0) || (i == (def.locNx) - 1) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	//{
		int local1=set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((j != 0) && (j != (def.locNy) - 1))
		{
			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1];
		}
		//else if(j == 0)
		//	HostArraysPtr.P_w[local] = def.P_atm;
		else
		{
			double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[local]
							  + def.ro0_w * (1 - HostArraysPtr.S_n[local])) * (HostArraysPtr.m[i + j * def.locNx + k * def.locNx * def.locNy]) * (def.g_const) * (def.hy);

			HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local1] + ro_g_dy;//HostArraysPtr.ro_w[local1] * (def.g_const) * (def.hy);
		}

	
		// В центре резервуара находится нагнетающая скважина
		if (is_injection_well(i, j, k, def))
		//if (((i == 0) && (j == 0)) || ((i == 1) && (j == 0)) || ((i == 0) && (j == 1)))
		{
			HostArraysPtr.P_w[local] = def.InjWell_Pw;
		}

		// В центре резервуара находится добывающая скважина
		if (is_output_well(i, j, k, def))
		//if (((i == def.Nx - 1) && (j == def.Ny - 1)) || ((i == def.Nx - 1) && (j == def.Ny - 2)) || ((i == def.Nx - 2) && (j == def.Ny - 1)))
		{
			HostArraysPtr.P_w[local] = def.OutWell_Pw;
		}
	

		test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
	//}
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
					int local = i + j * def.locNx + k * def.locNx * def.locNy;
					HostArraysPtr.m[local]=def.porosity[0];
					HostArraysPtr.S_n[local] = def.Background_Sn;

					double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[local]
					                  + def.ro0_w * (1 - HostArraysPtr.S_n[local])) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);

					if (j == 0)
					{
						HostArraysPtr.P_w[local] = def.Background_Pw; //def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[local] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_g_dy;
					}

					/*
					// нагнетательная скважина
					if (is_injection_well(i, j, k, def))
					{
						HostArraysPtr.P_w[local] = Injection_well_P(HostArraysPtr, i, j, k, def);
					}

					// добывающая скважина
					if (is_output_well(i, j, k, def))
					{
						HostArraysPtr.P_w[local] = Production_well_P(HostArraysPtr, i, j, k, def);
					}
					*/

					HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
					HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[local] - def.P_atm));

					test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.m[local], __FILE__, __LINE__);
				}
}
