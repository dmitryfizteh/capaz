#include "defines.h"

// явл€етс€ ли точка нагнетательной скважиной
int is_injection_well(int i, int j, int k, consts def)
{
#ifdef B_L
	if (((i == 1) && (j == 1)) || ((i == 0) && (j == 0)) || ((i == 1) && (j == 0)) || ((i == 0) && (j == 1)))
#endif
#ifdef THREE_PHASE
	if (j == (def.Ny) - 2)
#endif
#ifndef TWO_PHASE
		return 1;
	else
#endif
		return 0;
}

// явл€етс€ ли точка добывающей скважиной
int is_output_well(int i, int j, int k, consts def)
{
#ifdef B_L
	if (((i == def.Nx - 2) && (j == def.Ny - 2)) || ((i == def.Nx - 1) && (j == def.Ny - 1)) || ((i == def.Nx - 1) && (j == def.Ny - 2)) || ((i == def.Nx - 2) && (j == def.Ny - 1)))
#endif
#ifdef THREE_PHASE
	if (j == 1)
#endif
#ifndef TWO_PHASE
		return 1;
	else
#endif
		return 0;
}

// ”станавливает значени€ втекаемых/вытекаемых жидкостей q_i на скважинах
void wells_q(ptr_Arrays HostArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g, consts def)
{
#ifdef B_L
	// нагнетательна€ скважина
	if (is_injection_well(i, j, k, def))
	{
		*q_w = def.Q;
		*q_n = 0.;
		*q_g = 0.;
	}

	// добывающа€ скважина
	if (is_output_well(i, j, k, def))
	{
		*q_g = 0;
		double k_w=0., k_n=0.;
		assing_k(&k_w, &k_n, 1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

		double F_bl = (k_w / def.mu_w) / (k_w / def.mu_w + k_n / def.mu_n);
		*q_w = -1. * def.Q * F_bl;
		*q_n = -1. * def.Q * (1. - F_bl);
	}
#endif

#ifdef THREE_PHASE
/*
	double q = 0;
	if (is_output_well(i, j, k, def))
	{
	*q_w = 0.02;
	*q_g = 0.005;
	*q_n = 0;

	q = *q_w + *q_n + *q_g;

	*q_w = -q * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	*q_g = -q * (1 - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
	*q_n = -q * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	}

	if (is_injection_well(i, j, k, def))
	{
	*q_w = 0.02;
	*q_g = 0.005;
	*q_n = 0;
	}
*/
#endif
}

void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				if (is_active_point(i, j, k, def))
				{
					assign_P_Xi(HostArraysPtr, i, j, k, def);
					assign_ro(HostArraysPtr, i, j, k, def);
#ifdef THREE_PHASE
					assign_border_P(HostArraysPtr.P_w, HostArraysPtr.ro_w, i, j, k, def);
					assign_border_P(HostArraysPtr.P_n, HostArraysPtr.ro_n, i, j, k, def);
					assign_border_P(HostArraysPtr.P_g, HostArraysPtr.ro_g, i, j, k, def);
#endif
				}
			}
}

void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					assign_u(HostArraysPtr, i, j, k, def);
				}
}

void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
#ifdef NR
					assign_roS_nr(HostArraysPtr, t, i, j, k, def);
#else
					assign_roS(HostArraysPtr, t, i, j, k, def);
#endif
				}
}

void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					Newton(HostArraysPtr, i, j, k, def);
				}
}

void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
				if (is_active_point(i, j, k, def))
				{
					Border_S(HostArraysPtr, i, j, k, def);
					Border_P(HostArraysPtr, i, j, k, def);
				}
}

// ¬ычисление координаты точки, через которую будет вычисл€тьс€ значение на границе (i1, j1, k1)
void set_boundary_basic_coordinate(int i, int j, int k, int* i1, int* j1, int* k1, consts def)
{
	if (i == 0)
	{
		(*i1) ++;
	}
	if (i == (def.locNx) - 1)
	{
		(*i1) --;
	}
	if (j == 0)
	{
		(*j1) ++;
	}
	if (j == (def.locNy) - 1)
	{
		(*j1) --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		(*k1) ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		(*k1) --;
	}
}

void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
	HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
#ifdef THREE_PHASE
	HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_g * HostArraysPtr.P_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] / def.P_atm;
	test_positive(HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
#endif
	test_positive(HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_positive(HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

// –асчет центральной разности
double central_difference (double* ptr, char axis, consts def)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *(ptr-1) )/ (2. * (def.hx));	
		}
	case 'y':
		{
			return (*(ptr+def.locNx) - *(ptr-def.locNx) )/ (2. * (def.hy));
		}
	case 'z':
		{
			return (*(ptr + def.locNx * (def.locNy)) - *(ptr - def.locNx * (def.locNy)) )/ (2. * (def.hz));
		}
	default:
		{
			print_error("Axis of [central_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// –асчет левой разности
double left_difference (double* ptr, char axis, consts def)
{
	switch (axis)
	{
	case 'x':
		{
			return (*ptr - *(ptr-1) )/ def.hx;	
		}
	case 'y':
		{
			return (*ptr - *(ptr-def.locNx) )/ def.hy;
		}
	case 'z':
		{
			return (*ptr - *(ptr - def.locNx * (def.locNy)) )/ def.hz;
		}
	default:
		{
			print_error("Axis of [left_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// –асчет правой разности
double right_difference (double* ptr, char axis, consts def)
{
	switch (axis)
	{
	case 'x':
		{
			return (*(ptr+1) - *ptr )/ def.hx;	
		}
	case 'y':
		{
			return (*(ptr+def.locNx) - *ptr )/ def.hy;
		}
	case 'z':
		{
			return (*(ptr + def.locNx * (def.locNy)) - *ptr )/ def.hz;
		}
	default:
		{
			print_error("Axis of [right_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// –асчет скоростей в точке
void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
	if ((def.Nx) > 2)
	{
		if (i == 0)
		{
			HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * right_difference(HostArraysPtr.P_w+local, 'x', def);
			HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * right_difference(HostArraysPtr.P_n+local, 'x', def);
#ifdef THREE_PHASE
			HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * right_difference(HostArraysPtr.P_g+local, 'x', def);
#endif
		}
		else
		{
			if (i == (def.locNx) - 1)
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * left_difference(HostArraysPtr.P_w+local, 'x', def);
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * left_difference(HostArraysPtr.P_n+local, 'x', def);
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * left_difference(HostArraysPtr.P_g+local, 'x', def);
#endif
			}
			else
			{
				HostArraysPtr.ux_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'x', def);
				HostArraysPtr.ux_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'x', def);
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'x', def);
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.ux_w[local] = 0.;
		HostArraysPtr.ux_n[local] = 0.;
#ifdef THREE_PHASE
		HostArraysPtr.ux_g[local] = 0.;
#endif
	}

	if ((def.Ny) > 2)
	{
		if (j == 0)
		{
			HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (right_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local + 1] * (def.g_const));
			HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (right_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local + 1] * (def.g_const));
#ifdef THREE_PHASE
			HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (right_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local + 1] * (def.g_const));
#endif
		}
		else
		{
			if (j == (def.locNy) - 1)
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (left_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local - 1] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (left_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local - 1] * (def.g_const));
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (left_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local - 1] * (def.g_const));
#endif
			}
			else
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (central_difference (HostArraysPtr.P_w+local, 'y', def)	- HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (central_difference (HostArraysPtr.P_n+local, 'y', def)	- HostArraysPtr.ro_n[local] * (def.g_const));
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (central_difference (HostArraysPtr.P_g+local, 'y', def)	- HostArraysPtr.ro_g[local] * (def.g_const));
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uy_w[local] = 0.;
		HostArraysPtr.uy_n[local] = 0.;
#ifdef THREE_PHASE
		HostArraysPtr.uy_g[local] = 0.;
#endif
	}

	if ((def.Nz) > 2)
	{
		if (k == 0)
		{
			HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * right_difference (HostArraysPtr.P_w+local, 'z', def);
			HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * right_difference (HostArraysPtr.P_n+local, 'z', def);
#ifdef THREE_PHASE
			HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * right_difference (HostArraysPtr.P_g+local, 'z', def);
#endif
		}
		else
		{
			if (k == (def.locNz) - 1)
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * left_difference (HostArraysPtr.P_w+local, 'z', def);
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * left_difference (HostArraysPtr.P_n+local, 'z', def);
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * left_difference (HostArraysPtr.P_g+local, 'z', def);
#endif
			}
			else
			{
				HostArraysPtr.uz_w[local] = HostArraysPtr.Xi_w[local] * central_difference (HostArraysPtr.P_w+local, 'z', def);
				HostArraysPtr.uz_n[local] = HostArraysPtr.Xi_n[local] * central_difference (HostArraysPtr.P_n+local, 'z', def);
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[local] = HostArraysPtr.Xi_g[local] * central_difference (HostArraysPtr.P_g+local, 'z', def);
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uz_w[local] = 0.;
		HostArraysPtr.uz_n[local] = 0.;
#ifdef THREE_PHASE
		HostArraysPtr.uz_g[local] = 0.;
#endif
	}

	test_u(HostArraysPtr.ux_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.ux_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
#ifdef THREE_PHASE
	test_u(HostArraysPtr.ux_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
#endif
}

void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.Nz) < 2)))
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1 = 0, A2 = 0;

		double S_n = HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
#ifdef THREE_PHASE
		double divgrad3, Tx3, Ty3, Tz3, A3 = 0;

		double S_g = (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		double S_w = HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
#else
		double S_w = 1 - S_n;
#endif
		HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * S_w;
		HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * S_n;
#ifdef THREE_PHASE
		HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * S_g;
#endif

		if ((def.Nz) < 2)
		{
			divgrad1 = 0;
			divgrad2 = 0;
			Tz1 = 0;
			Tz2 = 0;
#ifdef THREE_PHASE
			divgrad3 = 0;
			Tz3 = 0;
#endif
		}
		else
		{
#ifdef THREE_PHASE
			divgrad1 = (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.)
			           * (HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			              + HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / ((def.hz) * (def.hz));
			divgrad3 = (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.)
			           * (HostArraysPtr.ro_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)])
			              - 2 * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			              * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
			              + HostArraysPtr.ro_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]
			              * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)])) / ((def.hz) * (def.hz));

			Tz3 = (HostArraysPtr.ro_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			       - HostArraysPtr.ro_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / (2. * (def.hz));
#else
			divgrad1 = (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.)
			           * (HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)])
			              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
			              + HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)])) / ((def.hz) * (def.hz));
#endif
			divgrad2 = (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.)
			           * (HostArraysPtr.ro_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              - 2 * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			              + HostArraysPtr.ro_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)])) / ((def.hz) * (def.hz));

			Tz1 = (HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			       - HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / (2. * (def.hz));
			Tz2 = (HostArraysPtr.ro_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			       - HostArraysPtr.ro_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / (2. * (def.hz));
		}

#ifdef THREE_PHASE
		divgrad1 += (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.) *
		            ((HostArraysPtr.ro_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hx) * (def.hx)) +
		             (HostArraysPtr.ro_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hy) * (def.hy)));

		divgrad3 += (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.) *
		            ((HostArraysPtr.ro_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hx) * (def.hx)) +
		             (HostArraysPtr.ro_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hy) * (def.hy)));

		Tx3 = (HostArraysPtr.ro_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hx));
		Ty3 = (HostArraysPtr.ro_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hy));
#else
		divgrad1 += (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.) *
		            ((HostArraysPtr.ro_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hx) * (def.hx)) +
		             (HostArraysPtr.ro_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              - 2 * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])
		              + HostArraysPtr.ro_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hy) * (def.hy)));
#endif
		divgrad2 += (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.) *
		            ((HostArraysPtr.ro_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		              - 2 * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		              + HostArraysPtr.ro_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hx) * (def.hx)) +
		             (HostArraysPtr.ro_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]
		              - 2 * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		              + HostArraysPtr.ro_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * (HostArraysPtr.S_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)])) / ((def.hy) * (def.hy)));

		Tx1 = (HostArraysPtr.ro_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hx));
		Ty1 = (HostArraysPtr.ro_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hy));
		Tx2 = (HostArraysPtr.ro_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ux_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hx));
		Ty2 = (HostArraysPtr.ro_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]
		       - HostArraysPtr.ro_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / (2. * (def.hy));

		test_arrowhead(Tx1 + Ty1 + Tz1, divgrad1, __FILE__, __LINE__);
		test_arrowhead(Tx2 + Ty2 + Tz2, divgrad2, __FILE__, __LINE__);

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// «начени€ q на скважинах
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

		if ((t < 2 * (def.dt)) || TWO_LAYERS)
		{
			A1 = HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + ((def.dt) / HostArraysPtr.m[local]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + ((def.dt) / HostArraysPtr.m[local]) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2);
#ifdef THREE_PHASE
			A3 = HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + ((def.dt) / HostArraysPtr.m[local]) * (divgrad3 - Tx3 - Ty3 - Tz3);
#endif
		}
		else
		{
			A1 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_w_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        + 4 * (def.tau) * HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
			A2 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_n_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        + 4 * (def.tau) * HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

			test_tau(HostArraysPtr.roS_w_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], A1, i + j * (def.locNx) + k * (def.locNx) * (def.locNy), def, __FILE__, __LINE__);
			test_tau(HostArraysPtr.roS_n_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], A2, i + j * (def.locNx) + k * (def.locNx) * (def.locNy), def, __FILE__, __LINE__);

#ifdef THREE_PHASE
			A3 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_g_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        + 4 * (def.tau) * HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
#endif
		}

		HostArraysPtr.roS_w_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_n_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A1;
		HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A2;

		test_positive(HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);

#ifdef THREE_PHASE
		HostArraysPtr.roS_g_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A3;
		test_positive(HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
#endif
		//		double delta_roS_w = HostArraysPtr.roS_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - HostArraysPtr.roS_w_old[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
		//		double delta_roS_n = HostArraysPtr.roS_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - HostArraysPtr.roS_n_old[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
		//		double delta_roS_g = HostArraysPtr.roS_g[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - HostArraysPtr.roS_g_old[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
	}
}


void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if(! HostArraysPtr.m[local])
			return;

		double q_w = 0.;
		double q_n = 0.;
		double q_g = 0.;

		// «начени€ q на скважинах
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

#ifdef THREE_PHASE
		HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
		        * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

		double Pg = HostArraysPtr.P_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double b1, b2, b3, A3 = 0;
#else
		HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
#endif
		HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double Pw = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double Pn = HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3, A1 = 0, A2 = 0;

		if ((def.Nz) < 2)
		{
			f3 = 0;
			g3 = 0;
#ifdef THREE_PHASE
			b3 = 0;
#endif
		}
		else
		{
			z2 = -(HostArraysPtr.P_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] - Pw) / def.hz;
			z1 = -(Pw - HostArraysPtr.P_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / def.hz;

			f3 = (((z2 + fabs(z2)) / 2.0 - (z1 - fabs(z1)) / 2.0) * (-1) * HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
			      (z1 + fabs(z1)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] +
			      (z2 - fabs(z2)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]) / def.hz;

			z2 = -(HostArraysPtr.P_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] - Pn) / def.hz;
			z1 = -(Pn - HostArraysPtr.P_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / def.hz;

			g3 = (((z2 + fabs(z2)) / 2.0 - (z1 - fabs(z1)) / 2.0) * (-1) * HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
			      (z1 + fabs(z1)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] +
			      (z2 - fabs(z2)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]) / def.hz;
#ifdef THREE_PHASE
			z2 = -(HostArraysPtr.P_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] - Pg) / def.hz;
			z1 = -(Pg - HostArraysPtr.P_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / def.hz;

			b3 = (((z2 + fabs(z2)) / 2.0 - (z1 - fabs(z1)) / 2.0) * (-1) * HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
			      (z1 + fabs(z1)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] +
			      (z2 - fabs(z2)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]) / def.hz;
#endif
		}

		x2 = -(HostArraysPtr.P_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pw) / def.hx;
		x1 = -(Pw - HostArraysPtr.P_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		y2 = -(HostArraysPtr.P_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pw) / def.hy + (def.g_const) * (HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		y1 = -(Pw - HostArraysPtr.P_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy + (def.g_const) * (HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

		f1 = (((x2 + fabs(x2)) / 2.0 - (x1 - fabs(x1)) / 2.0) * (-1) * HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (x1 + fabs(x1)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (x2 - fabs(x2)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		f2 = (((y2 + fabs(y2)) / 2.0 - (y1 - fabs(y1)) / 2.0) * (-1) * HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (y1 + fabs(y1)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (y2 - fabs(y2)) / 2.0 * (-1) * HostArraysPtr.Xi_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_w[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy;


		x2 = -(HostArraysPtr.P_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pn) / def.hx;
		x1 = -(Pn - HostArraysPtr.P_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		y2 = -(HostArraysPtr.P_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pn) / def.hy + (def.g_const) * (HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		y1 = -(Pn - HostArraysPtr.P_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy + (def.g_const) * (HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

		g1 = (((x2 + fabs(x2)) / 2.0 - (x1 - fabs(x1)) / 2.0) * (-1) * HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (x1 + fabs(x1)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (x2 - fabs(x2)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		g2 = (((y2 + fabs(y2)) / 2.0 - (y1 - fabs(y1)) / 2.0) * (-1) * HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (y1 + fabs(y1)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (y2 - fabs(y2)) / 2.0 * (-1) * HostArraysPtr.Xi_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_n[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy;

		A1 = HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - (def.dt / HostArraysPtr.m[local]) * (q_w + f1 + f2 + f3);
		A2 = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - (def.dt / HostArraysPtr.m[local]) * (q_n + g1 + g2 + g3);

		HostArraysPtr.roS_w_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_n_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A1;
		HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A2;

		test_positive(HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);

#ifdef THREE_PHASE
		x2 = -(HostArraysPtr.P_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pg) / def.hx;
		x1 = -(Pg - HostArraysPtr.P_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		y2 = -(HostArraysPtr.P_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pg) / def.hy + (def.g_const) * (HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		y1 = -(Pg - HostArraysPtr.P_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy + (def.g_const) * (HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

		b1 = (((x2 + fabs(x2)) / 2.0 - (x1 - fabs(x1)) / 2.0) * (-1) * HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (x1 + fabs(x1)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (x2 - fabs(x2)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;

		b2 = (((y2 + fabs(y2)) / 2.0 - (y1 - fabs(y1)) / 2.0) * (-1) * HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] -
		      (y1 + fabs(y1)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + (j - 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] +
		      (y2 - fabs(y2)) / 2.0 * (-1) * HostArraysPtr.Xi_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)] * HostArraysPtr.ro_g[i + (j + 1) * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hy;

		A3 = HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - (def.dt / HostArraysPtr.m[local]) * (q_g + b1 + b2 + b3);

		HostArraysPtr.roS_g_old[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = A3;

		test_positive(HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
#endif
	}
}

// ‘ункци€ загрузки данных в пам€ть хоста
void load_data_to_host(double *HostArrayPtr, double *DevArrayPtr, consts def)
{
}

// ‘ункци€ загрузки данных типа double в пам€ть ускорител€
void load_data_to_device(double *HostArrayPtr, double *DevArrayPtr, consts def)
{
}

// ‘ункци€ загрузки данных типа int в пам€ть ускорител€
void load_data_to_device_int(int *HostArrayPtr, int *DevArrayPtr, consts def)
{
}

// ¬ыделение пам€ти ускорител€ под массив точек расчетной области
void device_memory_allocation(ptr_Arrays *ArraysPtr, double **DevBuffer, consts def)
{
}

// ќсвобожение пам€ти ускорител€ из под массива точек расчетной области
void device_memory_free(ptr_Arrays ptDev, double *DevBuffer)
{
}

// »нициализаци€ ускорител€
void device_initialization(consts *def)
{
}

// ‘инализаци€ ускорител€
void device__finalization(void)
{
}

// «агрузка на хост данных дл€ обмена на границе
void load_exchange_data(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArrayPtr[1 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			HostBuffer[j + (def.locNy)*k + (def.locNy) * (def.locNz)] = HostArrayPtr[(def.locNx) - 2 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];

			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
			test_nan(HostBuffer[j + (def.locNy)*k + (def.locNy) * (def.locNz)], __FILE__, __LINE__);
		}

	/*for(int j=0;j<(def.locNy);j++)
		for(int k=0;k<(def.locNz);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(def.locNy)*k]);*/
}

// «агрузка на device данных обмена на границе
void save_exchange_data(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	//printf("\nSave\n");
	/*for(int j=0;j<(def.locNy);j++)
		for(int k=0;k<(def.locNz);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(def.locNy)*k]);*/

	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			if ((def.rank) != def.sizex - 1)
			{
				HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k + (def.locNy) * (def.locNz)];
				test_nan(HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
			}
			if ((def.rank) != 0)
			{
				HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
				test_nan(HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
			}
		}
}

