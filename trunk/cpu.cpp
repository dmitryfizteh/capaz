#include "defines.h"

// �������� �� ����� �������������� ���������
int is_injection_well(int i, int j, int k, consts def)
{
#ifdef B_L
	if (((i == 1) && (j == 1)) || ((i == 0) && (j == 0)) || ((i == 1) && (j == 0)) || ((i == 0) && (j == 1)))
#endif
#ifdef THREE_PHASE
//	if ((j == 1) && (i > 0) && (i < (def.locNx) / 3 - 1) && (((def.locNz) < 2) || (k > 0) && (k < (def.locNz) / 3 - 1)))
if ((i <= (def.Nx) / 4 + 1 && i >= (def.Nx) / 4 - 1 && j <= (def.Ny) / 4 + 1 && j >= (def.Ny) / 4 - 1) 
	|| (i <= 3 * (def.Nx) / 4 + 1 && i >= 3 * (def.Nx) / 4 - 1 && j <= (def.Ny) / 4 + 1 && j >= (def.Ny) / 4 - 1)
	|| (i <= (def.Nx) / 4 + 1 && i >= (def.Nx) / 4 - 1 && j <= 3 * (def.Ny) / 4 + 1 && j >= 3 * (def.Ny) / 4 - 1)
	|| (i <= 3 * (def.Nx) / 4 + 1 && i >= 3 * (def.Nx) / 4 - 1 && j <= 3 * (def.Ny) / 4 + 1 && j >= 3 * (def.Ny) / 4 - 1))
#endif
#ifndef TWO_PHASE
		return 1;
	else
#endif
		return 0;
}

// �������� �� ����� ���������� ���������
int is_output_well(int i, int j, int k, consts def)
{
#ifdef B_L
	if (((i == def.Nx - 2) && (j == def.Ny - 2)) || ((i == def.Nx - 1) && (j == def.Ny - 1)) || ((i == def.Nx - 1) && (j == def.Ny - 2)) || ((i == def.Nx - 2) && (j == def.Ny - 1)))
#endif
#ifdef THREE_PHASE
//	if ((j == 1) && (i > 0) && (i < (def.locNx) / 3 - 1) && (((def.locNz) < 2) || (k > 0) && (k < (def.locNz) / 3 - 1)))
if ((i <= (def.Nx) / 4 + 1 && i >= (def.Nx) / 4 - 1 && j <= (def.Ny) / 2 + 1 && j >= (def.Ny) / 2 - 1) 
	|| (i <= (def.Nx) / 2 + 1 && i >= (def.Nx) / 2 - 1 && j <= (def.Ny) / 4 + 1 && j >= (def.Ny) / 4 - 1)
	|| (i <= (def.Nx) / 2 + 1 && i >= (def.Nx) / 2 - 1 && j <= 3 * (def.Ny) / 4 + 1 && j >= 3 * (def.Ny) / 4 - 1)
	|| (i <= 3 * (def.Nx) / 4 + 1 && i >= 3 * (def.Nx) / 4 - 1 && j <= (def.Ny) / 2 + 1 && j >= (def.Ny) / 2 - 1))
#endif
#ifndef TWO_PHASE
		return 1;
	else
#endif
		return 0;
}

// ������������� �������� ���������/���������� ��������� q_i �� ���������
void wells_q(ptr_Arrays HostArraysPtr, int i, int j, int k, double* q_w, double* q_n, double* q_g, consts def)
{
#ifdef B_L
	// �������������� ��������
	if (is_injection_well(i, j, k, def))
	{
		*q_w = def.Q;
		*q_n = 0.;
		*q_g = 0.;
	}

	// ���������� ��������
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

	double q = 0.;
	*q_w = 0.0;
	*q_g = 0.0;
	*q_n = 0.0;

/*	if (is_injection_well(i, j, k, def))
	{
		*q_w = 0.01;
		*q_g = 0.005;
		*q_n = 0.0;
	}
	if (is_output_well(i, j, k, def))
	{
		q = 0.015;

		*q_w = -q * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		*q_g = -q * (1 - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
		*q_n = -q * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	}*/
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

// ���������� ���������� �����, ����� ������� ����� ����������� �������� �� ������� (i1, j1, k1)
int set_boundary_basic_coordinate(int i, int j, int k, consts def)
{
	int i1, j1, k1;

	i1 = i;
	j1 = j;
	k1 = k;

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

	return (i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy));
}

void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

	HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
	HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[local] - def.P_atm));
#ifdef THREE_PHASE
	HostArraysPtr.ro_g[local] = def.ro0_g * HostArraysPtr.P_g[local] / def.P_atm;
	test_ro(HostArraysPtr.ro_g[local], __FILE__, __LINE__);
#endif
	test_ro(HostArraysPtr.ro_w[local], __FILE__, __LINE__);
	test_ro(HostArraysPtr.ro_n[local], __FILE__, __LINE__);
}

// ������ ����������� ��������
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

// ������ ������������ ��������
double directed_difference (double x1, double x2, double* Xi, double* ro, char axis, consts def)
{
	switch (axis)
	{
	case 'x':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-1)) * (*(ro-1)) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+1)) * (*(ro+1))) / def.hx;
		}
	case 'y':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-def.locNx)) * (*(ro-def.locNx)) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+def.locNx)) * (*(ro+def.locNx))) / def.hy;
		}
	case 'z':
		{
			return (((x2 + fabs(x2)) / 2. - (x1 - fabs(x1)) / 2.) * (-1.) * (*Xi) * (*ro) -
		      (x1 + fabs(x1)) / 2. * (-1.) * (*(Xi-def.locNx * (def.locNy))) * (*(ro-def.locNx * (def.locNy))) +
		      (x2 - fabs(x2)) / 2. * (-1.) * (*(Xi+def.locNx * (def.locNy))) * (*(ro+def.locNx * (def.locNy)))) / def.hz;
		}
	default:
		{
			print_error("Axis of [directed_difference] conversation is empty", __FILE__, __LINE__);
			return -1;
		}
	}
}

// ������ ����� ��������
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

// ������ ������ ��������
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

// ������ ��������� � �����
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
			HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (right_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local] * (def.g_const));
			HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (right_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local] * (def.g_const));
#ifdef THREE_PHASE
			HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (right_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local] * (def.g_const));
#endif
		}
		else
		{
			if (j == (def.locNy) - 1)
			{
				HostArraysPtr.uy_w[local] = HostArraysPtr.Xi_w[local] * (left_difference (HostArraysPtr.P_w+local, 'y', def) - HostArraysPtr.ro_w[local] * (def.g_const));
				HostArraysPtr.uy_n[local] = HostArraysPtr.Xi_n[local] * (left_difference (HostArraysPtr.P_n+local, 'y', def) - HostArraysPtr.ro_n[local] * (def.g_const));
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[local] = HostArraysPtr.Xi_g[local] * (left_difference (HostArraysPtr.P_g+local, 'y', def) - HostArraysPtr.ro_g[local] * (def.g_const));
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

	test_u(HostArraysPtr.ux_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.ux_n[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_n[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_w[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_n[local], __FILE__, __LINE__);
#ifdef THREE_PHASE
	test_u(HostArraysPtr.ux_g[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uy_g[local], __FILE__, __LINE__);
	test_u(HostArraysPtr.uz_g[local], __FILE__, __LINE__);
#endif
}

void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.Nz) < 2)))
	{
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);
		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1 = 0, A2 = 0;

#ifdef THREE_PHASE
		double divgrad3, Tx3, Ty3, Tz3, A3 = 0;

		double S_g = (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local]);
		double S_w = HostArraysPtr.S_w[local];
#else
		double S_w = 1. - HostArraysPtr.S_n[local];
#endif
		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * S_w;
		HostArraysPtr.roS_n[local] = HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local];
#ifdef THREE_PHASE
		HostArraysPtr.roS_g[local] = HostArraysPtr.ro_g[local] * S_g;
#endif

		if ((def.Nz) < 2)
		{
			divgrad1 = 0.;
			divgrad2 = 0.;
			Tz1 = 0.;
			Tz2 = 0.;
#ifdef THREE_PHASE
			divgrad3 = 0.;
			Tz3 = 0.;
#endif
		}
		else
		{
#ifdef THREE_PHASE
			divgrad1 = (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.)
			           * (HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              - 2 * HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local]
			              + HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / ((def.hz) * (def.hz));
			divgrad3 = (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.)
			           * (HostArraysPtr.ro_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)])
			              - 2 * HostArraysPtr.ro_g[local] * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
			              + HostArraysPtr.ro_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]
			              * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)])) / ((def.hz) * (def.hz));

			Tz3 = (HostArraysPtr.ro_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_g[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			       - HostArraysPtr.ro_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.uz_g[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)]) / (2. * (def.hz));
#else
			divgrad1 = (HostArraysPtr.m[local] * (def.l) * (def.c_w) / 2.)
			           * (HostArraysPtr.ro_w[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)])
			              - 2 * HostArraysPtr.ro_w[local] * (1. - HostArraysPtr.S_n[local])
			              + HostArraysPtr.ro_w[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)] * (1. - HostArraysPtr.S_n[i + j * (def.locNx) + (k - 1) * (def.locNx) * (def.locNy)])) / ((def.hz) * (def.hz));
#endif
			divgrad2 = (HostArraysPtr.m[local] * (def.l) * (def.c_n) / 2.)
			           * (HostArraysPtr.ro_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + (k + 1) * (def.locNx) * (def.locNy)]
			              - 2 * HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local]
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

		// �������� q �� ���������
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

		if ((t < 2 * (def.dt)) || TWO_LAYERS)
		{
			A1 = HostArraysPtr.roS_w[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2);
#ifdef THREE_PHASE
			A3 = HostArraysPtr.roS_g[local] + ((def.dt) / HostArraysPtr.m[local]) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3);
#endif
		}
		else
		{
			A1 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_w + divgrad1 - Tx1 - Ty1 - Tz1)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_w_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
			A2 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_n + divgrad2 - Tx2 - Ty2 - Tz2)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_n_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);

			test_tau(HostArraysPtr.roS_w_old[local], HostArraysPtr.roS_w[local], A1, local, def, __FILE__, __LINE__);
			test_tau(HostArraysPtr.roS_n_old[local], HostArraysPtr.roS_n[local], A2, local, def, __FILE__, __LINE__);

#ifdef THREE_PHASE
			A3 = (1. / ((HostArraysPtr.m[local]) * (def.dt) + 2. * (def.tau))) * (2. * (def.dt) * (def.dt) * (q_g + divgrad3 - Tx3 - Ty3 - Tz3)
			        + ((HostArraysPtr.m[local]) * (def.dt) - 2. * (def.tau)) * HostArraysPtr.roS_g_old[local]
			        + 4. * (def.tau) * HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]);
#endif
		}

		HostArraysPtr.roS_w_old[local] = HostArraysPtr.roS_w[local];
		HostArraysPtr.roS_n_old[local] = HostArraysPtr.roS_n[local];
		HostArraysPtr.roS_w[local] = A1;
		HostArraysPtr.roS_n[local] = A2;

		test_positive(HostArraysPtr.roS_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[local], __FILE__, __LINE__);

#ifdef THREE_PHASE
		HostArraysPtr.roS_g_old[local] = HostArraysPtr.roS_g[local];
		HostArraysPtr.roS_g[local] = A3;
		test_positive(HostArraysPtr.roS_g[local], __FILE__, __LINE__);
#endif
		//		double delta_roS_w = HostArraysPtr.roS_w[local] - HostArraysPtr.roS_w_old[local];
		//		double delta_roS_n = HostArraysPtr.roS_n[local] - HostArraysPtr.roS_n_old[local];
		//		double delta_roS_g = HostArraysPtr.roS_g[local] - HostArraysPtr.roS_g_old[local];
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

		// �������� q �� ���������
		wells_q(HostArraysPtr, i, j, k, &q_w, &q_n, &q_g, def);

#ifdef THREE_PHASE
		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * HostArraysPtr.S_w[local];
		HostArraysPtr.roS_g[local] = HostArraysPtr.ro_g[local]
		        * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local]);

		double Pg = HostArraysPtr.P_g[local];
		double fx_g, fy_g, fz_g, A3 = 0.;
#else
		HostArraysPtr.roS_w[local] = HostArraysPtr.ro_w[local] * (1. - HostArraysPtr.S_n[local]);
#endif
		HostArraysPtr.roS_n[local] = HostArraysPtr.ro_n[local] * HostArraysPtr.S_n[local];
		double Pw = HostArraysPtr.P_w[local];
		double Pn = HostArraysPtr.P_n[local];

		double x1, x2, y1, y2, z1, z2, fx_w, fy_w, fz_w, fx_n, fy_n, fz_n, A1 = 0., A2 = 0.;

		if ((def.Nz) < 2)
		{
			fz_w = 0.;
			fz_n = 0.;
#ifdef THREE_PHASE
			fz_g = 0.;
#endif
		}
		else
		{
			z2 = -1. * right_difference (HostArraysPtr.P_w+local, 'z', def);
			z1 = -1. * left_difference (HostArraysPtr.P_w+local, 'z', def);
			fz_w = directed_difference (z1, z2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'z', def);

			z2 = -1. * right_difference (HostArraysPtr.P_n+local, 'z', def); 
			z1 = -1. * left_difference (HostArraysPtr.P_n+local, 'z', def); 
			fz_n = directed_difference (z1, z2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'z', def);
#ifdef THREE_PHASE
			z2 = -1. * right_difference (HostArraysPtr.P_g+local, 'z', def); 
			z1 = -1. * left_difference (HostArraysPtr.P_g+local, 'z', def); 
			fz_g = directed_difference (z1, z2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'z', def);
#endif
		}

		x2 = -1. * right_difference (HostArraysPtr.P_w+local, 'x', def); //-(HostArraysPtr.P_w[i + 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - Pw) / def.hx;
		x1 = -1. * left_difference (HostArraysPtr.P_w+local, 'x', def); //-(Pw - HostArraysPtr.P_w[i - 1 + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) / def.hx;
		y2 = -1. * right_difference (HostArraysPtr.P_w+local, 'y', def) + def.g_const * (HostArraysPtr.ro_w[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_w+local, 'y', def) + def.g_const * (HostArraysPtr.ro_w[local]);
		fx_w = directed_difference (x1, x2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'x', def);
		fy_w = directed_difference (y1, y2, HostArraysPtr.Xi_w+local, HostArraysPtr.ro_w+local, 'y', def);

		x2 = -1. * right_difference (HostArraysPtr.P_n+local, 'x', def); 
		x1 = -1. * left_difference (HostArraysPtr.P_n+local, 'x', def); 
		y2 = -1. * right_difference (HostArraysPtr.P_n+local, 'y', def) + def.g_const * (HostArraysPtr.ro_n[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_n+local, 'y', def) + def.g_const * (HostArraysPtr.ro_n[local]);
		fx_n = directed_difference (x1, x2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'x', def);
		fy_n = directed_difference (y1, y2, HostArraysPtr.Xi_n+local, HostArraysPtr.ro_n+local, 'y', def);

		A1 = HostArraysPtr.roS_w[local] - (def.dt / HostArraysPtr.m[local]) * (-q_w + fx_w + fy_w + fz_w);
		A2 = HostArraysPtr.roS_n[local] - (def.dt / HostArraysPtr.m[local]) * (-q_n + fx_n + fy_n + fz_n);

		HostArraysPtr.roS_w_old[local] = HostArraysPtr.roS_w[local];
		HostArraysPtr.roS_n_old[local] = HostArraysPtr.roS_n[local];
		HostArraysPtr.roS_w[local] = A1;
		HostArraysPtr.roS_n[local] = A2;

		test_positive(HostArraysPtr.roS_w[local], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[local], __FILE__, __LINE__);

#ifdef THREE_PHASE
		x2 = -1. * right_difference (HostArraysPtr.P_g+local, 'x', def); 
		x1 = -1. * left_difference (HostArraysPtr.P_g+local, 'x', def); 
		y2 = -1. * right_difference (HostArraysPtr.P_g+local, 'y', def) + def.g_const * (HostArraysPtr.ro_g[local]);
		y1 = -1. * left_difference (HostArraysPtr.P_g+local, 'y', def) + def.g_const * (HostArraysPtr.ro_g[local]);
		fx_g = directed_difference (x1, x2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'x', def);
		fy_g = directed_difference (y1, y2, HostArraysPtr.Xi_g+local, HostArraysPtr.ro_g+local, 'y', def);

		A3 = HostArraysPtr.roS_g[local] - (def.dt / HostArraysPtr.m[local]) * (q_g + fx_g + fy_g + fz_g);

		HostArraysPtr.roS_g_old[local] = HostArraysPtr.roS_g[local];
		HostArraysPtr.roS_g[local] = A3;

		test_positive(HostArraysPtr.roS_g[local], __FILE__, __LINE__);
#endif
	}
}

// ������� �������� ������ � ������ �����
void load_data_to_host(double *HostArrayPtr, double *DevArrayPtr, consts def)
{
}

// ������� �������� ������ ���� double � ������ ����������
void load_data_to_device(double *HostArrayPtr, double *DevArrayPtr, consts def)
{
}

// ������� �������� ������ ���� int � ������ ����������
void load_data_to_device_int(int *HostArrayPtr, int *DevArrayPtr, consts def)
{
}

// ��������� ������ ���������� ��� ������ ����� ��������� �������
void device_memory_allocation(ptr_Arrays *ArraysPtr, double **DevBuffer, consts def)
{
}

// ����������� ������ ���������� �� ��� ������� ����� ��������� �������
void device_memory_free(ptr_Arrays ptDev, double *DevBuffer)
{
}

// ������������� ����������
void device_initialization(consts *def)
{
}

// ����������� ����������
void device_finalization(void)
{
}

// �������� � ����� ������ ��� ������ �� �������. ��� ������� �� ����������� ���� �������. ����������� - ��� ��� ��������� � ����/�����.
void load_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArrayPtr[1 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[j + (def.locNy)*k] = HostArrayPtr[(def.locNx) - 2 + (def.locNx) * j + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[j + (def.locNy)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArrayPtr[i + (def.locNx) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostBuffer[i + (def.locNx)*k] = HostArrayPtr[i + (def.locNx) * ((def.locNy) - 2) + (def.locNx) * (def.locNy) * k];
			test_nan(HostBuffer[i + (def.locNx)*k], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

void load_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostBuffer[i + (def.locNx)*j] = HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 2)];
			test_nan(HostBuffer[i + (def.locNx)*j], __FILE__, __LINE__);
		}
}

// �������� �� ������ ������ ������ �� �������. ��� ������� �� ����������� ���� �������. ����������� - ��� ��� ��������� � ����/�����.
void save_exchange_data_part_xl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArrayPtr[(def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_xr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int j = 0; j < (def.locNy); j++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k] = HostBuffer[j + (def.locNy) * k];
			test_nan(HostArrayPtr[(def.locNx) - 1 + (def.locNx)*j + (def.locNx) * (def.locNy)*k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[i + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArrayPtr[i + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_yr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int k = 0; k < (def.locNz); k++)
		{
			HostArrayPtr[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k] = HostBuffer[i + (def.locNx)*k];
			test_nan(HostArrayPtr[i + (def.locNx) * ((def.locNy) - 1) + (def.locNx) * (def.locNy) * k], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zl(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArrayPtr[i + (def.locNx) * j] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArrayPtr[i + (def.locNx) * j], __FILE__, __LINE__);
		}
}

void save_exchange_data_part_zr(double *HostArrayPtr, double *DevArrayPtr, double *HostBuffer, double *DevBuffer, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
		{
			HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)] = HostBuffer[i + (def.locNx)*j];
			test_nan(HostArrayPtr[i + (def.locNx) * j + (def.locNx) * (def.locNy) * ((def.locNz) - 1)], __FILE__, __LINE__);
		}
}

