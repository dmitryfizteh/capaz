#include "../gpu.h"
#include "three-phase.h"

// ������� ���������� ����������� �������� �������������
__device__ double device_assign_S_w_e(ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	int media = 0;
	return (DevArraysPtr.S_w[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - gpu_def->S_wr[media]) / (1. - gpu_def->S_wr[media] - gpu_def->S_nr[media] - gpu_def->S_gr[media]);
}

__device__ double device_assign_S_n_e(ptr_Arrays DevArraysPtr, int i, int j, int k)
{
	int media = 0;
	return (DevArraysPtr.S_n[i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy)] - gpu_def->S_nr[media]) / (1. - gpu_def->S_wr[media] - gpu_def->S_nr[media] - gpu_def->S_gr[media]);
}

// ���������� ����������� ��������

// ������� ������� ���������� ����� ������������� ��� ���������� ����������� ��������
__device__ void device_assign_S_w_range(double* S1, double* S2)
{
	(*S1) = 0.1;
	(*S2) = 0.99;
}

__device__ void device_assign_S_g_range(double* S1, double* S2)
{
	(*S1) = 0.005;
	(*S2) = 0.95;
}

// ������� ���. �������� � �� ����������� ��� ����������� ����� ���������
__device__ double device_P_k_nw(double S)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	return gpu_def->P_d_nw[media] * pow((pow(S, A / (1. - A)) - 1.), 1. / A);
}

__device__ double device_P_k_gn(double S)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	return gpu_def->P_d_gn[media] * pow(pow((1. - S), A / (1. - A)) - 1., 1. / A);
}

__device__ double device_P_k_nw_S(double S)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	return gpu_def->P_d_nw[media] * pow(pow(S, A / (1. - A)) - 1., 1. / A - 1.) * pow(S, (A / (1. - A) - 1.)) / (1. - A)
		/ (1. - gpu_def->S_wr[media] - gpu_def->S_nr[media] - gpu_def->S_gr[media]);
}

__device__ double device_P_k_gn_S(double S)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	return gpu_def->P_d_gn[media] * pow(pow(1. - S, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S, A / (1. - A) - 1.) / (1. - A)
		/ (1. - gpu_def->S_wr[media] - gpu_def->S_nr[media] - gpu_def->S_gr[media]);
}

// ������� ���������� ����������� �������� � ����������� �� ���� ���������
// �� ����� ��������� [0, 1] ������� ����������� �������� ������ �������� ���������, ����������� �������� ��������������.
// �������� ����� ���������� � ����� mathcad.
__device__ double device_assign_P_k_nw(double S_w_e)
{
	double Pk_nw = 0;
/*	double S_w_range[2];
	device_assign_S_w_range(S_w_range, S_w_range + 1);

	if (S_w_e <= S_w_range[0])
	{
		Pk_nw = device_P_k_nw_S(S_w_range[0]) * (S_w_e - S_w_range[0]) + device_P_k_nw(S_w_range[0]);
	}
	else if (S_w_e >= S_w_range[1])
	{
		Pk_nw = device_P_k_nw_S(S_w_range[1]) * (S_w_e - S_w_range[1]) + device_P_k_nw(S_w_range[1]);;
	}
	else
	{
		Pk_nw = device_P_k_nw(S_w_e);
	}
*/
	return Pk_nw;
}

__device__ double device_assign_P_k_gn(double S_g_e)
{
	double Pk_gn = 0;
/*	double S_g_range[2];
	device_assign_S_g_range(S_g_range, S_g_range + 1);

	if (S_g_e <= S_g_range[0])
	{
		Pk_gn = device_P_k_gn_S(S_g_range[0]) * (S_g_e - S_g_range[0]) + device_P_k_gn(S_g_range[0]);
	}
	else if (S_g_e >= S_g_range[1])
	{
		Pk_gn = device_P_k_gn_S(S_g_range[1]) * (S_g_e - S_g_range[1]) + device_P_k_gn(S_g_range[1]);
	}
	else
	{
		Pk_gn = device_P_k_gn(S_g_e);
	}
*/
	return Pk_gn;
}

// ������� ���������� ����������� ����������� �������� �� �������������
__device__ double device_assign_P_k_nw_S(double S_w_e)
{
	double PkSw = 0;
/*	double S_w_range[2];
	device_assign_S_w_range(S_w_range, S_w_range + 1);

	if (S_w_e <= S_w_range[0])
	{
		PkSw = device_P_k_nw_S(S_w_range[0]);
	}
	else if (S_w_e >= S_w_range[1])
	{
		PkSw = device_P_k_nw_S(S_w_range[1]);
	}
	else
	{
		PkSw = device_P_k_nw_S(S_w_e);
	}
*/
	return PkSw;
}

__device__ double device_assign_P_k_gn_S(double S_g_e)
{
	double PkSn = 0;
/*	double S_g_range[2];
	device_assign_S_g_range(S_g_range, S_g_range + 1);

	if (S_g_e <= S_g_range[0])
	{
		PkSn = (-1) * device_P_k_gn_S(S_g_range[0]);
	}
	else if (S_g_e >= S_g_range[1])
	{
		PkSn = (-1) * device_P_k_gn_S(S_g_range[1]);
	}
	else
	{
		PkSn = device_P_k_gn_S(S_g_e);
	}
*/
	return PkSn;
}

// ������� ���������� ������������� ��������������
__device__ double device_assign_k_w(double S_w_e)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	double k_w = 0;

	if (S_w_e >= 1e-3)
	{
		k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)), (A - 1.) / A), 2.);
	}

	return k_w;
}

__device__ double device_assign_k_g(double S_g_e)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	double k_g = 0;

	if (S_g_e >= 1e-3)
	{
		k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);
	}

	return k_g;
}

__device__ double device_assign_k_n(double S_w_e, double S_n_e)
{
	int media = 0;
	double A = gpu_def->lambda[media];
	double k_n = 0;
	double S_g_e = 1. - S_w_e - S_n_e;

	if (S_n_e >= 1e-3)
	{
		double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);
		double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A / (A - 1.)), (A - 1.) / A), 2.);
		k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e);
	}

	return k_n;
}

//������� ���������� �������� ��������, ���������� � ������������� � ������ ����� � ����� (i,j,k) ����� media,
//������ �� ��������� �������� �������� ���������� (Pw,Sw,Sn)
//1. ����������, � ����� ������ �� ���� ��������
//2. ���������� �������� ������������ ���� n �� ������� ��������� ������������� � ����� �������
//3. ���������� ����������� ������������� �� �������� ������ ���������� ����������
//4. ���������� ������������� ������� �������������� � ������������ � ������������ ������ � ����������� ����� � �������
//5. ���������� ����������� �������� � ������������ � ������������ ������� �������
//6. ���������� ������� �������� c ������� �����������
//7. ���������� ������������� ������ �����

__global__ void assign_P_Xi_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int media = 0;
		double k_w, k_g, k_n, Pk_nw, Pk_gn;
		double S_w_e = device_assign_S_w_e(DevArraysPtr, i, j, k);
		double S_n_e = device_assign_S_n_e(DevArraysPtr, i, j, k);
		double S_g_e = 1. - S_w_e - S_n_e;

		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		k_w = device_assign_k_w(S_w_e);
		k_g = device_assign_k_g(S_g_e);
		k_n = device_assign_k_n(S_w_e, S_n_e);

		DevArraysPtr.Xi_w[local] = (-1.) * (gpu_def->K[media]) * k_w / gpu_def->mu_w;
		DevArraysPtr.Xi_n[local] = (-1.) * (gpu_def->K[media]) * k_n / gpu_def->mu_n;
		DevArraysPtr.Xi_g[local] = (-1.) * (gpu_def->K[media]) * k_g / gpu_def->mu_g;

		if ((i != 0) && (i != (gpu_def->locNx) - 1) && (j != 0) && (j != (gpu_def->locNy) - 1) && (((k != 0) && (k != (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2)))
		{
			Pk_nw = device_assign_P_k_nw(S_w_e);
			Pk_gn = device_assign_P_k_gn(S_g_e);

			DevArraysPtr.P_n[local] = DevArraysPtr.P_w[local] + Pk_nw;
			DevArraysPtr.P_g[local] = DevArraysPtr.P_n[local] + Pk_gn;
		}

		device_test_positive(DevArraysPtr.P_n[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.P_g[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_w[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_n[local], __FILE__, __LINE__);
		device_test_nan(DevArraysPtr.Xi_g[local], __FILE__, __LINE__);
	}
}

// ��������������� ������� ��� ������ �������:
// ���������� �������� ������� 3*3;
__device__ void device_reverse_matrix(double* a)
{
	int n = 3;
	double b[9], det = 0;

	// ���������� �������������� ������� �������
	for(int j = 0; j < n; j++)
		for(int i = 0; i < n; i++)
		{
			b[i + n * j] = a[(i + 1) % n + n * ((j + 1) % n)] * a[(i + 2) % n + n * ((j + 2) % n)]
			- a[(i + 2) % n + n * ((j + 1) % n)] * a[(i + 1) % n + n * ((j + 2) % n)];
		}

		// ���������� ������������ ������� 3*3;
		for(int i = 0; i < n; i++)
		{
			det += a[i] * b[i];
		}
		device_test_nan(det, __FILE__, __LINE__);

		// ���������������� � ������� �� �����������
		for(int j = 0; j < n; j++)
			for(int i = 0; i < n; i++)
			{
				a[i + n * j] = b[j + n * i] / det;
				device_test_nan(a[i + n * j], __FILE__, __LINE__);
			}
}

//������� ������� ������� 3*3 �� �������� ��������� (Pn,Sw,Sg) ������� ������� � ����� (i,j,k) ����� media
//1. ���������� ����������� �������������
//2. ��������������� ���������� ������������
//3. ���������� ����������� ��������
//4. ���������� ������������ ���� n
//5. ���������� �������� ���� ������� �������
//6. ���������� ������� ����������� ����������� ����������� �������� �� �������������
//7. ���������� ������� ������� �����������
//8. ���������� ������������ ������� ������� �����������
//9. ��������� ������� ������� ������� ������� � ����� ����
__global__ void Newton_method_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i != 0) && (i < (gpu_def->locNx) - 1) && (j != 0) && (j < (gpu_def->locNy) - 1) && (((k != 0) && (k < (gpu_def->locNz) - 1)) || ((gpu_def->locNz) < 2))
		&& (device_is_active_point(i, j, k) == 1))
	{
		double S_w_e, S_g_e, S_n_e, Pk_nw, Pk_gn, PkSw, PkSn, Sg, F1, F2, F3;
		double dF[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);

		for (int w = 1; w <= gpu_def->newton_iterations; w++)
		{
			S_w_e = device_assign_S_w_e(DevArraysPtr, i, j, k);
			S_n_e = device_assign_S_n_e(DevArraysPtr, i, j, k);
			S_g_e = 1. - S_w_e - S_n_e;

			Pk_nw = device_assign_P_k_nw(S_w_e);
			Pk_gn = device_assign_P_k_gn(S_g_e);
			PkSw = device_assign_P_k_nw_S(S_w_e);
			PkSn = device_assign_P_k_gn_S(S_g_e);

			Sg = 1. - DevArraysPtr.S_w[local] - DevArraysPtr.S_n[local];

			F1 = gpu_def->ro0_w * (1. + (gpu_def->beta_w) * (DevArraysPtr.P_w[local] - gpu_def->P_atm))
			     * DevArraysPtr.S_w[local] - DevArraysPtr.roS_w[local];
			F2 = gpu_def->ro0_n * (1. + (gpu_def->beta_n) * (DevArraysPtr.P_w[local] + Pk_nw - gpu_def->P_atm))
			     * DevArraysPtr.S_n[local] - DevArraysPtr.roS_n[local];
			F3 = gpu_def->ro0_g * (DevArraysPtr.P_w[local] + Pk_nw + Pk_gn) / gpu_def->P_atm
			     * Sg - DevArraysPtr.roS_g[local];

			// ������� ������� �����������. ������� �� 0 �� 8 ������������� F1P, F1Sw, F1Sn, F2P, F2Sw, F2Sn, F3P, F3Sw, F3Sn
			dF[0] = gpu_def->ro0_w * gpu_def->beta_w * DevArraysPtr.S_w[local];
			dF[3] = gpu_def->ro0_n * gpu_def->beta_n * DevArraysPtr.S_n[local];
			dF[6] = gpu_def->ro0_g * Sg / gpu_def->P_atm;
			dF[1] = gpu_def->ro0_w * (1 + gpu_def->beta_w * (DevArraysPtr.P_w[local] - gpu_def->P_atm));
			dF[4] = gpu_def->ro0_n * (1. + (gpu_def->beta_n) * PkSw) * DevArraysPtr.S_n[local];
			dF[7] = (-1) * gpu_def->ro0_g * (DevArraysPtr.P_w[local] + Pk_nw + Pk_gn - Sg * (PkSn + PkSw)) / gpu_def->P_atm;
			dF[2] = 0;
			dF[5] = gpu_def->ro0_n * (1. + gpu_def->beta_n * (DevArraysPtr.P_w[local] + Pk_nw - gpu_def->P_atm));
			dF[8] = (-1) * gpu_def->ro0_g * (DevArraysPtr.P_w[local] + Pk_nw + Pk_gn - Sg * PkSn) / gpu_def->P_atm;

			device_reverse_matrix(dF);

			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[local]
			        - (dF[0] * F1 + dF[1] * F2 + dF[2] * F3);
			DevArraysPtr.S_w[local] = DevArraysPtr.S_w[local]
			        - (dF[3] * F1 + dF[4] * F2 + dF[5] * F3);
			DevArraysPtr.S_n[local] = DevArraysPtr.S_n[local]
			        - (dF[6] * F1 + dF[7] * F2 + dF[8] * F3);
		}

		device_test_S(DevArraysPtr.S_w[local], __FILE__, __LINE__);
		device_test_S(DevArraysPtr.S_n[local], __FILE__, __LINE__);
		device_test_positive(DevArraysPtr.P_w[local], __FILE__, __LINE__);
	}
}

//������� ��������� ������� �������� ��� (Sw,Sg),Pn

// ������� ��������� ������� � ������� ������ ��������, �� � ��������� �������������� ����������
__global__ void Border_S_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int i1 = i, j1 = j, k1 = k;

		device_set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1);

		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		int local1 = i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy);

		if ((j != 0) || ((gpu_def->source) <= 0))
		{
			DevArraysPtr.S_w[local] = DevArraysPtr.S_w[local1];
			DevArraysPtr.S_n[local] = DevArraysPtr.S_n[local1];
		}

		if ((j == 0) && ((gpu_def->source) > 0))
		{
			DevArraysPtr.S_w[local] = gpu_def->S_w_gr;
			DevArraysPtr.S_n[local] = gpu_def->S_n_gr;
		}
	}
}

__global__ void Border_P_kernel(ptr_Arrays DevArraysPtr)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int k = threadIdx.z + blockIdx.z * blockDim.z;

	if ((i < (gpu_def->locNx)) && (j < (gpu_def->locNy)) && (k < (gpu_def->locNz)) && (device_is_active_point(i, j, k) == 1))
	{
		int i1 = i, j1 = j, k1 = k;

		device_set_boundary_basic_coordinate(i, j, k, &i1, &j1, &k1);

		int local = i + j * (gpu_def->locNx) + k * (gpu_def->locNx) * (gpu_def->locNy);
		int local1 = i1 + j1 * (gpu_def->locNx) + k1 * (gpu_def->locNx) * (gpu_def->locNy);

		double S_w_e = device_assign_S_w_e(DevArraysPtr, i1, j1, k1);
		double S_n_e = device_assign_S_n_e(DevArraysPtr, i1, j1, k1);
		double S_g_e = 1. - S_w_e - S_n_e;

		double P_k_nw = device_assign_P_k_nw(S_w_e);
		double P_k_gn = device_assign_P_k_gn(S_g_e);

		// ���� �������� ������ �������� �� �������� ����� �������� (������� ������������)
		if ((j != 0) && (j != (gpu_def->locNy) - 1))
		{
			DevArraysPtr.P_w[local] = DevArraysPtr.P_w[local1];
			DevArraysPtr.P_n[local] = DevArraysPtr.P_w[local1] + P_k_nw;
			DevArraysPtr.P_g[local] = DevArraysPtr.P_w[local1] + P_k_nw + P_k_gn;
		
		}
		else if (j == 0)
		{
			DevArraysPtr.P_w[local] = (DevArraysPtr.P_w[local1]
			- (gpu_def->ro0_w) * (gpu_def->g_const) * (gpu_def->hy) * (1. - (gpu_def->beta_w) * (gpu_def->P_atm))) 
				/ (1. + (gpu_def->beta_w) * (gpu_def->ro0_w) * (gpu_def->g_const) * (gpu_def->hy));
			DevArraysPtr.P_n[local] = (DevArraysPtr.P_w[local1]
			+ P_k_nw - (gpu_def->ro0_n) * (gpu_def->g_const) * (gpu_def->hy) * (1. - (gpu_def->beta_n) * (gpu_def->P_atm))) 
				/ (1. + (gpu_def->beta_n) * (gpu_def->ro0_n) * (gpu_def->g_const) * (gpu_def->hy));
			DevArraysPtr.P_g[local] = (DevArraysPtr.P_w[local1]
			+ P_k_nw + P_k_gn) / (1. + (gpu_def->ro0_g) * (gpu_def->g_const) * (gpu_def->hy) / (gpu_def->P_atm));
		}
		else
		{
			DevArraysPtr.P_w[local] = (DevArraysPtr.P_w[local1]
			+ (gpu_def->ro0_w) * (gpu_def->g_const) * (gpu_def->hy) * (1. - (gpu_def->beta_w) * (gpu_def->P_atm))) 
				/ (1. - (gpu_def->beta_w) * (gpu_def->ro0_w) * (gpu_def->g_const) * (gpu_def->hy));
			DevArraysPtr.P_n[local] = (DevArraysPtr.P_w[local1]
			+ P_k_nw + (gpu_def->ro0_n) * (gpu_def->g_const) * (gpu_def->hy) * (1. - (gpu_def->beta_n) * (gpu_def->P_atm))) 
				/ (1. - (gpu_def->beta_n) * (gpu_def->ro0_n) * (gpu_def->g_const) * (gpu_def->hy));
			DevArraysPtr.P_g[local] = (DevArraysPtr.P_w[local1]
			+ P_k_nw + P_k_gn) / (1. - (gpu_def->ro0_g) * (gpu_def->g_const) * (gpu_def->hy) / (gpu_def->P_atm));
		}
	}
}

void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					// �������������� ��������� ��������� ���������� � ����������
					int I = local_to_global(i, 'x', def);
					int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

					HostArraysPtr.m[local]=def.porosity[0];
					// �������� ��������� ������������� � ��������� �������������
	/*				int j1 = def.locNy / 2;

					if (j < j1)
					{
						HostArraysPtr.S_w[local] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
						HostArraysPtr.S_n[local] = def.S_n_gr + (def.S_n_init - def.S_n_gr) * j / j1;
					}
					else
					*/
					if ((j == 0) && ((def.source) > 0))
					{
						HostArraysPtr.S_w[local] = def.S_w_gr;
						HostArraysPtr.S_n[local] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_w[local] = def.S_w_init;
						HostArraysPtr.S_n[local] = def.S_n_init;
					}

					double ro_g_dy = (def.ro0_g * (1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local])
					+ def.ro0_w * HostArraysPtr.S_w[local]
					+ def.ro0_n * HostArraysPtr.S_n[local]) * (HostArraysPtr.m[local]) * (def.g_const) * (def.hy);

					// ���� �������� ������ �������� �� �������� ����� ��������
					if (j == 0)
					{
						HostArraysPtr.P_w[local] = def.P_atm;
						HostArraysPtr.P_n[local] = def.P_atm;
						HostArraysPtr.P_g[local] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[local] = HostArraysPtr.P_w[local - (def.locNx)] + ro_g_dy;
						HostArraysPtr.P_n[local] = HostArraysPtr.P_n[local - (def.locNx)] + ro_g_dy;
						HostArraysPtr.P_g[local] = HostArraysPtr.P_g[local - (def.locNx)] + ro_g_dy;
					}

					HostArraysPtr.ro_w[local] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[local] - def.P_atm));
					HostArraysPtr.ro_n[local] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[local] - def.P_atm));
					HostArraysPtr.ro_g[local] = def.ro0_g * HostArraysPtr.P_g[local] / def.P_atm;

					test_S(HostArraysPtr.S_n[local], __FILE__, __LINE__);
					test_S(HostArraysPtr.S_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_n[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.P_g[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_w[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_n[local], __FILE__, __LINE__);
					test_positive(HostArraysPtr.ro_g[local], __FILE__, __LINE__);
				}
}

