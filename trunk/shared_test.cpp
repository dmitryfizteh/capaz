#include "defines.h"

// ������������

// �������� �� ��������������� ������� ������
#ifndef THREE_PHASE
#ifndef TWO_PHASE
#ifndef B_L
#error "Task DEFINE is not found."
#endif
#endif
#endif


// �������, ���������� ��� ������
void print_error(char *error, char *file, int line)
{
	printf("Error: %s\nFile: \"%s\"\nLine: %d\n\n", error, file, line);
	fflush(stdout);
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif
	exit(1);
}

// ������� �������� �� ����� �� ����������� ��������� ��������
// �� ���� ������ ��������� ������� ����������
#ifndef THREE_PHASE
void test_correct_P_S(ptr_Arrays HostArraysPtr, consts def)
{
	for (int i = 0; i < (def.locNx); i++)
		for (int j = 0; j < (def.locNy); j++)
			for (int k = 0; k < (def.locNz); k++)
			{
				test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
				test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ux_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uy_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uz_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
			}
}
#endif

// ���� �� NaN
// ��������� ������ test_nan(x, __FILE__, __LINE__);
void test_nan(double x, char *file, int line)
{
#ifdef MY_TEST
	if (isnan(x))
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	// �������� ����� ������� ����������� ������ ��� ���� ������
	if (x > 1e12)
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ������������� � �� NaN
// ��������� ������ test_positive(x, __FILE__, __LINE__);
void test_positive(double x, char *file, int line)
{
#ifdef MY_TEST
	if (isnan(x))
	{
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (x < 0)
	{
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ��������� ������������� � [0;1]
// ��������� ������ test_S(x, __FILE__, __LINE__);
void test_S(double S, char *file, int line)
{
#ifdef MY_TEST
	if (isnan(S))
	{
		printf("Error: S=NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S < 0)
	{
		printf("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S > 1)
	{
		printf("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ��������� ��������� � [-100;100]
// ��������� ������ test_u(x, __FILE__, __LINE__);
void test_u(double u, char *file, int line)
{
#ifdef MY_TEST
	if (isnan(u))
	{
		printf("Error: u=NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u < -10)
	{
		printf("Error: u<-100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u > 10)
	{
		printf("Error: u>100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ������� ���������, ��� ������ �������� ����� ������ (�� ������) �������
// ���� ��� �� ���, ���������� ��������������
void test_arrowhead(double big, double small, char *file, int line)
{
#ifdef MY_TEST_1
	if (fabs(big / 30) < fabs(small))
	{
		printf("Warning: See task parameters.\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ������� ���������, ��� ������ ��������� ��������� ������������� ����� ������ (�� ������) �������
// ���� ��� �� ���, ���������� ��������������
void test_tau(double S_old, double S_now, double S_new, int media, consts def, char *file, int line)
{
#ifdef MY_TEST_1
	double L = def.m[media] * (S_new - S_old) / (2 * (def.dt));
	double R = def.tau * (S_new - 2 * S_now + S_old) / ((def.dt) * (def.dt));

	if (fabs(L / 30) < fabs(R))
	{
		printf("Warning: parameter tau is very much.\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ������������ ���������� ������
void read_defines_test(consts def)
{
#ifdef MY_TEST
	test_positive(def.hx, __FILE__, __LINE__);
	test_positive(def.hy, __FILE__, __LINE__);
	test_positive(def.hz, __FILE__, __LINE__);
	test_positive(def.tau, __FILE__, __LINE__);
	test_positive(def.dt, __FILE__, __LINE__);
	test_positive(def.c_w, __FILE__, __LINE__);
	test_positive(def.c_n, __FILE__, __LINE__);
	test_positive(def.l, __FILE__, __LINE__);
	test_positive(def.beta_w, __FILE__, __LINE__);
	test_positive(def.beta_n, __FILE__, __LINE__);
	test_positive(def.ro0_w, __FILE__, __LINE__);
	test_positive(def.ro0_n, __FILE__, __LINE__);
	test_positive(def.mu_w, __FILE__, __LINE__);
	test_positive(def.mu_n, __FILE__, __LINE__);
	test_positive(def.g_const, __FILE__, __LINE__);
	test_positive(def.P_atm, __FILE__, __LINE__);

	test_positive(def.source, __FILE__, __LINE__);
	test_positive(def.newton_iterations, __FILE__, __LINE__);
	test_positive(def.timeX, __FILE__, __LINE__);
	test_positive(def.save_plots, __FILE__, __LINE__);
	test_positive(def.print_screen, __FILE__, __LINE__);
	test_positive(def.Nx, __FILE__, __LINE__);
	test_positive(def.Ny, __FILE__, __LINE__);
	test_positive(def.Nz, __FILE__, __LINE__);

#ifdef B_L
	test_positive(def.Q, __FILE__, __LINE__);
#endif
#ifdef THREE_PHASE
	test_positive(def.c_g, __FILE__, __LINE__);
	test_positive(def.beta_g, __FILE__, __LINE__);
	test_positive(def.ro0_g, __FILE__, __LINE__);
	test_positive(def.mu_g, __FILE__, __LINE__);
	test_positive(def.S_w_gr, __FILE__, __LINE__);
#endif
	test_positive(def.S_n_gr, __FILE__, __LINE__);
#endif

}

