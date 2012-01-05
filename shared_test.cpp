#include "defines.h"

// Тестирование

// Проверка на предопределение макроса задачи
#ifndef THREE_PHASE
#ifndef TWO_PHASE
#ifndef B_L
#error "Task DEFINE is not found."
#endif
#endif
#endif

// Функция проверки на выход из допустимого диапазона значений
// во всех точках расчетной области процессора
#ifndef THREE_PHASE
void test_correct_P_S(ptr_Arrays HostArraysPtr, localN locN, int rank, consts def)
{
	for(int i=0;i<(locN.x);i++)
		for(int j=0;j<(locN.y);j++)
			for(int k=0;k<(locN.z);k++)
			{
				test_S(HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
				test_positive(HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
			}
}
#endif

// Тест на NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_nan (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Тест на положительное и не NaN
// Синтаксис вызова test_positive(x, __FILE__, __LINE__);
void test_positive (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( x < 0 )
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Тест на вхождение насыщенностей в [0;1]
// Синтаксис вызова test_S(x, __FILE__, __LINE__);
void test_S (double S, char *file, int line)
{
#ifdef MY_TEST
	if ( isnan(S) )
		printf("Error: S=NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( S < 0 )
		printf("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( S > 1 )
		printf("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Функция проверяет, что первый аргумент много больше (по модулю) второго
// Если это не так, печатается предупреждение
void test_arrowhead(double big, double small, char *file, int line)
{
#ifdef MY_TEST_1
	if (abs(big/30) < abs(small))
		printf("Warning: See task parameters.\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Функция проверяет, что первое слагаемое уравнения неразрывности много больше (по модулю) второго
// Если это не так, печатается предупреждение
void test_tau(double S_old, double S_now, double S_new, int media, consts def, char *file, int line)
{
#ifdef MY_TEST_1
	double L=def.m[media] * (S_new - S_old) / (2 * (def.dt));
	double R=def.tau * (S_new - 2*S_now + S_old) / ((def.dt)*(def.dt));

	if (abs(L/30) < abs(R))
		printf("Warning: parameter tau is very much.\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Тест на корректность параметров задачи
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