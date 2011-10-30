#include "defines.h"

// Тестирование

// Функция проверки на выход из допустимого диапазона значений
// во всех точках расчетной области процессора
#ifndef THREE_PHASE
void test_correct_P_S(ptr_Arrays HostArraysPtr, int localNx, int rank, consts def)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
			{
				if (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]<0)
					printf ("\nWarning! S2<0 in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]<=0)
					printf ("\nWarning! P<=0 in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);

				test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
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

// Тест на положительное целое и не NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_positive (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( x < 0 )
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
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
	
#ifdef THREE_PHASE
	test_positive(def.c_g, __FILE__, __LINE__);
	test_positive(def.beta_g, __FILE__, __LINE__);
	test_positive(def.ro0_g, __FILE__, __LINE__);
	test_positive(def.mu_g, __FILE__, __LINE__);
	test_positive(def.S_w_gr, __FILE__, __LINE__);
	test_positive(def.S_g_gr, __FILE__, __LINE__);
#else
	test_positive(def.S_n_gr, __FILE__, __LINE__);
#endif

#endif
	
}