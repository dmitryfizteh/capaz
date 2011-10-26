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
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
}

// Тест на положительное целое и не NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_pint (int x, char *file, int line)
{
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( x < 0 )
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
}
