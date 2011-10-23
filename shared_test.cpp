#include "defines.h"

// Тестовая функция
double Factorial (double x)
{
	return x;
}
//
#ifdef MY_TEST
namespace {

	// Проверить факториал от 0.
	TEST(FactorialTest, HandlesZeroInput) {
		EXPECT_EQ(0, Factorial(0));
	}

	// Проверить факториал некоторых положительных значений.
	TEST(FactorialTest, HandlesPositiveInput) {
		EXPECT_EQ(1, Factorial(1));
		EXPECT_EQ(2, Factorial(2));
		EXPECT_EQ(3, Factorial(3));
		EXPECT_EQ(8, Factorial(8));
	}

}  // namespace
#endif

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
				if (isnan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! S_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! P_w is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! ux_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! uy_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! uz_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.P_n [i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! P_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.ro_n [i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! ro_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (isnan(HostArraysPtr.Xi_n  [i+j*localNx+k*localNx*(def.Ny)]))
					printf ("\nWarning! Xi_n is Nan in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
			}
}
#endif

#ifndef THREE_PHASE
void NaN_test (ptr_Arrays DevArraysPtr, ptr_Arrays HostArraysPtr, int localNx, int rank, consts def)
{

}
#endif