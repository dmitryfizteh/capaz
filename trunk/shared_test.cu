#include "gpu.h"
#include "cuPrintf.cu"

// �������� ������ GPU
void checkErrors(char *label, char *file, int line)
{
#ifdef MY_TEST
	cudaError_t err;

	err = cudaThreadSynchronize();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\nFile:\"%s\"\nLine:\"%d\"\n\n", e, label, file, line);
	}

	err = cudaGetLastError();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\nFile:\"%s\"\nLine:\"%d\"\n\n", e, label, file, line);
		fflush(stdout);
	}
#endif
}

// �������, ���������� ��� ������
__device__ void device_print_error(char *error, char *file, int line)
{
	CUPRINTF("Error: %s\nFile: \"%s\"\nLine: %d\n\n", error, file, line);
}

// ���� �� NaN
// ��������� ������ device_test_nan(x, __FILE__, __LINE__);
__device__ void device_test_nan(double x, char *file, int line)
{
#ifdef MY_TEST
	if ((x > 1e+30) || (x < -1 * 1e+40))
	{
		CUPRINTF("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ������������� � �� NaN
// ��������� ������ device_test_positive(x, __FILE__, __LINE__);
__device__ void device_test_positive(double x, char *file, int line)
{
#ifdef MY_TEST
	if ((x > 1e+30) || (x < 0))
	{
		CUPRINTF("Error: NaN or X<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ��������� ������������� � [0;1]
// ��������� ������ device_test_S(x, __FILE__, __LINE__);
__device__ void device_test_S(double S, char *file, int line)
{
#ifdef MY_TEST
	if (S < 0)
	{
		CUPRINTF("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (S > 1)
	{
		CUPRINTF("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ��������� ��������� � [-100;100]
// ��������� ������ test_u(x, __FILE__, __LINE__);
__device__ void device_test_u(double u, char *file, int line)
{
#ifdef MY_TEST
	if (u < -1e8)
	{
		CUPRINTF("Error: u<-100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (u > 1e8)
	{
		CUPRINTF("Error: u>100\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}

// ���� �� ��������� ���������� � [0;3000]
// ��������� ������ test_ro(x, __FILE__, __LINE__);
__device__ void device_test_ro(double ro, char *file, int line)
{
#ifdef MY_TEST
	if (ro < 0)
	{
		CUPRINTF("Error: ro < 0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
	if (ro > 3000)
	{
		CUPRINTF("Error: ro > 5000\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	}
#endif
}