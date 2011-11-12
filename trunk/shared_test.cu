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

	err=cudaGetLastError();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\nFile:\"%s\"\nLine:\"%d\"\n\n", e, label, file, line);
		fflush(stdout);
	}
#endif
}

// ���� �� NaN
// ��������� ������ device_test_nan(x, __FILE__, __LINE__);
__device__ void device_test_nan (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( (x>1e+30) || (x<-1*1e+40))
		CUPRINTF("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// ���� �� ������������� � �� NaN
// ��������� ������ device_test_positive(x, __FILE__, __LINE__);
__device__ void device_test_positive (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( (x>1e+30) || (x<0))
		CUPRINTF("Error: NaN or X<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// ���� �� ��������� ������������� � [0;1]
// ��������� ������ device_test_S(x, __FILE__, __LINE__);
__device__ void device_test_S (double S, char *file, int line)
{
#ifdef MY_TEST
	if ( S < 0 )
		CUPRINTF("Error: S<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( S > 1 )
		CUPRINTF("Error: S>1\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}