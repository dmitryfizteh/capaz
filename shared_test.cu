#include "gpu.h"
#include "cuPrintf.cu"

// Проверка ошибок GPU
void checkErrors(char *label) 
{
#ifdef MY_TEST
	cudaError_t err;

	err = cudaThreadSynchronize();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\n", e, label);
	}

	err=cudaGetLastError();
	if (err != cudaSuccess)
	{
		char *e = (char*) cudaGetErrorString(err);
		printf("CUDA Error: %s (at %s)\n", e, label);
		fflush(stdout);
	}
#endif
}

// Тест на NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
__device__ void device_test_nan (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( (x>1e+30) || (x<-1*1e+40))
		CUPRINTF("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}

// Тест на положительное и не NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
__device__ void device_test_positive (double x, char *file, int line)
{
#ifdef MY_TEST
	if ( (x>1e+30) || (x<0))
		CUPRINTF("Error: NaN or X<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
#endif
}