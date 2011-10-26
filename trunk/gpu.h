#include "defines.h"
#include <cuda.h>

#ifdef MY_TEST
__device__ void test_nan (double x, char *file, int line);
__device__ void test_pint (int x, char *file, int line);
void checkErrors(char *label);
#endif


