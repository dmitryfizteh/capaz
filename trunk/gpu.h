#include "defines.h"
#include <cuda.h>

void gpu_test_nan(double* HostArrayPtr, double* DevArrayPtr, int rank, int size, int localNx, consts def, char* file, int line);
void checkErrors(char *label);


