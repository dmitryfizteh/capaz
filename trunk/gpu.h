#include "defines.h"
#include <cuda.h>

__constant__ consts gpu_def [1];

void checkErrors(char *label);
