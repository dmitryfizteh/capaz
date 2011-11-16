#include "defines.h"

void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
}

void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
}

void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
}

void communication_initialization(int argc, char* argv[], int* size, int* rank, consts def)
{
	*size=1;
	*rank=0;
}

void communication_finalization(void)
{
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
}