#include "defines.h"

void Pn_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
}

void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
}

void Pw_Sn_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
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