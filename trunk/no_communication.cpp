#include "defines.h"

void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
}

void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
}

void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
}

void communication_initialization(int argc, char* argv[], consts *def)
{
	(*def).size = 1;
	(*def).rank = 0;
}

void communication_finalization(void)
{
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
}

