#include <shmem++.h>
#include "defines.h"
//#include <shmem.h>

#define MX 20000
static double *pf, *cf[MX];

// ����� ������� �� �������� ����� ����� ������������
// 0. ��������� ������ � ��������� � ������ �����
// 1.1 �������� ������ �������, 
// 1.2 �������� ����� �������.
// ��� ������� ����������� ��������������� ������ �� ���������
// 3. ��������� ���������� ������ � ������ ����������
void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
	load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def); // (0)

	if (rank!=size-1)
		shmem_double_put (cf[rank+1],HostBuffer+(def.Ny)*(def.Nz),sizeof(double)*(def.Ny)*(def.Nz),rank+1); // (1.1)
	if (rank!=0)
		shmem_double_put (cf[rank-1]+(def.Ny)*(def.Nz),HostBuffer,sizeof(double)*(def.Ny)*(def.Nz),rank-1); // (1.2)
	shmem_barrier_all();
	HostBuffer=cf[rank];

	save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def); // (3)
}

// ����� ���������� ���������� �������� P2, ���������� ro1 � ro2, Xi ����� ������������
void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ro_g, DevArraysPtr.ro_g, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.Xi_g, DevArraysPtr.Xi_g, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.P_g, DevArraysPtr.P_g, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#else
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#endif
	exchange(HostArraysPtr.ro_w, DevArraysPtr.ro_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.ro_n, DevArraysPtr.ro_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.Xi_w, DevArraysPtr.Xi_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.Xi_n, DevArraysPtr.Xi_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
}


// ����� ���������� ���������� ��������� ����� ������������
// � ������ ������������� ��������� ������� ����� ������������ �� ��� X
// �������� u1y � u2y �� ���������
void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
	exchange(HostArraysPtr.ux_w, DevArraysPtr.ux_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.ux_n, DevArraysPtr.ux_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ux_g, DevArraysPtr.ux_g, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#endif
}

// ����� ���������� ���������� �������� ���� P1 � ������������ NAPL S2 ����� ������������
void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.S_g, DevArraysPtr.S_g, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#else
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def);
#endif
}

void communication_initialization(int argc, char* argv[], int* size, int* rank, consts def)
{
	shmem_init(&argc, &argv);
	*size=shmem_n_pes(); // The amount of processors
	*rank=shmem_my_pe(); // The number of processor
	std::cout << "size=" <<*size<<"  "<<"rank= "<<*rank<<"\n";

	pf = (double*)emalloc( sizeof(double)*(def.Ny)*(def.Nz));
	shmem_coarray_all( (void*)pf, (long)(def.Ny*(def.Nz)*sizeof(*pf)), (void**)cf );
}

void communication_finalization(void)
{
	shmem_finalize();
}

// ���������� ������ Barrier ��� ��������� ������������
void barrier(void)
{
	shmem_barrier_all();
}