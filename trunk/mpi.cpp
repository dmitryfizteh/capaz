#include <mpi.h>
#include "defines.h"

// �������� � ����� ������ ������ �������
void right_send_recv(double* HostBuffer, int destination_rank, int send_recv_id, int localNx, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer+(def.Ny)*(def.Nz),(def.Ny)*(def.Nz),MPI_DOUBLE,destination_rank,send_recv_id,destination_rank,send_recv_id+1,MPI_COMM_WORLD,&status)==MPI_SUCCESS)
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\n");
}

// ��������� � �������� ������ �� ����� �������
void left_recv_send(double* HostBuffer, int destination_rank, int send_recv_id, int localNx, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer,(def.Ny)*(def.Nz),MPI_DOUBLE,destination_rank,send_recv_id+1,destination_rank,send_recv_id,MPI_COMM_WORLD,&status)==MPI_SUCCESS)
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\n");
}

// ����� ������� �� �������� ����� ����� ������������
// 0. ��������� ������ � ��������� � ������ �����
// 1.  ��� ���� ������ ����������� 
// 1.1 ��������/�������� ������ �������, 
// 1.2 ��������/�������� ����� �������.
// 2.2 ��� �������� - ��������/�������� ����� �������,
// 2.2 ��������/�������� ������.
// ��� ������� ����������� ��������������� ������ �� ���������
// 3. ��������� ���������� ������ � ������ ����������
void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
	load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, localNx, blocksY, blocksZ, rank, size, def); // (0)

	if(rank%2 == 0) // (1)
	{
		if (rank!=size-1)
			right_send_recv(HostBuffer, rank+1, 500, localNx, def); // (1.1)

		if (rank!=0)
			left_recv_send(HostBuffer, rank-1, 502, localNx, def); // (1.2)
	}
	else
	{
		if (rank!=0) // � ��������, ������ ��������
			left_recv_send(HostBuffer, rank-1, 500, localNx, def); // (2.1)

		if (rank!=size-1)
			right_send_recv(HostBuffer, rank+1, 502, localNx, def); // (2.2)
	}

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
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,size); // The amount of processors
    MPI_Comm_rank(MPI_COMM_WORLD,rank); // The number of processor
	std::cout << "size=" <<*size<<"  "<<"rank= "<<*rank<<"\n";
}

void communication_finalization(void)
{
	MPI_Finalize();
}

// ���������� ������ Barrier ��� ��������� ������������
void barrier(void)
{
	MPI_Barrier(MPI_COMM_WORLD);
}