#include <mpi.h>
#include "defines.h"

// Передача и прием данных правой границе
void right_send_recv(double* HostBuffer, int destination_rank, int send_recv_id, localN locN, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer+(locN.y)*(locN.z),(locN.y)*(locN.z),MPI_DOUBLE,destination_rank,send_recv_id,destination_rank,send_recv_id+1,MPI_COMM_WORLD,&status)==MPI_SUCCESS)
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
}

// Получение и передача данных на левой границе
void left_recv_send(double* HostBuffer, int destination_rank, int send_recv_id, localN locN, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer,(locN.y)*(locN.z),MPI_DOUBLE,destination_rank,send_recv_id+1,destination_rank,send_recv_id,MPI_COMM_WORLD,&status)==MPI_SUCCESS)
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
}

// Обмен данными на границах между всеми процессорами
// 0. Загружаем данные с усорителя в память хоста
// 1.  Для всех четных процессоров 
// 1.1 передаем/получаем правую границу, 
// 1.2 получаем/передаем левую границу.
// 2.2 Для нечетных - получаем/передаем левую границу,
// 2.2 передаем/получаем правую.
// Для крайних процессоров соответствующие обмены не требуются
// 3. Загружаем полученные данные в память ускорителя
void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def); // (0)

	if(rank%2 == 0) // (1)
	{
		if (rank!=parts.x-1)
			right_send_recv(HostBuffer, rank+1, 500, locN, def); // (1.1)

		if (rank!=0)
			left_recv_send(HostBuffer, rank-1, 502, locN, def); // (1.2)
	}
	else
	{
		if (rank!=0) // В принципе, лишняя проверка
			left_recv_send(HostBuffer, rank-1, 500, locN, def); // (2.1)

		if (rank!=parts.x-1)
			right_send_recv(HostBuffer, rank+1, 502, locN, def); // (2.2)
	}

	save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def); // (3)
}

// Обмен граничными значениями давления P2, плотностей ro1 и ro2, Xi между процессорами
void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ro_g, DevArraysPtr.ro_g, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.Xi_g, DevArraysPtr.Xi_g, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.P_g, DevArraysPtr.P_g, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
#endif
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.ro_w, DevArraysPtr.ro_w, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.ro_n, DevArraysPtr.ro_n, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.Xi_w, DevArraysPtr.Xi_w, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.Xi_n, DevArraysPtr.Xi_n, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
}


// Обмен граничными значениями скоростей между процессорами
// В случае распределения расчетной области между процессорами по оси X
// передача u1y и u2y не требуется
void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	exchange(HostArraysPtr.ux_w, DevArraysPtr.ux_w, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.ux_n, DevArraysPtr.ux_n, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ux_g, DevArraysPtr.ux_g, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
#endif
}

// Обмен граничными значениями давления воды P1 и насыщенности NAPL S2 между процессорами
void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
#endif
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n, HostBuffer, DevBuffer, locN, blocksY, blocksZ, rank, parts, def);
}

void communication_initialization(int argc, char* argv[], int* size, int* rank, consts def)
{
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,size); // The amount of processors
    MPI_Comm_rank(MPI_COMM_WORLD,rank); // The number of processor
	//std::cout << "size=" <<*size<<"  "<<"rank= "<<*rank<<"\n";
}

void communication_finalization(void)
{
	MPI_Finalize();
}

// Реализация фунции Barrier для различных коммуникаций
void barrier(void)
{
	MPI_Barrier(MPI_COMM_WORLD);
}