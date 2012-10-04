#include <mpi.h>
#include "defines.h"

/*
// Недописанная функция оптимального деления сетки по процессорам
void kak_delit(void)
{
	unsigned int Nx=250, Ny=240;
	unsigned int size=4;
	double t_l=0.001;
	double t_p=1;
	double t_r=2;
	double T=0;
	unsigned int N1=1;
	unsigned int N2=min(Nx, Ny);
	unsigned int N3=max(Nx, Ny);

	double T_max=2e50;
	unsigned int s2_max, s3_max;
	unsigned int flag=0;

	unsigned int s1=1;
	unsigned int s2, s3;
	for(s2=1;s2<size && s2<N2;s2++)
		for(s3=1;s3<size/s2 && s3<N3;s3++)
		{
			T=t_l + N2*t_p + t_r*(N1*N2*N3/(s1*s2*s3));
			if (T<T_max)
			{
				T_max=T;
				flag=1;
				s2_max=s2;
				s3_max=s3;
			}

			T=2*t_l + (N2 / s2 + N3 / s3)*t_p + t_r*(N1*N2*N3/(s1*s2*s3));
			if (T<T_max)
			{
				T_max=T;
				flag=2;
				s2_max=s2;
				s3_max=s3;
			}
		}

		std::cout<<"s2="<<s2<<"  s3="<<s3<<"  flag="<<flag<<"\n";
}
*/

// Передача и прием данных правой границе
void right_send_recv(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id, destination_rank, send_recv_id + 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}

// Получение и передача данных на левой границе
void left_recv_send(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id, consts def)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id + 1, destination_rank, send_recv_id, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
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
void exchange(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	if(def.sizex > 1) 
	{
		if ((def.rankx) % 2 == 0) // (1)
		{
			if ((def.rankx) != (def.sizex) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNy) * (def.locNz), (def.rank) + 1, 500, def);    // (1.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'r'); // (3)
			}

			if ((def.rankx) != 0)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNy) * (def.locNz), (def.rank) - 1, 502, def);    // (1.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'l'); // (3)
			}
		}
		else
		{
			if ((def.rankx) != 0) // В принципе, лишняя проверка
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNy) * (def.locNz), (def.rank) - 1, 500, def);    // (2.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'l'); // (3)
			}

			if ((def.rankx) != (def.sizex) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNy) * (def.locNz), (def.rank) + 1, 502, def);    // (2.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'x', 'r'); // (3)
			}
		}
	}
	if(def.sizey > 1) 
	{
		if ((def.ranky) % 2 == 0) // (1)
		{
			if ((def.ranky) != (def.sizey) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNx) * (def.locNz), (def.rank) + (def.sizex), 504, def);    // (1.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'r'); // (3)
			}

			if ((def.ranky) != 0)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNx) * (def.locNz), (def.rank) - (def.sizex), 506, def);    // (1.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'l'); // (3)
			}
		}
		else
		{
			if ((def.ranky) != 0) // В принципе, лишняя проверка
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNx) * (def.locNz), (def.rank) - (def.sizex), 504, def);    // (2.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'l'); // (3)
			}

			if ((def.ranky) != (def.sizey) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNx) * (def.locNz), (def.rank) + (def.sizex), 506, def);    // (2.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'y', 'r'); // (3)
			}
		}
	}
	if(def.sizez > 1) 
	{
		if ((def.rankz) % 2 == 0) // (1)
		{
			if ((def.rankz) != (def.sizez) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 508, def);    // (1.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'r'); // (3)
			}

			if ((def.rankz) != 0)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 510, def);    // (1.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'l'); // (3)
			}
		}
		else
		{
			if ((def.rankz) != 0) // В принципе, лишняя проверка
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'l'); // (0)
				left_recv_send(HostBuffer, (def.locNx) * (def.locNy), (def.rank) - (def.sizex) * (def.sizey), 500, def);    // (2.1)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'l'); // (3)
			}

			if ((def.rankz) != (def.sizez) - 1)
			{
				load_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'r'); // (0)
				right_send_recv(HostBuffer, (def.locNx) * (def.locNy), (def.rank) + (def.sizex) * (def.sizey), 502, def);    // (2.2)
				save_exchange_data(HostArrayPtr, DevArrayPtr, HostBuffer, DevBuffer, def, 'z', 'r'); // (3)
			}
		}
	}
}

// Обмен граничными значениями давления P2, плотностей ro1 и ro2, Xi между процессорами
void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ro_g, DevArraysPtr.ro_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_g, DevArraysPtr.Xi_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.P_g, DevArraysPtr.P_g, HostBuffer, DevBuffer, def);
#endif
	exchange(HostArraysPtr.P_n, DevArraysPtr.P_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ro_w, DevArraysPtr.ro_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ro_n, DevArraysPtr.ro_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_w, DevArraysPtr.Xi_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.Xi_n, DevArraysPtr.Xi_n, HostBuffer, DevBuffer, def);
}


// Обмен граничными значениями скоростей между процессорами
// В случае распределения расчетной области между процессорами по оси X
// передача u1y и u2y не требуется
void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
	exchange(HostArraysPtr.ux_w, DevArraysPtr.ux_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.ux_n, DevArraysPtr.ux_n, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.uy_w, DevArraysPtr.uy_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.uy_n, DevArraysPtr.uy_n, HostBuffer, DevBuffer, def);
#ifdef THREE_PHASE
	exchange(HostArraysPtr.ux_g, DevArraysPtr.ux_g, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.uy_g, DevArraysPtr.uy_g, HostBuffer, DevBuffer, def);
#endif

	if((def).Nz >= 2)
	{
		exchange(HostArraysPtr.uz_w, DevArraysPtr.uz_w, HostBuffer, DevBuffer, def);
		exchange(HostArraysPtr.uz_n, DevArraysPtr.uz_n, HostBuffer, DevBuffer, def);
#ifdef THREE_PHASE
		exchange(HostArraysPtr.uz_g, DevArraysPtr.uz_g, HostBuffer, DevBuffer, def);
#endif
	}
}

// Обмен граничными значениями давления воды P1 и насыщенности NAPL S2 между процессорами
void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def)
{
#ifdef THREE_PHASE
	exchange(HostArraysPtr.S_w, DevArraysPtr.S_w, HostBuffer, DevBuffer, def);
#endif
	exchange(HostArraysPtr.P_w, DevArraysPtr.P_w, HostBuffer, DevBuffer, def);
	exchange(HostArraysPtr.S_n, DevArraysPtr.S_n, HostBuffer, DevBuffer, def);
}

void communication_initialization(int argc, char* argv[], consts* def)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &((*def).size)); // The amount of processors
	MPI_Comm_rank(MPI_COMM_WORLD, &((*def).rank)); // The number of processor
	//std::cout << "size =" <<(*def).size<<"  "<<"rank = "<<(*def).rank<<"\n";
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

