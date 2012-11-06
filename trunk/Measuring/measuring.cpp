#include <mpi.h>
#include <float.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
#define isnan _isnan
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#define MEASURE_COUNT 10
#define MAX_BUFFER_SIZE 20000000
#define MIN_BUFFER_SIZE 10000000

void exchange(double* HostBuffer, int buffer_size, int size, int rank);
void right_send_recv(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id);
void left_recv_send(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id);

struct params
{
	clock_t task_time;
	int buffer_size;
};
typedef struct params params;

int main(int argc, char* argv[])
{
	double *HostBuffer, latency, send_double_time; 
	double sum_x, sum_xy, sum_x_2, sum_y; // переменные для метода наименьших квадратов
	int size, rank;
	params result[MEASURE_COUNT];
	char fname[] = "times.txt";
	FILE *fp;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 

	if (!(HostBuffer = new double[MAX_BUFFER_SIZE]))
		printf("Memory for *HostBuffer is not allocated in function host_memory_alloc");

	for(int i = 0; i < MAX_BUFFER_SIZE; i++)
	{
		HostBuffer[i] = rand();
	}

	if (!(fp = fopen(fname, "w")))
		printf("Not open file: %s", fname);
	fprintf(fp, "Buffer size\tExchange time, s\n");

	for(int i = 0; i < MEASURE_COUNT; i++)
	{
		result[i].buffer_size = MIN_BUFFER_SIZE + (int)((MAX_BUFFER_SIZE - MIN_BUFFER_SIZE) * i / (MEASURE_COUNT - 1));
		fprintf(fp, "%d \t", result[i].buffer_size);
		result[i].task_time = clock();
		exchange(HostBuffer, result[i].buffer_size, size, rank);
		MPI_Barrier(MPI_COMM_WORLD);
		result[i].task_time = clock() - result[i].task_time;
		fprintf(fp, "%.5f\n", (double)result[i].task_time / CLOCKS_PER_SEC);
	}

	fclose(fp);

	send_double_time = 0;
	latency = 0;
	sum_x = sum_y = sum_x_2 = sum_xy = 0;

	for(int i = 0; i < MEASURE_COUNT; i++)
	{
		sum_x += result[i].buffer_size;
		sum_y += result[i].task_time;
		sum_x_2 += result[i].buffer_size * result[i].buffer_size;
		sum_xy += result[i].buffer_size * result[i].task_time;
	}

	// Метод наименьших квадратов для нахождения коэффициентов прямой
	send_double_time = (sum_xy - sum_x * sum_y) / (sum_x_2 - sum_x * sum_x);
	latency = (sum_y - send_double_time * sum_x) / MEASURE_COUNT;

	// Переход от clock_t к секундам
	send_double_time /= CLOCKS_PER_SEC;
	latency /= CLOCKS_PER_SEC;

	if (!(rank))
	{
		printf("double_send_time = %.12f\tlatency = %.5f\n", send_double_time, latency);
	}

	delete HostBuffer;

	MPI_Finalize();

#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush(stdout);
	getchar();
#endif

	return 0;
}


void exchange(double* HostBuffer, int buffer_size, int size, int rank)
{
	if(size > 1) 
	{
		if ((rank) % 2 == 0)
		{
			if ((rank) != (size) - 1)
				right_send_recv(HostBuffer, buffer_size, (rank) + 1, 500);

			if ((rank) != 0)
				left_recv_send(HostBuffer, buffer_size, (rank) - 1, 502);
		}
		else
		{
			if ((rank) != 0) // В принципе, лишняя проверка
				left_recv_send(HostBuffer, buffer_size, (rank) - 1, 500); 

			if ((rank) != (size) - 1)
				right_send_recv(HostBuffer, buffer_size, (rank) + 1, 502);
		}
	}
}

// Передача и прием данных правой границе
void right_send_recv(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id, destination_rank, send_recv_id + 1, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}

// Получение и передача данных на левой границе
void left_recv_send(double* HostBuffer, int buffer_size, int destination_rank, int send_recv_id)
{
	MPI_Status status;

	if (! MPI_Sendrecv_replace(HostBuffer, buffer_size, MPI_DOUBLE, destination_rank, send_recv_id + 1, destination_rank, send_recv_id, MPI_COMM_WORLD, &status) == MPI_SUCCESS)
	{
		printf("MPI Error: MPI_Sendrecv_replace returned an error.\nFile:\"%s\"\nLine:\"%d\"\n\n", __FILE__, __LINE__);
	}
}