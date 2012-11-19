#include "measuring.h"

struct params
{
	clock_t task_time;
	int buffer_size;
};
typedef struct params params;

int main(int argc, char* argv[])
{
	double *HostBuffer, *DevBuffer = 0, latency, send_double_time; 
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
		printf("Host-host: double_send_time = %e\tlatency = %.5f\n", send_double_time, latency);
	}

#ifdef USE_GPU
	if (!(rank))
	{
		measuring_host_device_exchange(HostBuffer, DevBuffer, rank);
	}
#endif

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

void measuring_host_device_exchange(double *HostBuffer, double *DevBuffer, int rank) 
{
	params host_device_result[MEASURE_COUNT];
	double sum_x, sum_xy, sum_x_2, sum_y; // переменные для метода наименьших квадратов
	double latency, send_double_time; 

	device_initialization(rank);
	device_memory_allocation(&DevBuffer, MAX_BUFFER_SIZE);
	set_devbuffer_values(DevBuffer, MAX_BUFFER_SIZE);

	for(int i = 0; i < MEASURE_COUNT; i++)
	{
		host_device_result[i].buffer_size = MIN_BUFFER_SIZE + (int)((MAX_BUFFER_SIZE - MIN_BUFFER_SIZE) * i / (MEASURE_COUNT - 1));
		host_device_result[i].task_time = clock();
		load_exchange_data_part(HostBuffer, DevBuffer, host_device_result[i].buffer_size);
		save_exchange_data_part(HostBuffer, DevBuffer, host_device_result[i].buffer_size);
		host_device_result[i].task_time = clock() - host_device_result[i].task_time;
	}

	send_double_time = 0;
	latency = 0;
	sum_x = sum_y = sum_x_2 = sum_xy = 0;

	for(int i = 0; i < MEASURE_COUNT; i++)
	{
		sum_x += host_device_result[i].buffer_size;
		sum_y += host_device_result[i].task_time;
		sum_x_2 += host_device_result[i].buffer_size * host_device_result[i].buffer_size;
		sum_xy += host_device_result[i].buffer_size * host_device_result[i].task_time;
	}

	// Метод наименьших квадратов для нахождения коэффициентов прямой
	send_double_time = (sum_xy - sum_x * sum_y) / (sum_x_2 - sum_x * sum_x);
	latency = (sum_y - send_double_time * sum_x) / MEASURE_COUNT;

	// Переход от clock_t к секундам
	send_double_time /= CLOCKS_PER_SEC;
	latency /= CLOCKS_PER_SEC;

	if (!(rank))
	{
		printf("Host-device: double_send_time = %e\tlatency = %.5f\n", send_double_time, latency);
	}

	device_memory_free(DevBuffer);
	device_finalization();
}