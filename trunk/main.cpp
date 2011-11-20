#include "defines.h"

// Счетчики времени
clock_t start_time, finish_time;  
// Буферные массивы для обмена между процессорами
double *HostBuffer;
double *DevBuffer;

int main(int argc, char* argv[])
{
	consts def;
	read_defines(argc, argv, &def);

	// Переменные количества процессоров и номера текущего процессора
	int size=0, rank=0;
	// Хостовый массив данных расчетной области процессора
	ptr_Arrays HostArraysPtr;
	// GPU-массив данных расчетной области процессора
	ptr_Arrays DevArraysPtr;
	// Размеры (локальные) расчетной области процессора
	localN locN;
	// Количество блоков ускорителя
	int blocksX=0, blocksY=0, blocksZ=0;
	// Счетчик шагов по времени
	int j=0;
	// Количество частей, на которые делим область
	parts_sizes parts;
	
	// Инициализация коммуникаций, перевод глобальных параметров в локальные процессора, 
	// выделение памяти, загрузка начальных/сохраненных данных
	initialization(&HostArraysPtr, &DevArraysPtr, &j, &locN, &size, &parts, &rank, &blocksX, &blocksY, &blocksZ, argc, argv, def);

	// Нулевой процессор выводит название запускаемой задачи
	if (!rank)
	{
#ifdef TWO_PHASE
		std::cout << "Two phase filtration by CAPAZ on "<<size<<" node(s).\n\n";
#endif
#ifdef THREE_PHASE
		std::cout << "Three phase filtration by CAPAZ on "<<size<<" node(s).\n\n";
#endif
#ifdef B_L
		std::cout << "Backley-Leverett filtration by CAPAZ on "<<size<<" node(s).\n\n";
#endif
	}

	// Тест
	//save_data_plots(HostArraysPtr, DevArraysPtr, 0, size, rank, localNx, def);
	
	start_time=clock();

	// Цикл шагов по времени (каждая итерация - новый слой по времени)
	// 1. Проводятся расчеты P1 и S2 на следующем временном слое
	// 2. Каждые (def.print_screen) раз на экран выводится информация о временном слое
	// 3. Каждые save_plots раз данные выгружаются в память хоста и 
	//    сохраняются в файлы графиков (**), в файл сохраняется состояние задачи (***)
	for (j++; j <= def.timeX/(def.dt); j++)
	{
		if ((j % (def.print_screen) == 0) && rank==0) // (2)
		{
			printf ("t=%.3f\n",j*(def.dt)); 
			fflush(stdout);
		}

		time_step_function(HostArraysPtr, DevArraysPtr, DevBuffer, def,j*(def.dt),locN,rank,parts,blocksX,blocksY,blocksZ); // (1)

		if ((j % (def.save_plots)) == 0) // (3)
		{
			// Следующие 2 функции вызываются строго в таком порядке,
			// т.к. save использует данные, загруженные save_data_plots
			save_data_plots(HostArraysPtr, DevArraysPtr, j*(def.dt), parts, rank, locN, def); // (**)
			//save(HostArraysPtr, DevArraysPtr, j, rank, size, localNx); // (***)
		}
	}

	// Вывод информации о времени работы программы в секундах
	finish_time=clock();
	finish_time-=start_time;
	printf( "Task time in seconds:\t%.2f\n", (double) finish_time/CLOCKS_PER_SEC);

	// Завершение работы и освобождение памяти
	finalization(HostArraysPtr, DevArraysPtr, DevBuffer);

	// При запуске в Windows после работы программы оставлять окно консоли
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush( stdout);
	getchar();
#endif
	return 0;
}

// Функция полного цикла расчетов на следующем временном слое
// 1. Расчет плотности жидкостей ro, давлений NAPL P2, переменных Xi
// 2. Обмен между процессорами пограничными значениями P2, ro и Xi
// 3. Расчет скоростей жидкостей
// 4. Обмен между процессорами пограничными значениями скоростей жидкостей
// 5. Расчет переменной roS на следующем временном слое
// 6. Расчет методом Ньютона давления воды P1 и насыщенности DNAPL S2
// 7. Применение граничных условий для P1 и S2
// 8. Обмен между процессорами пограничными значениями P1 и S2
void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t, localN locN, int rank,parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	P_S_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, locN,blocksY, blocksZ, rank, parts, def); // (8)
	ro_P_Xi_calculation(HostArraysPtr,DevArraysPtr,def,locN,rank,parts,blocksX,blocksY, blocksZ); // (1)
	P_ro_Xi_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, locN,blocksY, blocksZ, rank, parts, def); // (2)
	u_calculation(HostArraysPtr,DevArraysPtr,locN,rank,parts,blocksX,blocksY, blocksZ, def); // (3)
	u_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, locN,blocksY, blocksZ, rank, parts, def); // (4)
	roS_calculation(HostArraysPtr,DevArraysPtr,def,t,locN,rank,parts,blocksX,blocksY, blocksZ); // (5)
	P_S_calculation(HostArraysPtr,DevArraysPtr,def,locN,rank,parts,blocksX,blocksY, blocksZ); // (6)
	boundary_conditions(HostArraysPtr,DevArraysPtr,locN,rank,parts,blocksX,blocksY, blocksZ, def); // (7)
	
}

// Преобразование локальных координат процессора к глобальным
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа 
// (если 2 соседа с обеих сторон,то +2 точки). 
// Глобальные границы хранятся как обычные точки (отсюда и условие на rank==0)
int i_to_I(int i, int rank, parts_sizes parts, consts def)
{
	int I;
	if (rank <= (def.Nx)%parts.x)
	{
		if(rank==0)
			I=i;
		else
			I=((def.Nx)/parts.x+1)*rank+i-1;
	}
	else
		I=((def.Nx)/parts.x+1)*rank-(rank-(def.Nx)%parts.x)+i-1;

	test_positive(I, __FILE__, __LINE__);
	return I;
}

// Вычисление расчетной области (нагрузки) процессора
// Если поровну не распределяется, то первые (def.NX)%size получают +1 в нагрузку.
// Каждый процессор содержит дополнительную точку в массиве для
// обмена данными, если имеет соседа 
// (если 2 соседа с обеих сторон,то +2 точки). 
// Глобальные границы хранятся как обычные точки (отсюда и условие на rank==0)
void global_to_local_vars (localN* locN, parts_sizes parts, int rank, consts def)
{
	(*locN).x = (def.Nx) / parts.x;

	if (rank % parts.x < (def.Nx) % parts.x)
		((*locN).x) ++;

	// Крайние процессоры получают по 1 точке для граничных данных,
	// остальные - по 2 на обе границы
	// Если процессор один, то границ у него нет и дополнительные точки не нужны
	if (parts.x > 1)
	{
		if ((rank % parts.x == 0) || (rank % parts.x  == parts.x - 1))
			((*locN).x) ++;
		else
			((*locN).x) += 2;
	}

	(*locN).y = (def.Ny) / parts.y;

	if ((rank % ((parts.x) * (parts.y))) / parts.x < (def.Ny) % parts.y)
		((*locN).y) ++;

	if (parts.y > 1)
	{
		if (((rank % ((parts.x) * (parts.y))) / parts.x == 0) || ((rank % ((parts.x) * (parts.y))) / parts.x == parts.y - 1))
			((*locN).y) ++;
		else
			((*locN).y) += 2;
	}

	(*locN).z = (def.Nz) / parts.z;

	if (rank / (parts.x) / (parts.y) < (def.Nz) % parts.z)
		((*locN).z) ++;

	if (parts.z > 1)
	{
		if ((rank / (parts.x) / (parts.y) == 0) || (rank / (parts.x) / (parts.y) == parts.z - 1))
			((*locN).z) ++;
		else
			((*locN).z) += 2;
	}

	test_positive((*locN).x, __FILE__, __LINE__);
	test_positive((*locN).y, __FILE__, __LINE__);
	test_positive((*locN).z, __FILE__, __LINE__);
}

// Является ли точка активной (т.е. не предназначенной только для обмена на границах)
int is_active_point(int i, int j, int k, localN locN, int rank, parts_sizes parts)
{
	if((rank % parts.x != 0 && i == 0) || (rank % parts.x != parts.x - 1 && i == locN.x - 1)
		|| ((rank % ((parts.x) * (parts.y))) / parts.x != 0 && j == 0)	|| ((rank % ((parts.x) * (parts.y))) / parts.x != parts.y - 1 && j == locN.y - 1)
		|| (((rank / (parts.x) / (parts.y) != 0 && k == 0) || (rank / (parts.x) / (parts.y) == parts.z - 1 && k == locN.z - 1)) && parts.z > 1))
		return 0;
	else
		return 1;
}

// Применение начальных данных во всех точках
void parts_initialization(int size, parts_sizes *parts)
{
	(*parts).x = size;
	(*parts).y = 1;
	(*parts).z = 1;
}

void data_initialization(ptr_Arrays HostArraysPtr, int* t, localN locN, int rank, parts_sizes parts, consts def)
{
	*t=0;
	for(int i=0;i<locN.x;i++)
		for(int j=0;j<locN.y;j++)
			for(int k=0;k<locN.z;k++)
				if(is_active_point(i, j, k, locN, rank, parts))
					{
						// Преобразование локальных координат процессора к глобальным
						int I=i_to_I(i,rank,parts,def);

#ifdef THREE_PHASE
						int j1 = locN.y / 2;

						if(j < j1)
						{
							HostArraysPtr.S_w[i+j*locN.x+k*locN.x*locN.y] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
							HostArraysPtr.S_g[i+j*locN.x+k*locN.x*locN.y] = def.S_g_gr + (def.S_g_init - def.S_g_gr) * j / j1;
						}
						else
						{
							HostArraysPtr.S_w[i+j*locN.x+k*locN.x*locN.y] = def.S_w_init;
							HostArraysPtr.S_g[i+j*locN.x+k*locN.x*locN.y] = def.S_g_init;
						}

						if(j == 0)
							HostArraysPtr.P_n[i+j*locN.x+k*locN.x*locN.y] = def.P_atm;
						else
							HostArraysPtr.P_n[i+j*locN.x+k*locN.x*locN.y] = HostArraysPtr.P_n[i + (j - 1) * locN.x + k * locN.x * locN.y]
								+ (def.ro0_n * (1. - HostArraysPtr.S_w[i + (j - 1) * locN.x + k * locN.x * locN.y] - HostArraysPtr.S_g[i + (j - 1) * locN.x + k * locN.x * locN.y]) 
								+ def.ro0_w * HostArraysPtr.S_w[i + (j - 1) * locN.x + k * locN.x * locN.y]
								+ def.ro0_g * HostArraysPtr.S_g[i + (j - 1) * locN.x + k * locN.x * locN.y]) * (def.g_const) * (def.hy);

						HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y] = 0;						
#else
						// Если точка на верхней границе, не далее (def.source) точек от центра,
						// то в ней начальная насыщенность. Иначе, нулевая
						if ((j==0) && (I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
							HostArraysPtr.S_n[i+j*locN.x+k*locN.x*locN.y]=def.S_n_gr;
						else
							HostArraysPtr.S_n[i+j*locN.x+k*locN.x*locN.y]=0;

						HostArraysPtr.P_w[i+j*locN.x+k*locN.x*locN.y]=def.P_atm+j * (def.ro0_w) * (def.g_const)*(def.hy);
						HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y]=0;
#endif
					
						/*
						if ((HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]>=(def.NX)/2.*(def.h1)) && (HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]<=4.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]<=2./5.*locN.y*(def.h2)) && (HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]>=(-1.)*HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]/4.+2./5.*locN.y*(def.h2)))
								HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y]=1;

						if ((HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]>=(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]<=2.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]<=4./5.*locN.y*(def.h2)) && (HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]>=3./5.*locN.y*(def.h2)))
								HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y]=1;
								*/
					
						/*
						if ((HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]>=2.*(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*locN.x+k*locN.x*locN.y]<=3.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]>=1./10.*locN.y*(def.h2)) && (HostArraysPtr.y[i+j*locN.x+k*locN.x*locN.y]<=3./10.*locN.y*(def.h2)))
								HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y]=1;
						*/
#ifndef THREE_PHASE 
					test_nan(HostArraysPtr.S_n[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
#else
					test_nan(HostArraysPtr.P_n[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_g[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_w[i+j*locN.x+k*locN.x*locN.y], __FILE__, __LINE__);
#endif
					}
}

//----------------------------------------------------------------------------------------------------
// Служебные функции

// Инициализация коммуникаций (1), перевод глобальных параметров в локальные процессора (2), 
// инициализация ускорителя (2.5), выделение памяти (3), загрузка начальных/сохраненных данных (4)
// Для задачи Б-Л загрузка проницаемостей из файла.
void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int* j, localN* locN, int* size, parts_sizes *parts, int* rank, int* blocksX, int* blocksY, int* blocksZ, int argc, char* argv[], consts def)
{
	FILE *f_save;

	communication_initialization(argc, argv, size, rank, def); // (1)

	parts_initialization(*size, parts);

	global_to_local_vars(locN, *parts, *rank, def); // (2)

	device_initialization(*rank, blocksX, blocksY, blocksZ, *locN, def); // (2.5)

	memory_allocation(HostArraysPtr, DevArraysPtr, *locN, def); // (3)

#ifdef B_L
	load_permeability((*HostArraysPtr).K, *locN); // (5)
	load_data_to_device((*HostArraysPtr).K, (*DevArraysPtr).K, *locN, def);
#endif

	// Если процессор может открыть файл сохраненного состояния,
	// то восстанавливаем состояние, иначе применяем начальные условия
	if (f_save=fopen("save/save.dat","rb"))
	{
		fclose(f_save);
		restore(*HostArraysPtr, j, *rank, *parts, *locN, def);
	}
	else
		data_initialization (*HostArraysPtr, j, *locN, *rank, *parts, def); // (4)

#ifdef THREE_PHASE
	load_data_to_device((*HostArraysPtr).P_n, (*DevArraysPtr).P_n, *locN, def);
	load_data_to_device((*HostArraysPtr).S_w, (*DevArraysPtr).S_w, *locN, def);
	load_data_to_device((*HostArraysPtr).S_g, (*DevArraysPtr).S_g, *locN, def);
	load_data_to_device((*HostArraysPtr).roS_g_old, (*DevArraysPtr).roS_g_old, *locN, def);
#else
	load_data_to_device((*HostArraysPtr).P_w, (*DevArraysPtr).P_w, *locN, def);
	load_data_to_device((*HostArraysPtr).S_n, (*DevArraysPtr).S_n, *locN, def);
#endif
	load_data_to_device((*HostArraysPtr).roS_w_old, (*DevArraysPtr).roS_w_old, *locN, def);
	load_data_to_device((*HostArraysPtr).roS_n_old, (*DevArraysPtr).roS_n_old, *locN, def);
	load_data_to_device_int((*HostArraysPtr).media, (*DevArraysPtr).media, *locN, def);
}

// Завершение работы (1), освобождение памяти (2)
void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	memory_free(HostArraysPtr, DevArraysPtr); // (2)
	communication_finalization(); // (1)
	device__finalization(); // (1)
}

// Выделение памяти хоста (1) и ускорителя (2) под массив точек расчетной области
void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, localN locN, consts def)
{
	host_memory_allocation(HostArraysPtr, locN, def); // (1)
	device_memory_allocation(DevArraysPtr, &DevBuffer, locN, def); // (2)
}

// Освобожение памяти хоста (1) и ускорителя (2) из под массива точек расчетной области
void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr)
{
	host_memory_free(HostArraysPtr); // (1)
	device_memory_free(DevArraysPtr, DevBuffer); // (2)
}

// Выделение памяти хоста под массив точек расчетной области
void host_memory_allocation(ptr_Arrays* ArraysPtr, localN locN, consts def)		
{	
	if (!(HostBuffer=new double[2 * ((locN.y) * (locN.z))]))
		printf ("\nWarning! Memory for *HostBuffer is not allocated in function host_memory_alloc\n");

	try
	{
		(*ArraysPtr).P_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).P_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ro_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ro_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ux_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uy_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uz_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ux_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uy_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uz_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).Xi_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).Xi_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_w_old=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_n=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_n_old=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).media=new int [(locN.x)*(locN.y)*(locN.z)];
#ifdef B_L
		(*ArraysPtr).K=new double [(locN.x)*(locN.y)*(locN.z)];
#endif
#ifdef THREE_PHASE
		(*ArraysPtr).P_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).S_w=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).S_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ro_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).ux_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uy_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).uz_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).Xi_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_g=new double [(locN.x)*(locN.y)*(locN.z)];
		(*ArraysPtr).roS_g_old=new double [(locN.x)*(locN.y)*(locN.z)];
#else
		(*ArraysPtr).S_n=new double [(locN.x)*(locN.y)*(locN.z)];
#endif
	}
	catch(...)
	{
		printf ("\nError! Not enough host memory\n");
		exit(0);
	}
}

// Освобожение памяти хоста из под массива точек расчетной области
void host_memory_free(ptr_Arrays ArraysPtr)
{ 
	delete HostBuffer;
	delete[] ArraysPtr.P_w;
	delete[] ArraysPtr.P_n;
	delete[] ArraysPtr.ro_w;
	delete[] ArraysPtr.ro_n;
	delete[] ArraysPtr.ux_w;
	delete[] ArraysPtr.uy_w;
	delete[] ArraysPtr.uz_w;
	delete[] ArraysPtr.ux_n;
	delete[] ArraysPtr.uy_n;
	delete[] ArraysPtr.uz_n;
	delete[] ArraysPtr.Xi_w;
	delete[] ArraysPtr.Xi_n;
	delete[] ArraysPtr.roS_w;
	delete[] ArraysPtr.roS_w_old;
	delete[] ArraysPtr.roS_n;
	delete[] ArraysPtr.roS_n_old;
	delete[] ArraysPtr.media;
#ifdef B_L
	delete[] ArraysPtr.K;
#endif
#ifdef THREE_PHASE
	delete[] ArraysPtr.P_g;
	delete[] ArraysPtr.S_w;
	delete[] ArraysPtr.S_g;
	delete[] ArraysPtr.ro_g;
	delete[] ArraysPtr.ux_g;
	delete[] ArraysPtr.uy_g;
	delete[] ArraysPtr.uz_g;
	delete[] ArraysPtr.Xi_g;
	delete[] ArraysPtr.roS_g;
	delete[] ArraysPtr.roS_g_old;
#else
	delete[] ArraysPtr.S_n;
#endif
}

// Функция сохранения графиков в файлы
void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, parts_sizes parts, int rank, localN locN, consts def)
{
	// Загрузка в память хоста результатов расчета
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.P_n, DevArraysPtr.P_n , locN, def);
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtr.S_w , locN, def);
	load_data_to_host(HostArraysPtr.S_g, DevArraysPtr.S_g , locN, def);
	load_data_to_host(HostArraysPtr.ux_w, DevArraysPtr.ux_w , locN, def);
	load_data_to_host(HostArraysPtr.uy_w, DevArraysPtr.uy_w , locN, def);
	load_data_to_host(HostArraysPtr.uz_w, DevArraysPtr.uz_w , locN, def);
	load_data_to_host(HostArraysPtr.ux_g, DevArraysPtr.ux_g , locN, def);
	load_data_to_host(HostArraysPtr.uy_g, DevArraysPtr.uy_g , locN, def);
	load_data_to_host(HostArraysPtr.uz_g, DevArraysPtr.uz_g , locN, def);
#else
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtr.P_w , locN, def);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtr.S_n , locN, def);
#endif
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtr.ux_n , locN, def);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtr.uy_n , locN, def);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtr.uz_n , locN, def);

#ifndef THREE_PHASE
	// Проверка на выход из допустимого диапазона значений P и S
#ifdef MY_TEST
	test_correct_P_S(HostArraysPtr, locN, rank, def);
#endif
#endif
	
	// Нулевой процессор создает директории, файлы и прописывает заголовки файлов
	if (rank==0)
		print_plots_top (t, def);

	// По очереди для каждого из процессоров вызываем функцию вывода на график
	// своей части массива.
	for (int cpu = 0; cpu < (parts.x * parts.y * parts.z); cpu ++)
	{
		// Реализация фунции Barrier для различных коммуникаций
		barrier();
		if (rank==cpu)
			print_plots(HostArraysPtr, t, rank, parts, locN,def);	
	}
}

// Функция создания директорий, файлов для графиков и сохранения заголовков в них 
void print_plots_top (double t, consts def)
{
	char fname[30];
	FILE *fp;

	sprintf(fname,"plots/S=%012.4f.dat",t);

#ifdef _WIN32
	_mkdir("plots");
#else
	mkdir("plots",0000777);
#endif

	// Создание (или перезапись) файла с графиками
	// 1. Для распределения насыщенностей NAPL S_n
	// 2. Для распределения давлений воды P_w
	// 3. Для распределения скоростей {u_x, u_y, u_z}
	// 4. Для распределения типов грунтов
	if(!(fp=fopen(fname,"wt")))
		std::cout << "Not open file(s) in function SAVE_DATA_PLOTS! \n";

	fprintf(fp,"TITLE =  \"Filtration in time=%5.2f\" \n", t); 

	if((def.Nz) < 2)
	{
#ifdef THREE_PHASE
//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_w\",\"S_g\",\"S_n\",\"P_n\",\"u_x\",\"u_y\",\"media\" \n");
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_w\",\"S_g\",\"S_n\",\"P_n\",\"uw_x\",\"uw_y\",\"ug_x\",\"ug_y\",\"un_x\",\"un_y\",\"media\" \n");
#else
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_n\",\"P_w\",\"u_x\", \"u_y\",\"media\" \n");
#endif
		fprintf(fp,"ZONE T = \"BIG ZONE\", K=%d,J=%d, F = POINT\n", (def.Nx), (def.Ny));
	}
	else
	{
#ifdef THREE_PHASE
//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_g\",\"S_n\",\"P_n\",\"u_x\", \"u_y\",\"u_z\",\"media\" \n");
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_g\",\"S_n\",\"P_n\",\"uw_x\",\"uw_y\",\"uw_z\",\"ug_x\",\"ug_y\",\"ug_z\",\"un_x\",\"un_y\",\"un_z\",\"media\" \n");
#else
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_n\",\"P_w\",\"u_x\", \"u_y\", \"u_z\", \"media\" \n");
#endif
		fprintf(fp,"ZONE T = \"BIG ZONE\", K=%d,J=%d,I=%d, F = POINT\n", (def.Nx), (def.Ny), (def.Nz));
	}

	fclose(fp);
}


// Функция сохранения данных в файлы графиков (!3D)
void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, parts_sizes parts, localN locN, consts def)
{
	char fname[30];
	FILE *fp;
	int local;

	sprintf(fname,"plots/S=%012.4f.dat",t);
		
	// Открытие на дозапись и сохранение графиков
	// 1. Для распределения насыщенностей NAPL S_n
	// 2. Для распределения давлений воды P_w
	// 3. Для распределения скоростей {u_x, u_y, u_z}
	// 4. Для распределения типов грунтов
	if(!(fp=fopen(fname,"at")))
		std::cout << "Not open file(s) in function SAVE_DATA_PLOTS! \n";

	for(int i=0; i<locN.x; i++)
		for(int j=0; j<locN.y; j++)
			for(int k=0; k<locN.z; k++)
				if(is_active_point(i, j, k, locN, rank, parts))
				{
					local=i+j*locN.x+k*locN.x*locN.y;

					// Преобразование локальных координат процессора к глобальным
					int I=i_to_I(i,rank,parts,def);
#ifdef THREE_PHASE
					if(def.Nz < 2)
					{
/*						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/					
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.ux_n[local], 
							(-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
		
					}

					else
					{
/*						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
							HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local],
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);

					}
#else
					if(def.Nz < 2)
					{
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)

					}
					else
					{
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)
					}
#endif
				}

	/* Не очень хорошо, так как запуск функции, условный оператор в цикле
	for(int i=0; i<locN.x; i++)
		for(int k=0; k<locN.z; k++)
			if ((i_to_I(i,rank,size, def)==(def.Nx)/2) && is_active_point(i,locN.x,rank,size))
				for(int j=0; j<locN.y; j++)
					fprintf(fp_S2y,"%.2e %.3e\n", HostArraysPtr.y[locN.x/2+j*locN.x+k*locN.x*locN.y], HostArraysPtr.S_n[locN.x/2+j*locN.x+k*locN.x*locN.y]); 
	

	for(int i=0; i<locN.x; i++)
		for(int k=0; k<locN.z; k++)
			fprintf(fp_S2x,"%.2e %.3e\n", I*(def.hx), HostArraysPtr.S_n[i+locN.x*locN.y/2+k*locN.x*locN.y]); 
	*/
		
	fclose(fp);
}

// Функция загрузки файла проницаемостей
void load_permeability(double* K, localN locN)
{
	FILE *input;
	char *file="../noise.dat";

	if(!(input=fopen(file,"rt")))
	{
		file="noise.dat";
		if(!(input=fopen(file,"rt")))
		{
			printf("Not open file \"%s\"!\nError in file \"%s\" at line %d\n", file,__FILE__,__LINE__);
			fflush(stdout);
		}
	}

	int Nx, Ny;
	fscanf(input, "%d %d\n", &Nx, &Ny);

	if((Nx!=locN.x) || (Ny!=locN.y))
	{
		printf("Nx/Ny from noise.dat not equal\nError in file \"%s\" at line %d\n",__FILE__,__LINE__);
		fflush(stdout);
	}

	char* str=new char[30*Nx];

	char value[30];
	for(int j=0; j<Ny; j++)
	{
		int n=0;
		fgets(str,30*Nx,input);
		for(int i=0; i<Nx; i++)			
		{
			int iter=0;
			if(str[n]==' ')
				n++;
			for (n;str[n]!=' ';n++,iter++)
			{
				value[iter]=str[n];
			}
			value[iter]='\0';
			n++;

			for (int k=0;k<locN.z;k++)
				K[i+j*locN.x+k*locN.x*locN.y]=1e-10 * exp(atof(value));
		}
	}

	fclose(input);
}

// Сохранение состояния в файл
void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, int rank, parts_sizes parts, localN locN, consts def)
{
	// Загружаем в память хоста данные по roS_old
	// P1 и S2 загружены уже для функции сохранения графиков,
	// x,y и media не изменились и загрузки не требуют.
	//load_data_to_host(HostArraysPtr.P1, DevArraysPtr.P1 , localNx);
	//load_data_to_host(HostArraysPtr.S2, DevArraysPtr.S2 , localNx);
	load_data_to_host(HostArraysPtr.roS_w_old, DevArraysPtr.roS_w_old , locN, def);
	load_data_to_host(HostArraysPtr.roS_n_old, DevArraysPtr.roS_n_old , locN, def);
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.roS_g_old, DevArraysPtr.roS_n_old , locN, def);
#endif

	FILE *f_save;

	if (rank==0)
	{

#ifdef _WIN32
		_mkdir("save");
#else
		mkdir("save",0000777);
#endif
	
		if(!(f_save=fopen("save/save.dat","wb")))
		{
			printf("\nError: Not open file \"save.dat\"!\n");
			exit(0);
		}
		fclose(f_save);
	}

	for (int cpu=0; cpu<(parts.x * parts.y * parts.z);cpu++)
	{
		// Реализация фунции Barrier для различных коммуникаций
		barrier();
		if (rank==cpu)
		{
			if(!(f_save=fopen("save/save.dat","ab")))
			{
				printf("\nError: Not open file \"save.dat\"!\n");
				exit(0);
			}
			fwrite(&j, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.P_n, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			fwrite(HostArraysPtr.S_w, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			fwrite(HostArraysPtr.S_g, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
#else
			fwrite(HostArraysPtr.P_w, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			fwrite(HostArraysPtr.S_n, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
#endif
			//fwrite(HostArraysPtr.x, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			//fwrite(HostArraysPtr.y, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			//fwrite(HostArraysPtr.z, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			fwrite(HostArraysPtr.roS_w_old, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
			fwrite(HostArraysPtr.roS_n_old, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.roS_g_old, sizeof(double), (locN.x) * (locN.y) * (locN.z), f_save);
#endif
			fwrite(HostArraysPtr.media, sizeof(int), (locN.x) * (locN.y) * (locN.z), f_save);
			fclose(f_save);
		}
	}
}

// Восстановление состояния из файла
void restore (ptr_Arrays HostArraysPtr, int* j, int rank, parts_sizes parts, localN locN, consts def)
{
	FILE *f_save;
	for (int cpu=0; cpu<(parts.x * parts.y * parts.z);cpu++)
	{
		// Реализация фунции Barrier для различных коммуникаций
		barrier();
		localN lN;

		lN.x = 0;
		lN.y = 0;
		lN.z = 0;

		if (rank==cpu)
		{
			if(!(f_save=fopen("save/save.dat","rb")))
			{
					printf("\nError: Not open file \"save.dat\"!\n");
					exit(0);
			}
			for (int queue=0;queue<=rank;queue++)
			{
				global_to_local_vars(&lN, parts, queue, def);
				fread(j, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.P_n, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				fread(HostArraysPtr.S_w, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				fread(HostArraysPtr.S_g, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
#else
				fread(HostArraysPtr.P_w, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				fread(HostArraysPtr.S_n, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
#endif
				//fread(HostArraysPtr.x, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				//fread(HostArraysPtr.y, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				//fread(HostArraysPtr.z, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				fread(HostArraysPtr.roS_w_old, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
				fread(HostArraysPtr.roS_n_old, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.roS_g_old, sizeof(double), (lN.x) * (lN.y) * (lN.y), f_save);
#endif
				fread(HostArraysPtr.media, sizeof(int), (lN.x) * (lN.y) * (lN.y), f_save);
			}
			fclose(f_save);
		}
	}
}



//------------------------------------------------------------------------------------------

// Считывание параметров задачи из файла
void read_defines(int argc, char *argv[], consts* def)
{
#ifdef THREE_PHASE
	(*def).aw[0]=aw[0];
	(*def).aw[1]=aw[1];
	(*def).bw[0]=bw[0];
	(*def).bw[1]=bw[1];
	(*def).ag[0]=ag[0];
	(*def).ag[1]=ag[1];
	(*def).bg[0]=bg[0];
	(*def).bg[1]=bg[1];
	(*def).S_w_range[0]=S_w_range[0];
	(*def).S_w_range[1]=S_w_range[1];
	(*def).S_w_range[2]=S_w_range[2];
	(*def).S_g_range[0]=S_g_range[0];
	(*def).S_g_range[1]=S_g_range[1];
	(*def).S_g_range[2]=S_g_range[2];
#endif

	FILE *defs;
	char *file;
	char str[250]="", attr_name[50]="", attr_value[50]="";

	file=DEFINES_FILE;

	if(!(defs=fopen(file,"rt")))
	{
		file="defines.ini";
		if(!(defs=fopen(file,"rt")))
		{
			printf("Not open file \"%s\"!\nError in file \"%s\" at line %d\n", file,__FILE__,__LINE__);
			fflush(stdout);
			return;
		}
	}

	while (!feof(defs))
	{
		unsigned int i,j,a;
		fgets(str,250,defs);
		if(str[0]=='#')
			continue;
		for (i=0;str[i]!='=';i++)
		{
			if (i>=strlen(str))
				continue;
			attr_name[i]=str[i];
		}

		attr_name[i]='\0';
		a=strlen(str);
		for(j=i+1;str[j]!=' ' && (j < a);j++)
			attr_value[j-i-1]=str[j];
		attr_value[j-i-1]='\0';

		if(!strcmp(attr_name,"HX")) 
			{(*def).hx = atof(attr_value); continue;}
		if(!strcmp(attr_name,"HY")) 
			{(*def).hy = atof(attr_value); continue;}
		if(!strcmp(attr_name,"HZ")) 
			{(*def).hz = atof(attr_value); continue;}
		if(!strcmp(attr_name,"TAU")) 
			{(*def).tau = atof(attr_value); continue;}
		if(!strcmp(attr_name,"DT")) 
			{(*def).dt = atof(attr_value); continue;}
		if(!strcmp(attr_name,"C_W")) 
			{(*def).c_w = atof(attr_value); continue;}
		if(!strcmp(attr_name,"C_N")) 
			{(*def).c_n = atof(attr_value); continue;}
		if(!strcmp(attr_name,"L")) 
			{(*def).l = atof(attr_value); continue;}
		if(!strcmp(attr_name,"BETA_W")) 
			{(*def).beta_w = atof(attr_value); continue;}
		if(!strcmp(attr_name,"BETA_N")) 
			{(*def).beta_n = atof(attr_value); continue;}
		if(!strcmp(attr_name,"RO_W")) 
			{(*def).ro0_w = atof(attr_value); continue;}
		if(!strcmp(attr_name,"RO_N")) 
			{(*def).ro0_n = atof(attr_value); continue;}
		if(!strcmp(attr_name,"MU_W")) 
			{(*def).mu_w = atof(attr_value); continue;}
		if(!strcmp(attr_name,"MU_N")) 
			{(*def).mu_n = atof(attr_value); continue;}
		if(!strcmp(attr_name,"G_CONST")) 
			{(*def).g_const = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_ATM")) 
			{(*def).P_atm = atof(attr_value); continue;}

		if(!strcmp(attr_name,"LAMBDA_0")) 
			{(*def).lambda[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"LAMBDA_1")) 
			{(*def).lambda[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"M_0")) 
			{(*def).m[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"M_1")) 
			{(*def).m[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_WR_0")) 
			{(*def).S_wr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_WR_1")) 
			{(*def).S_wr[1] = atof(attr_value); continue;}

#ifndef B_L
		if(!strcmp(attr_name,"K_0")) 
			{(*def).K[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"K_1")) 
			{(*def).K[1] = atof(attr_value); continue;}
#endif
#ifdef TWO_PHASE
		if(!strcmp(attr_name,"P_D_0")) 
			{(*def).P_d[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_D_1")) 
			{(*def).P_d[1] = atof(attr_value); continue;}
#endif

#ifdef THREE_PHASE
		if(!strcmp(attr_name,"C_G")) 
			{(*def).c_g = atof(attr_value); continue;}
		if(!strcmp(attr_name,"BETA_G")) 
			{(*def).beta_g = atof(attr_value); continue;}
		if(!strcmp(attr_name,"RO_G")) 
			{(*def).ro0_g = atof(attr_value); continue;}
		if(!strcmp(attr_name,"MU_G")) 
			{(*def).mu_g = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_D_NW_0")) 
			{(*def).P_d_nw[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_D_NW_1")) 
			{(*def).P_d_nw[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_D_GN_0")) 
			{(*def).P_d_gn[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"P_D_GN_1")) 
			{(*def).P_d_gn[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_W_GR")) 
			{(*def).S_w_gr = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_G_GR")) 
			{(*def).S_g_gr = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_W_INIT")) 
			{(*def).S_w_init = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_G_INIT")) 
			{(*def).S_g_init = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_NR_0")) 
			{(*def).S_nr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_NR_1")) 
			{(*def).S_nr[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_GR_0")) 
			{(*def).S_gr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_GR_1")) 
			{(*def).S_gr[1] = atof(attr_value); continue;}
/*		if(!strcmp(attr_name,"S_W_RANGE_0")) 
			{(*def).S_w_range[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_W_RANGE_1")) 
			{(*def).S_w_range[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_G_RANGE_0")) 
			{(*def).S_g_range[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_G_RANGE_1")) 
			{(*def).S_g_range[1] = atof(attr_value); continue;}
*/
#else
		if(!strcmp(attr_name,"S_N_GR")) 
			{(*def).S_n_gr = atof(attr_value); continue;}
#endif

		if(!strcmp(attr_name,"SOURCE"))
			{(*def).source = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"ITERATIONS"))
			{(*def).newton_iterations = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"TIMEX"))
			{(*def).timeX = atof(attr_value); continue;}
		if(!strcmp(attr_name,"SAVE_PLOTS"))
			{(*def).save_plots = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"PRINT_SCREEN"))
			{(*def).print_screen = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"NX"))
			{(*def).Nx = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"NY"))
			{(*def).Ny = atoi(attr_value); continue;}
		if(!strcmp(attr_name,"NZ"))
			{(*def).Nz = atoi(attr_value); continue;}
	}

	fclose(defs);

	read_defines_test(*def);
}
