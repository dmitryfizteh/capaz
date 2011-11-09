#include "defines.h"

// �������� �������
clock_t start_time, finish_time;  
// �������� ������� ��� ������ ����� ������������
double *HostBuffer;
double *DevBuffer;

int main(int argc, char* argv[])
{
#ifdef TWO_PHASE
	std::cout << "Two phase filtration by CAPAZ\n\n";
#endif
#ifdef THREE_PHASE
	std::cout << "Three phase filtration by CAPAZ\n\n";
#endif
#ifdef B_L
	std::cout << "Backley-Leverett filtration by CAPAZ\n\n";
#endif

	consts def;
	read_defines(argc, argv, &def);

	// ���������� ���������� ����������� � ������ �������� ����������
	int size=0, rank=0;
	// �������� ������ ������ ��������� ������� ����������
	ptr_Arrays HostArraysPtr;
	// GPU-������ ������ ��������� ������� ����������
	ptr_Arrays DevArraysPtr;
	// ������� (���������) ��������� ������� ����������
	int localNx=0, localNy=0;
	// ���������� ������ ����������
	int blocksX=0, blocksY=0, blocksZ=0;
	// ������� ����� �� �������
	int j=0;
	
	// ������������� ������������, ������� ���������� ���������� � ��������� ����������, 
	// ��������� ������, �������� ���������/����������� ������
	initialization(&HostArraysPtr, &DevArraysPtr, &j, &localNx, &localNy, &size, &rank, &blocksX, &blocksY, &blocksZ, argc, argv, def);
	// ����
	save_data_plots(HostArraysPtr, DevArraysPtr, 0, size, rank, localNx, def);
	
	start_time=clock();

	// ���� ����� �� ������� (������ �������� - ����� ���� �� �������)
	// 1. ���������� ������� P1 � S2 �� ��������� ��������� ����
	// 2. ������ (def.print_screen) ��� �� ����� ��������� ���������� � ��������� ����
	// 3. ������ save_plots ��� ������ ����������� � ������ ����� � 
	//    ����������� � ����� �������� (**), � ���� ����������� ��������� ������ (***)
	for (j++; j <= def.timeX/(def.dt); j++)
	{
		time_step_function(HostArraysPtr, DevArraysPtr, DevBuffer, def,j*(def.dt),localNx,localNy,rank,size,blocksX,blocksY,blocksZ); // (1)

		if ((j % (def.print_screen)) == 0) // (2)
		{
			printf ("t=%.3f\n",j*(def.dt)); 
			fflush(stdout);
		}
		if ((j % (def.save_plots)) == 0) // (3)
		{
			// ��������� 2 ������� ���������� ������ � ����� �������,
			// �.�. save ���������� ������, ����������� save_data_plots
			save_data_plots(HostArraysPtr, DevArraysPtr, j*(def.dt), size, rank, localNx, def); // (**)
			//save(HostArraysPtr, DevArraysPtr, j, rank, size, localNx); // (***)
		}
	}

	// ����� ���������� � ������� ������ ��������� � ��������
	finish_time=clock();
	finish_time-=start_time;
	printf( "Task time in seconds:\t%.2f\n", (double) finish_time/CLOCKS_PER_SEC);

	// ���������� ������ � ������������ ������
	finalization(HostArraysPtr, DevArraysPtr, DevBuffer);

	// ��� ������� � Windows ����� ������ ��������� ��������� ���� �������
#ifdef _WIN32
	getchar();
#endif
	return 0;
}

// ������� ������� ����� �������� �� ��������� ��������� ����
// 1. ������ ��������� ��������� ro, �������� NAPL P2, ���������� Xi
// 2. ����� ����� ������������ ������������ ���������� P2, ro � Xi
// 3. ������ ��������� ���������
// 4. ����� ����� ������������ ������������ ���������� ��������� ���������
// 5. ������ ���������� roS �� ��������� ��������� ����
// 6. ������ ������� ������� �������� ���� P1 � ������������ DNAPL S2
// 7. ���������� ��������� ������� ��� P1 � S2
// 8. ����� ����� ������������ ������������ ���������� P1 � S2
void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t, int localNx, int localNy, int rank,int size, int blocksX, int blocksY, int blocksZ)
{
	P_S_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, localNx,blocksY, blocksZ, rank, size, def); // (8)
	ro_P_Xi_calculation(HostArraysPtr,DevArraysPtr,def,localNx,rank,size,blocksX,blocksY, blocksZ); // (1)
	P_ro_Xi_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, localNx,blocksY, blocksZ, rank, size, def); // (2)
	u_calculation(HostArraysPtr,DevArraysPtr,localNx,rank,size,blocksX,blocksY, blocksZ, def); // (3)
	u_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, localNx,blocksY, blocksZ, rank, size, def); // (4)
	roS_calculation(HostArraysPtr,DevArraysPtr,def,t,localNx,rank,size,blocksX,blocksY, blocksZ); // (5)
	P_S_calculation(HostArraysPtr,DevArraysPtr,def,localNx,rank,size,blocksX,blocksY, blocksZ); // (6)
	boundary_conditions(HostArraysPtr,DevArraysPtr,localNx,rank,size,blocksX,blocksY, blocksZ, def); // (7)
	
}

// �������������� ��������� ��������� ���������� � ����������
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������ 
// (���� 2 ������ � ����� ������,�� +2 �����). 
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� rank==0)
int i_to_I(int i, int rank, int size, consts def)
{
	int I;
	if (rank <= (def.Nx)%size)
	{
		if(rank==0)
			I=i;
		else
			I=((def.Nx)/size+1)*rank+i-1;
	}
	else
		I=((def.Nx)/size+1)*rank-(rank-(def.Nx)%size)+i-1;

	test_positive(I, __FILE__, __LINE__);
	return I;
}

// ���������� ��������� ������� (��������) ����������
// ���� ������� �� ��������������, �� ������ (def.NX)%size �������� +1 � ��������.
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������ 
// (���� 2 ������ � ����� ������,�� +2 �����). 
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� rank==0)
void global_to_local_vars (int* localNx, int* localNy, int size, int rank, consts def)
{
	*localNx=(def.Nx)/size;
	if (rank < (def.Nx)%size)
		(*localNx)++;

	// ������� ���������� �������� �� 1 ����� ��� ��������� ������,
	// ��������� - �� 2 �� ��� �������
	// ���� ��������� ����, �� ������ � ���� ��� � �������������� ����� �� �����
	if (size>1)
	{
		if ((rank==0) || (rank==size-1))
			(*localNx)++;
		else
			(*localNx)+=2;
	}

	*localNy=(def.Ny);
	/*
	*localNy=(def.Ny)/size;
	if (rank < (def.Ny)%size)
		*localNy++;
	*/
	test_positive(*localNx, __FILE__, __LINE__);
	test_positive(*localNy, __FILE__, __LINE__);
}

// �������� �� ����� �������� (�.�. �� ��������������� ������ ��� ������ �� ��������)
int is_active_point(int i, int localNx, int rank, int size)
{
	if((rank!=0 && i==0) || (rank!=size-1 && i==localNx-1))
		return 0;
	else
		return 1;
}

// ���������� ��������� ������ �� ���� ������
void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def)
{
	*t=0;
	for(int i=0;i<localNx;i++)
		for(int j=0;j<localNy;j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
					{
						// �������������� ��������� ��������� ���������� � ����������
						int I=i_to_I(i,rank,size,def);

						HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]=I*(def.hx);
						HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]=j*(def.hy);
						HostArraysPtr.z[i+j*localNx+k*localNx*(def.Ny)]=k*(def.hz);

#ifdef THREE_PHASE
						int j1 = (def.Ny);

						if(j < j1)
						{
							HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
							HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = def.S_g_gr + (def.S_g_init - def.S_g_gr) * j / j1;
						}
						else
						{
							HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = def.S_w_init;
							HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = def.S_g_init;
						}

						if(j == 0)
							HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
						else
							HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i + (j - 1) * localNx + k * localNx * (def.Ny)]
								+ (def.ro0_n * (1. - HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)]) 
								+ def.ro0_w * HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)]
								+ def.ro0_g * HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);

						HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)] = 0;						
#else
						// ���� ����� �� ������� �������, �� ����� (def.source) ����� �� ������,
						// �� � ��� ��������� ������������. �����, �������
						if ((j==0) && (I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
							HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]=def.S_n_gr;
						else
							HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]=0;

						HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]=def.P_atm+j * (def.ro0_w) * (def.g_const)*(def.hy);
						HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)]=0;
#endif
					
						/*
						if ((HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]>=(def.NX)/2.*(def.h1)) && (HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]<=4.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]<=2./5.*(def.Ny)*(def.h2)) && (HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]>=(-1.)*HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]/4.+2./5.*(def.Ny)*(def.h2)))
								HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)]=1;

						if ((HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]>=(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]<=2.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]<=4./5.*(def.Ny)*(def.h2)) && (HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]>=3./5.*(def.Ny)*(def.h2)))
								HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)]=1;
								*/
					
						/*
						if ((HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]>=2.*(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]<=3.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]>=1./10.*(def.Ny)*(def.h2)) && (HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]<=3./10.*(def.Ny)*(def.h2)))
								HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)]=1;
						*/
#ifndef THREE_PHASE 
					test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
#else
					test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
#endif
					}
}

//----------------------------------------------------------------------------------------------------
// ��������� �������

// ������������� ������������ (1), ������� ���������� ���������� � ��������� ���������� (2), 
// ������������� ���������� (2.5), ��������� ������ (3), �������� ���������/����������� ������ (4)
void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int* j, int* localNx, int* localNy, int* size, int* rank, int* blocksX, int* blocksY, int* blocksZ, int argc, char* argv[], consts def)
{
	FILE *f_save;

	communication_initialization(argc, argv, size, rank, def); // (1)

	global_to_local_vars(localNx, localNy, *size, *rank, def); // (2)

	device_initialization(*rank, blocksX, blocksY, blocksZ, *localNx, def); // (2.5)

	memory_allocation(HostArraysPtr, DevArraysPtr, (def.Nx)/(*size)+3, (def.Ny), def); // (3)

	// ���� ������� ��������� ����� ������� ���� ������������ ���������,
	// �� ��������������� ���������, ����� ��������� ��������� �������
	if (f_save=fopen("save/save.dat","rb"))
	{
		fclose(f_save);
		restore(*HostArraysPtr, j, *rank, *size, *localNx, def);
	}
	else
		data_initialization (*HostArraysPtr, j, *localNx, *localNy, *rank, *size, def); // (4)

#ifdef THREE_PHASE
	load_data_to_device((*HostArraysPtr).P_n, (*DevArraysPtr).P_n, *localNx, def);
	load_data_to_device((*HostArraysPtr).S_w, (*DevArraysPtr).S_w, *localNx, def);
	load_data_to_device((*HostArraysPtr).S_g, (*DevArraysPtr).S_g, *localNx, def);
	load_data_to_device((*HostArraysPtr).roS_g_old, (*DevArraysPtr).roS_g_old, *localNx, def);
#else
	load_data_to_device((*HostArraysPtr).P_w, (*DevArraysPtr).P_w, *localNx, def);
	load_data_to_device((*HostArraysPtr).S_n, (*DevArraysPtr).S_n, *localNx, def);
#endif
	load_data_to_device((*HostArraysPtr).roS_w_old, (*DevArraysPtr).roS_w_old, *localNx, def);
	load_data_to_device((*HostArraysPtr).roS_n_old, (*DevArraysPtr).roS_n_old, *localNx, def);
	load_data_to_device_int((*HostArraysPtr).media, (*DevArraysPtr).media, *localNx, def);
}

// ���������� ������ (1), ������������ ������ (2)
void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	memory_free(HostArraysPtr, DevArraysPtr); // (2)
	communication_finalization(); // (1)
}

// ��������� ������ ����� (1) � ���������� (2) ��� ������ ����� ��������� �������
void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int localNx, int nY, consts def)
{
	host_memory_allocation(HostArraysPtr, localNx, nY, def); // (1)
	device_memory_allocation(DevArraysPtr, &DevBuffer, localNx, def); // (2)
}

// ����������� ������ ����� (1) � ���������� (2) �� ��� ������� ����� ��������� �������
void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr)
{
	host_memory_free(HostArraysPtr); // (1)
	device_memory_free(DevArraysPtr, DevBuffer); // (2)
}

// ��������� ������ ����� ��� ������ ����� ��������� �������
void host_memory_allocation(ptr_Arrays* ArraysPtr, int localNx, int nY, consts def)		
{	
	if (!(HostBuffer=new double[2*(def.Ny)*(def.Nz)]))
		printf ("\nWarning! Memory for *HostBuffer is not allocated in function host_memory_alloc\n");

	try
	{
		(*ArraysPtr).x=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).y=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).z=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).P_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).P_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ro_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ro_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ux_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uy_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uz_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ux_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uy_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uz_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).Xi_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).Xi_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_w_old=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_n=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_n_old=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).media=new int [localNx*nY*(def.Nz)];
#ifdef THREE_PHASE
		(*ArraysPtr).P_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).S_w=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).S_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ro_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).ux_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uy_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).uz_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).Xi_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_g=new double [localNx*nY*(def.Nz)];
		(*ArraysPtr).roS_g_old=new double [localNx*nY*(def.Nz)];
#else
		(*ArraysPtr).S_n=new double [localNx*nY*(def.Nz)];
#endif
	}
	catch(...)
	{
		printf ("\nError! Not enough host memory\n");
		exit(0);
	}
}

// ����������� ������ ����� �� ��� ������� ����� ��������� �������
void host_memory_free(ptr_Arrays ArraysPtr)
{ 
	delete HostBuffer;
	delete[] ArraysPtr.x;
	delete[] ArraysPtr.y;
	delete[] ArraysPtr.z;
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

// ������� ���������� �������� � �����
void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, int size, int rank, int localNx, consts def)
{
	// �������� � ������ ����� ����������� �������
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.P_n, DevArraysPtr.P_n , localNx, def);
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtr.S_w , localNx, def);
	load_data_to_host(HostArraysPtr.S_g, DevArraysPtr.S_g , localNx, def);
#else
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtr.P_w , localNx, def);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtr.S_n , localNx, def);
#endif
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtr.ux_n , localNx, def);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtr.uy_n , localNx, def);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtr.uz_n , localNx, def);

#ifndef THREE_PHASE
	// �������� �� ����� �� ����������� ��������� �������� P � S
#ifdef MY_TEST
	test_correct_P_S(HostArraysPtr, localNx, rank, def);
#endif
#endif
	
	// ������� ��������� ������� ����������, ����� � ����������� ��������� ������
	if (rank==0)
		print_plots_top (t, def);

	// �� ������� ��� ������� �� ����������� �������� ������� ������ �� ������
	// ����� ����� �������.
	for (int cpu=0; cpu<size;cpu++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
		barrier();
		if (rank==cpu)
			print_plots(HostArraysPtr, t, rank, size, localNx,def);	
	}
}

// ������� �������� ����������, ������ ��� �������� � ���������� ���������� � ��� 
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

	// �������� (��� ����������) ����� � ���������
	// 1. ��� ������������� ������������� NAPL S_n
	// 2. ��� ������������� �������� ���� P_w
	// 3. ��� ������������� ��������� {u_x, u_y, u_z}
	// 4. ��� ������������� ����� �������
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


// ������� ���������� ������ � ����� �������� (!3D)
void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, int size, int localNx, consts def)
{
	char fname[30];
	FILE *fp;
	int local;

	sprintf(fname,"plots/S=%012.4f.dat",t);
		
	// �������� �� �������� � ���������� ��������
	// 1. ��� ������������� ������������� NAPL S_n
	// 2. ��� ������������� �������� ���� P_w
	// 3. ��� ������������� ��������� {u_x, u_y, u_z}
	// 4. ��� ������������� ����� �������
	if(!(fp=fopen(fname,"at")))
		std::cout << "Not open file(s) in function SAVE_DATA_PLOTS! \n";

	for(int i=0; i<localNx; i++)
		for(int j=0; j<(def.Ny); j++)
			for(int k=0; k<(def.Nz); k++)
				if(is_active_point(i, localNx, rank, size))
				{
					local=i+j*localNx+k*localNx*(def.Ny);
#ifdef THREE_PHASE
					if((def.Nz) < 2)
					{
/*						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/					
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.ux_n[local], 
							(-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
		
					}

					else
					{
/*						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_g[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], 
							HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
							HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local],
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);

					}
#else
					if((def.Nz) < 2)
					{
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local], HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)

					}
					else
					{
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %d\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny-1)*(def.hy)-HostArraysPtr.y[local], HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)
					}
#endif
				}

	/* �� ����� ������, ��� ��� ������ �������, �������� �������� � �����
	for(int i=0; i<localNx; i++)
		for(int k=0; k<(def.Nz); k++)
			if ((i_to_I(i,rank,size, def)==(def.Nx)/2) && is_active_point(i,localNx,rank,size))
				for(int j=0; j<(def.Ny); j++)
					fprintf(fp_S2y,"%.2e %.3e\n", HostArraysPtr.y[localNx/2+j*localNx+k*localNx*(def.Ny)], HostArraysPtr.S_n[localNx/2+j*localNx+k*localNx*(def.Ny)]); 
	

	for(int i=0; i<localNx; i++)
		for(int k=0; k<(def.Nz); k++)
			fprintf(fp_S2x,"%.2e %.3e\n", HostArraysPtr.x[i+localNx*(def.Ny)/2+k*localNx*(def.Ny)], HostArraysPtr.S_n[i+localNx*(def.Ny)/2+k*localNx*(def.Ny)]); 
	*/
		
	fclose(fp);
}


// ���������� ��������� � ����
void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, int rank, int size, int localNx, consts def)
{
	// ��������� � ������ ����� ������ �� roS_old
	// P1 � S2 ��������� ��� ��� ������� ���������� ��������,
	// x,y � media �� ���������� � �������� �� �������.
	//load_data_to_host(HostArraysPtr.P1, DevArraysPtr.P1 , localNx);
	//load_data_to_host(HostArraysPtr.S2, DevArraysPtr.S2 , localNx);
	load_data_to_host(HostArraysPtr.roS_w_old, DevArraysPtr.roS_w_old , localNx, def);
	load_data_to_host(HostArraysPtr.roS_n_old, DevArraysPtr.roS_n_old , localNx, def);
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.roS_g_old, DevArraysPtr.roS_n_old , localNx, def);
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

	for (int cpu=0; cpu<size;cpu++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
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
			fwrite(HostArraysPtr.P_n, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.S_w, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.S_g, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
#else
			fwrite(HostArraysPtr.P_w, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.S_n, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
#endif
			fwrite(HostArraysPtr.x, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.y, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.z, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.roS_w_old, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
			fwrite(HostArraysPtr.roS_n_old, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.roS_g_old, sizeof(double), localNx * (def.Ny) * (def.Nz), f_save);
#endif
			fwrite(HostArraysPtr.media, sizeof(int), localNx * (def.Ny) * (def.Nz), f_save);
			fclose(f_save);
		}
	}
}

// �������������� ��������� �� �����
void restore (ptr_Arrays HostArraysPtr, int* j, int rank, int size, int localNx, consts def)
{
	FILE *f_save;
	for (int cpu=0; cpu<size;cpu++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
		barrier();
		int lNx=0, lNy=0;
		if (rank==cpu)
		{
			if(!(f_save=fopen("save/save.dat","rb")))
			{
					printf("\nError: Not open file \"save.dat\"!\n");
					exit(0);
			}
			for (int queue=0;queue<=rank;queue++)
			{
				global_to_local_vars(&lNx,&lNy,size,queue, def);
				fread(j, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.P_n, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.S_w, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.S_g, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
#else
				fread(HostArraysPtr.P_w, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.S_n, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
#endif
				fread(HostArraysPtr.x, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.y, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.z, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.roS_w_old, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
				fread(HostArraysPtr.roS_n_old, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.roS_g_old, sizeof(double), lNx * (def.Ny) * (def.Nz), f_save);
#endif
				fread(HostArraysPtr.media, sizeof(int), lNx * (def.Ny) * (def.Nz), f_save);
			}
			fclose(f_save);
		}
	}
}



//------------------------------------------------------------------------------------------

// ���������� ���������� ������ �� �����
void read_defines(int argc, char *argv[], consts* def)
{
#ifdef TWO_PHASE
	(*def).P_d[0]=P_d[0];
	(*def).P_d[1]=P_d[1];
#endif
#ifndef THREE_PHASE
	(*def).K[0]=K[0];
	(*def).K[1]=K[1];
	(*def).lambda[0]=lambda[0];
	(*def).lambda[1]=lambda[1];
	(*def).S_wr[0]=S_wr[0];
	(*def).S_wr[1]=S_wr[1];
	(*def).m[0]=m[0];
	(*def).m[1]=m[1];
#else
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
#ifdef THREE_PHASE
		if(!strcmp(attr_name,"LAMBDA_0")) 
		{(*def).lambda[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"LAMBDA_1")) 
		{(*def).lambda[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"K_0")) 
		{(*def).K[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"K_1")) 
		{(*def).K[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"M_0")) 
		{(*def).m[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"M_1")) 
		{(*def).m[1] = atof(attr_value); continue;}

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
		if(!strcmp(attr_name,"S_WR_0")) 
		{(*def).S_wr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_WR_1")) 
		{(*def).S_wr[1] = atof(attr_value); continue;}
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
		{(*def).timeX = atoi(attr_value); continue;}
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
