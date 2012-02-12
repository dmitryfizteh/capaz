#include "defines.h"

// �������� ������� ��� ������ ����� ������������
double *HostBuffer;
double *DevBuffer;

int main(int argc, char* argv[])
{
	consts def;
	read_defines(argc, argv, &def);

	// �������� ������ ������ ��������� ������� ����������
	ptr_Arrays HostArraysPtr;
	// GPU-������ ������ ��������� ������� ����������
	ptr_Arrays DevArraysPtr;
	// ������� ����� �� �������
	int j=0;
	// ������� ������� ���������� �������������� ����� ���������
	clock_t task_time; 

	// ������������� ������������, ������� ���������� ���������� � ��������� ����������, 
	// ��������� ������, �������� ���������/����������� ������
	initialization(&HostArraysPtr, &DevArraysPtr, &j, argc, argv, &def);

	// ����
	//save_data_plots(HostArraysPtr, DevArraysPtr, 0, def);
	
	task_time=clock();

	// ���� ����� �� ������� (������ �������� - ����� ���� �� �������)
	// 1. ���������� ������� P1 � S2 �� ��������� ��������� ����
	// 2. ������ (def.print_screen) ��� �� ����� ��������� ���������� � ��������� ����
	// 3. ������ save_plots ��� ������ ����������� � ������ ����� � 
	//    ����������� � ����� �������� (**), � ���� ����������� ��������� ������ (***)
	for (j++; j <= def.timeX/(def.dt); j++)
	{
		if ((j % (def.print_screen) == 0) && (def.rank) == 0) // (2)
		{
			printf ("t=%.3f\n", j * (def.dt)); 
			fflush(stdout);
		}

		time_step_function(HostArraysPtr, DevArraysPtr, DevBuffer, def, j * (def.dt)); // (1)

		if ((j % (def.save_plots)) == 0) // (3)
		{
			// ��������� 2 ������� ���������� ������ � ����� �������,
			// �.�. save ���������� ������, ����������� save_data_plots
			save_data_plots(HostArraysPtr, DevArraysPtr, j * (def.dt), def); // (**)
			//save(HostArraysPtr, DevArraysPtr, j, def); // (***)
		}
	}

	// ���� ���� ��� ���������� �������� �������
	barrier();
	// ����� ���������� � ������� ������ ��������� � ��������
	task_time=clock()-task_time;
	if(!(def.rank))
		printf( "Task time in seconds:\t%.2f\n", (double) task_time/CLOCKS_PER_SEC);

	// ���������� ������ � ������������ ������
	finalization(HostArraysPtr, DevArraysPtr, DevBuffer);

	// ��� ������� � Windows ����� ������ ��������� ��������� ���� �������
#ifdef _WIN32
	printf("\nPress <Enter> to exit...\n");
	fflush( stdout);
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
void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t)
{
	P_S_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (8)
	ro_P_Xi_calculation(HostArraysPtr, DevArraysPtr, def); // (1)
	P_ro_Xi_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (2)
	u_calculation(HostArraysPtr, DevArraysPtr, def); // (3)
	u_exchange(HostArraysPtr, DevArraysPtr, HostBuffer, DevBuffer, def); // (4)
	roS_calculation(HostArraysPtr, DevArraysPtr, t, def); // (5)
	P_S_calculation(HostArraysPtr, DevArraysPtr, def); // (6)
	boundary_conditions(HostArraysPtr, DevArraysPtr, def); // (7)
	
}

// �������������� ��������� ��������� ���������� � ����������
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������ 
// (���� 2 ������ � ����� ������,�� +2 �����). 
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� (def.rank)==0)
int i_to_I(int i, consts def)
{
	int I;
	if ((def.rank) <= (def.Nx) % (def.sizex))
	{
		if((def.rank) == 0)
			I = i;
		else
			I = ((def.Nx) / (def.sizex) + 1) * (def.rank) + i - 1;
	}
	else
		I = ((def.Nx) / (def.sizex) + 1) * (def.rank) - ((def.rank) - (def.Nx) % (def.sizex)) + i - 1;

	test_positive(I, __FILE__, __LINE__);
	return I;
}

// ���������� ��������� ������� (��������) ����������
// ���� ������� �� ��������������, �� ������ (def.NX)%size �������� +1 � ��������.
// ������ ��������� �������� �������������� ����� � ������� ���
// ������ �������, ���� ����� ������ 
// (���� 2 ������ � ����� ������,�� +2 �����). 
// ���������� ������� �������� ��� ������� ����� (������ � ������� �� (def.rank)==0)
void global_to_local_vars (consts *def)
{
	(*def).locNx = (*def).Nx / (*def).sizex;

	if ((*def).rank % (*def).sizex < (*def).Nx % (*def).sizex)
		((*def).locNx) ++;

	// ������� ���������� �������� �� 1 ����� ��� ��������� ������,
	// ��������� - �� 2 �� ��� �������
	// ���� ��������� ����, �� ������ � ���� ��� � �������������� ����� �� �����
	if ((*def).sizex > 1)
	{
		if (((*def).rank % (*def).sizex == 0) || ((*def).rank % (*def).sizex  == (*def).sizex - 1))
			((*def).locNx) ++;
		else
			((*def).locNx) += 2;
	}

	(*def).locNy = (*def).Ny / (*def).sizey;

	if (((*def).rank % ((*def).sizex) * (*def).sizey) / (*def).sizex < (*def).Ny % (*def).sizey)
		((*def).locNy) ++;

	if ((*def).sizey > 1)
	{
		if ((((*def).rank % ((*def).sizex * (*def).sizey)) / (*def).sizex == 0) || (((*def).rank % ((*def).sizex * (*def).sizey)) / (*def).sizex == (*def).sizey - 1))
			((*def).locNy) ++;
		else
			((*def).locNy) += 2;
	}

	(*def).locNz = (*def).Nz / (*def).sizez;

	if ((*def).rank / (*def).sizex / (*def).sizey < (*def).Nz % (*def).sizez)
		((*def).locNz) ++;

	if ((*def).sizez > 1)
	{
		if (((*def).rank / (*def).sizex / (*def).sizey == 0) || ((*def).rank / (*def).sizex / (*def).sizey == (*def).sizez - 1))
			((*def).locNz) ++;
		else
			((*def).locNz) += 2;
	}

	test_positive((*def).locNx, __FILE__, __LINE__);
	test_positive((*def).locNy, __FILE__, __LINE__);
	test_positive((*def).locNz, __FILE__, __LINE__);
}

// �������� �� ����� �������� (�.�. �� ��������������� ������ ��� ������ �� ��������)
int is_active_point(int i, int j, int k, consts def)
{
	if(((def.rank) % (def.sizex) != 0 && i == 0) || ((def.rank) % (def.sizex) != (def.sizex) - 1 && i == def.locNx - 1)
		|| (((def.rank) % ((def.sizex) * (def.sizey))) / (def.sizex) != 0 && j == 0)	|| (((def.rank) % ((def.sizex) * (def.sizey))) / (def.sizex) != (def.sizey) - 1 && j == def.locNy - 1)
		|| ((((def.rank) / (def.sizex) / (def.sizey) != 0 && k == 0) || ((def.rank) / (def.sizex) / (def.sizey) == (def.sizez) - 1 && k == def.locNz - 1)) && (def.sizez) > 1))
		return 0;
	else
		return 1;
}

// ���������� ��������� ������ �� ���� ������
void sizes_initialization(consts *def)
{
	(*def).sizex = (*def).size;
	(*def).sizey = 1;
	(*def).sizez = 1;
}

void blocks_initialization(consts *def)
{
	(*def).blocksX = 0;
	(*def).blocksY = 0;
	(*def).blocksZ = 0;	
}

// ������� ���������� "�����������" ��������� * g * hy
double ro_eff_gdy(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

#ifdef THREE_PHASE
	double ro_g_dy = (HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) *(def.locNy)] * (1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) 
		+ HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) *(def.locNy)] * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
	+ HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) *(def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]) * (def.m[media]) * (def.g_const) * (def.hy);
#else
	double ro_g_dy = (HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) *(def.locNy)] * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] 
	+ HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) *(def.locNy)] * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[media]) * (def.g_const) * (def.hy);
#endif
	return ro_g_dy;
}

void data_initialization(ptr_Arrays HostArraysPtr, int* t, consts def)
{
	*t = 0;
	for(int i = 0; i < def.locNx; i++)
		for(int j = 0; j < def.locNy; j++)
			for(int k = 0; k < def.locNz; k++)
				if(is_active_point(i, j, k, def))
					{
						// �������������� ��������� ��������� ���������� � ����������
						int I = i_to_I(i, def);

#ifdef THREE_PHASE
						int media = HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy] = 0;	
						int j1 = def.locNy / 2;

						if(j < j1)
						{
							HostArraysPtr.S_w[i+j*def.locNx+k*def.locNx*def.locNy] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
							HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy] = def.S_n_gr + (def.S_n_init - def.S_n_gr) * j / j1;
						}
						else
						{
							HostArraysPtr.S_w[i+j*def.locNx+k*def.locNx*def.locNy] = def.S_w_init;
							HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy] = def.S_n_init;
						}

						if(j == 0)
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy] = def.P_atm;
						else
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_eff_gdy(HostArraysPtr, i, j-1, k, def);
	
						HostArraysPtr.ro_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

						///!!!! �� ����������� ����������� ����! ��� ���� ������� ����� ���� ����� P_w, P_g
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
						HostArraysPtr.ro_g[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));	
#endif
#ifdef TWO_PHASE
						// ���� ����� �� ������� �������, �� ����� (def.source) ����� �� ������,
						// �� � ��� ��������� ������������. �����, �������
						if ((j==0) && (I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
							HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy]=def.S_n_gr;
						else
							HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy]=0;

						if(j == 0)
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=def.P_atm;
						else
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=HostArraysPtr.P_w[i+(j-1)*def.locNx+k*def.locNx*def.locNy] + ro_eff_gdy(HostArraysPtr, i, j-1, k, def);
							
						HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=0;

						HostArraysPtr.ro_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

						///!!!! �� ����������� ����������� ����! ��� ���� ������� ����� ���� ����� P_n
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
#endif

#ifdef B_L
						HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=0;
						HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy]=0.5;

						// �� ����������� ���� �������
						HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=def.P_atm;

						// � ����� ������ ���� ����������� ��������
						if (((i==0) && (j==def.Ny-2)) || ((i==1) && (j==def.Ny-1)) || ((i==0) && (j==def.Ny-1)) || ((i==1) && (j==def.Ny-2)))
						{
							HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 1e6;
							//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0;
						}

						HostArraysPtr.ro_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

						///!!!! �� ����������� ����������� ����! ��� ���� ������� ����� ���� ����� P_n
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

#endif

						/*
						if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=(def.NX)/2.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=4.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=2./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=(-1.)*HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]/4.+2./5.*def.locNy*(def.h2)))
								HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;

						if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=2.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=4./5.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=3./5.*def.locNy*(def.h2)))
								HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;
								*/
					
						/*
						if ((HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]>=2.*(def.NX)/5.*(def.h1)) && (HostArraysPtr.x[i+j*def.locNx+k*def.locNx*def.locNy]<=3.*(def.NX)/5.*(def.h1)))
							if ((HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]>=1./10.*def.locNy*(def.h2)) && (HostArraysPtr.y[i+j*def.locNx+k*def.locNx*def.locNy]<=3./10.*def.locNy*(def.h2)))
								HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=1;
						*/

					test_nan(HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
#ifdef THREE_PHASE 
					test_nan(HostArraysPtr.S_w[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
#endif
					}
}

//----------------------------------------------------------------------------------------------------
// ��������� �������

// ����� ����������� ������
void print_task_name(consts def)
{
	// ������� ��������� ������� �������� ����������� ������
	if (!(def.rank))
	{
		#ifdef TWO_PHASE
				std::cout << "Two phase filtration by CAPAZ on "<<(def.size)<<" node(s).\n";
		#endif
		#ifdef THREE_PHASE
				std::cout << "Three phase filtration by CAPAZ on "<<(def.size)<<" node(s).\n";
		#endif
		#ifdef B_L
				std::cout << "Backley-Leverett filtration by CAPAZ on "<<(def.size)<<" node(s).\n";
		#endif
		read_version();
	    fflush(stdout);
	}
}

// ������������� ������������ (1), ������� ���������� ���������� � ��������� ���������� (2), 
// ������������� ���������� (2.5), ��������� ������ (3), �������� ���������/����������� ������ (4)
// ��� ������ �-� �������� �������������� �� �����.
void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int* j, int argc, char* argv[], consts* def)
{
	FILE *f_save;

	communication_initialization(argc, argv, def); // (1)

	print_task_name(*def);

	sizes_initialization(def);

	blocks_initialization(def);

	global_to_local_vars(def); // (2)

	device_initialization(def); // (2.5)

	memory_allocation(HostArraysPtr, DevArraysPtr, *def); // (3)

#ifdef B_L
	load_permeability((*HostArraysPtr).K, *def); // (5)
	load_data_to_device((*HostArraysPtr).K, (*DevArraysPtr).K, *def);
#endif

	// ���� ��������� ����� ������� ���� ������������ ���������,
	// �� ��������������� ���������, ����� ��������� ��������� �������
	if (f_save=fopen("save/save.dat","rb"))
	{
		fclose(f_save);
		restore(*HostArraysPtr, j, *def);
	}
	else
		data_initialization (*HostArraysPtr, j, *def); // (4)

#ifdef THREE_PHASE
	load_data_to_device((*HostArraysPtr).S_w, (*DevArraysPtr).S_w, *def);
	load_data_to_device((*HostArraysPtr).roS_g_old, (*DevArraysPtr).roS_g_old, *def);
#endif
	load_data_to_device((*HostArraysPtr).P_w, (*DevArraysPtr).P_w, *def);
	load_data_to_device((*HostArraysPtr).S_n, (*DevArraysPtr).S_n, *def);
	load_data_to_device((*HostArraysPtr).roS_w_old, (*DevArraysPtr).roS_w_old, *def);
	load_data_to_device((*HostArraysPtr).roS_n_old, (*DevArraysPtr).roS_n_old, *def);
	load_data_to_device_int((*HostArraysPtr).media, (*DevArraysPtr).media, *def);
}

// ���������� ������ (1), ������������ ������ (2)
void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer)
{
	memory_free(HostArraysPtr, DevArraysPtr); // (2)
	communication_finalization(); // (1)
	device__finalization(); // (1)
}

// ��������� ������ ����� (1) � ���������� (2) ��� ������ ����� ��������� �������
void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, consts def)
{
	host_memory_allocation(HostArraysPtr, def); // (1)
	device_memory_allocation(DevArraysPtr, &DevBuffer, def); // (2)
}

// ����������� ������ ����� (1) � ���������� (2) �� ��� ������� ����� ��������� �������
void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr)
{
	host_memory_free(HostArraysPtr); // (1)
	device_memory_free(DevArraysPtr, DevBuffer); // (2)
}

// ��������� ������ ����� ��� ������ ����� ��������� �������
void host_memory_allocation(ptr_Arrays* ArraysPtr, consts def)		
{	
	if (!(HostBuffer=new double[2 * ((def.locNy) * (def.locNz))]))
		printf ("\nWarning! Memory for *HostBuffer is not allocated in function host_memory_alloc\n");

	try
	{
		(*ArraysPtr).S_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).P_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).P_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ro_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ro_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ux_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uy_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uz_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ux_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uy_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uz_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).Xi_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).Xi_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_w_old=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_n=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_n_old=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).media=new int [(def.locNx)*(def.locNy)*(def.locNz)];
#ifdef B_L
		(*ArraysPtr).K=new double [(def.locNx)*(def.locNy)*(def.locNz)];
#endif
#ifdef THREE_PHASE
		(*ArraysPtr).P_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).S_w=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ro_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).ux_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uy_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).uz_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).Xi_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_g=new double [(def.locNx)*(def.locNy)*(def.locNz)];
		(*ArraysPtr).roS_g_old=new double [(def.locNx)*(def.locNy)*(def.locNz)];
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
	delete[] ArraysPtr.ro_g;
	delete[] ArraysPtr.ux_g;
	delete[] ArraysPtr.uy_g;
	delete[] ArraysPtr.uz_g;
	delete[] ArraysPtr.Xi_g;
	delete[] ArraysPtr.roS_g;
	delete[] ArraysPtr.roS_g_old;
#endif
	delete[] ArraysPtr.S_n;
}

// ������� ���������� �������� � �����
void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def)
{
	// �������� � ������ ����� ����������� �������
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.S_w, DevArraysPtr.S_w , def);
	load_data_to_host(HostArraysPtr.ux_w, DevArraysPtr.ux_w , def);
	load_data_to_host(HostArraysPtr.uy_w, DevArraysPtr.uy_w , def);
	load_data_to_host(HostArraysPtr.uz_w, DevArraysPtr.uz_w , def);
	load_data_to_host(HostArraysPtr.ux_g, DevArraysPtr.ux_g , def);
	load_data_to_host(HostArraysPtr.uy_g, DevArraysPtr.uy_g , def);
	load_data_to_host(HostArraysPtr.uz_g, DevArraysPtr.uz_g , def);
#endif
	load_data_to_host(HostArraysPtr.P_w, DevArraysPtr.P_w , def);
	load_data_to_host(HostArraysPtr.S_n, DevArraysPtr.S_n , def);
	load_data_to_host(HostArraysPtr.ux_n, DevArraysPtr.ux_n , def);
	load_data_to_host(HostArraysPtr.uy_n, DevArraysPtr.uy_n , def);
	load_data_to_host(HostArraysPtr.uz_n, DevArraysPtr.uz_n , def);

#ifndef THREE_PHASE
	// �������� �� ����� �� ����������� ��������� �������� P � S
#ifdef MY_TEST
	test_correct_P_S(HostArraysPtr, def);
#endif
#endif
	
	// ������� ��������� ������� ����������, ����� � ����������� ��������� ������
	if ((def.rank)==0)
		print_plots_top (t, def);

	// �� ������� ��� ������� �� ����������� �������� ������� ������ �� ������
	// ����� ����� �������.
	for (int cpu = 0; cpu < (def.sizex * (def.sizey) * (def.sizez)); cpu ++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
		barrier();
		if ((def.rank)==cpu)
			print_plots(HostArraysPtr, t, def);	
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
//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"u_x\",\"u_y\",\"media\" \n");
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"un_x\",\"un_y\",\"ug_x\",\"ug_y\",\"media\" \n");
#else
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"S_n\",\"P_w\",\"u_x\", \"u_y\",\"media\" \n");
#endif
		fprintf(fp,"ZONE T = \"BIG ZONE\", K=%d,J=%d, F = POINT\n", (def.Nx), (def.Ny));
	}
	else
	{
#ifdef THREE_PHASE
//		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"u_x\", \"u_y\",\"u_z\",\"media\" \n");
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_w\",\"S_n\",\"S_g\",\"P_w\",\"uw_x\",\"uw_y\",\"uw_z\",\"un_x\",\"un_y\",\"un_z\",\"ug_x\",\"ug_y\",\"ug_z\",\"media\" \n");
#else
		fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\",\"S_n\",\"P_w\",\"u_x\", \"u_y\", \"u_z\", \"media\" \n");
#endif
		fprintf(fp,"ZONE T = \"BIG ZONE\", K=%d,J=%d,I=%d, F = POINT\n", (def.Nx), (def.Ny), (def.Nz));
	}

	fclose(fp);
}


// ������� ���������� ������ � ����� ��������
void print_plots(ptr_Arrays HostArraysPtr, double t, consts def)
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

	for(int i=0; i<def.locNx; i++)
		for(int j=0; j<def.locNy; j++)
			for(int k=0; k<def.locNz; k++)
				if(is_active_point(i, j, k, def))
				{
					local=i+j*def.locNx+k*def.locNx*def.locNy;

					// �������������� ��������� ��������� ���������� � ����������
					int I=i_to_I(i, def);
#ifdef THREE_PHASE
					if(def.Nz < 2)
					{
/*						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], 
							HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/					
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], 
							HostArraysPtr.ux_w[local], (-1)*HostArraysPtr.uy_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.ux_g[local], 
							(-1)*HostArraysPtr.uy_g[local], HostArraysPtr.media[local]);
		
					}

					else
					{
/*						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1)*(def.hy)-HostArraysPtr.y[local],  
							HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], 
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]);
*/
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1-j)*(def.hy),  
							HostArraysPtr.S_w[local], HostArraysPtr.S_n[local], 1. - HostArraysPtr.S_w[local] - HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], 
							HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local],
							HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local],
							HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local], HostArraysPtr.media[local]);

					}
#endif
#ifdef TWO_PHASE
					if(def.Nz < 2)
					{
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)

					}
					else
					{
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %d\n", I*(def.hx), k*(def.hz), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.media[local]); // (1)
					}
#endif
#if B_L
					if(def.Nz < 2)
					{
						fprintf(fp,"%.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", I*(def.hx), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.K[local]); // (1)

					}
					else
					{
						fprintf(fp,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e %.3e\n", I*(def.hx), k*(def.hz), (def.Ny-1-j)*(def.hy), HostArraysPtr.S_n[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local], HostArraysPtr.K[local]); // (1)
					}
#endif
				}

	
	fclose(fp);
}

// ������� ���������� ������ � ����� �������� ������� BjnIO [http://lira.imamod.ru/BjnIO_3D.html]
void print_plots_BjnIO(ptr_Arrays HostArraysPtr, double t, consts def)
{
/*
	char *dname;
	char *targetfuncs="r2d.bjn";
	char *targetgrid="r2dgrid.bjn";
	LPLPLPFLOAT F1;
	double *F1D;

	F1 = alloc_float_mas_n1_n2_n3(def.locNx, def.locNy, def.locNz); 
	if(F1 == NULL) exit(0);

	F1D = (double *)malloc(def.locNx, def.locNy, def.locNz*sizeof(double)); 
	if(F1D == NULL) exit(0);

	// ������������� ����� � ������� �������
	int err = WriteBjnGzippedScalar8RecInit(targetfuncs, "r2dgrid.bjn", def.Nx, def.Ny, def.Nz);   
	if(err) 
		fprintf(stderr, "Can't create file %s.\n", targetfuncs); 

	// ������ ����� ������ �������� �������
	err = WriteBjnGzippedScalar8RecFuncByBlock(targetfuncs, "S_n", F1D, def.locNx*(def.rank), 0, 0, def.locNx, def.locNy, def.locNz, 5); // ��������� def.locNx*(def.rank) �� ������ ��������
	if(err) 
		fprintf(stderr, "Can't add func block data err %d\n", err);

	// ������ ������������ � ������������� �������� �������� ������� S_n
	err = WriteBjnGzippedScalar8RecFuncMinMax(targetfuncs, "S_n", 0, 1);
	if(err) 
		fprintf(stderr, "Can't add data about block err %d\n", err);

	// ������������� ����� �����
	err = WriteBjnGzippedScalar8RecInit(targetgrid, "", def.Nx, def.Ny, def.Nz);      
	if(err) 
		fprintf(stderr, "Can't create file %s.\n", targetgrid);

	// ��� ������� �� �����������
	for(direct=0;direct<3;direct++)
	{
		for(i1=0; i1<m1; i1++)
			for(i2=0; i2<m2; i2++)
				for(i3=0; i3<m3; i3++)
				{
					for(i=0; i<n1; i++)     
						for(j=0; j<n2; j++)     
							for(k=0; k<n3; k++)
							{
								float x=(float)i1*(float)n1+(float)i;
								float y=(float)i2*(float)n2+(float)j;
								float z=(float)i3*(float)n3+(float)k;
								switch(direct)
								{
									case 0: F1[i][j][k] = (float)x; dname="x"; break;
									case 1: F1[i][j][k] = (float)y; dname="y"; break;
									case 2: F1[i][j][k] = (float)(ffmin+k*(ffmax-ffmin)/3.f); dname="z"; break;
								}

								fmin=minab(fmin,F1[i][j][k]); 
								fmax=maxab(fmax,F1[i][j][k]);
							}

				// ������ ����� ������ �����
				err = WriteBjnGzippedScalar8RecFuncByBlock(targetgrid, dname, F1, def.locNx*(def.rank), 0, 0, def.locNx, def.locNy, def.locNz, 9); // ��������� def.locNx*(def.rank) �� ������ ��������
				if(err)
					fprintf(stderr, "Can't add grid `%s` block data err %d\n", dname, err);
				}

		// ������ ������������� � ������������ �������� �����
		err = WriteBjnGzippedScalar8RecFuncMinMax(targetgrid, dname, fmin, fmax);
		if(err) 
			fprintf(stderr, "Can't add data about block err %d\n", err);
	}
*/
}

// ������� �������� ����� ��������������
void load_permeability(double* K, consts def)
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

	if((Nx!=def.locNx) || (Ny!=def.locNy))
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

			for (int k=0;k<def.locNz;k++)
				K[i+j*def.locNx+k*def.locNx*def.locNy]=1e-10 * exp(atof(value));
		}
	}

	fclose(input);
}

// ���������� ��������� � ����
void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, consts def)
{
	// ��������� � ������ ����� ������ �� roS_old
	// P1 � S2 ��������� ��� ��� ������� ���������� ��������,
	// x,y � media �� ���������� � �������� �� �������.
	//load_data_to_host(HostArraysPtr.P1, DevArraysPtr.P1 , localNx);
	//load_data_to_host(HostArraysPtr.S2, DevArraysPtr.S2 , localNx);
	load_data_to_host(HostArraysPtr.roS_w_old, DevArraysPtr.roS_w_old , def);
	load_data_to_host(HostArraysPtr.roS_n_old, DevArraysPtr.roS_n_old , def);
#ifdef THREE_PHASE
	load_data_to_host(HostArraysPtr.roS_g_old, DevArraysPtr.roS_n_old , def);
#endif

	FILE *f_save;

	if ((def.rank)==0)
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

	for (int cpu=0; cpu<(def.sizex * (def.sizey) * (def.sizez));cpu++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
		barrier();
		if ((def.rank)==cpu)
		{
			if(!(f_save=fopen("save/save.dat","ab")))
			{
				printf("\nError: Not open file \"save.dat\"!\n");
				exit(0);
			}
			fwrite(&j, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.S_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#endif
			fwrite(HostArraysPtr.P_w, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.S_n, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.x, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.y, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			//fwrite(HostArraysPtr.z, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_w_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fwrite(HostArraysPtr.roS_n_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#ifdef THREE_PHASE
			fwrite(HostArraysPtr.roS_g_old, sizeof(double), (def.locNx) * (def.locNy) * (def.locNz), f_save);
#endif
			fwrite(HostArraysPtr.media, sizeof(int), (def.locNx) * (def.locNy) * (def.locNz), f_save);
			fclose(f_save);
		}
	}
}

// �������������� ��������� �� �����
void restore (ptr_Arrays HostArraysPtr, int* j, consts def)
{
	FILE *f_save;
	for (int cpu=0; cpu<(def.sizex * (def.sizey) * (def.sizez));cpu++)
	{
		// ���������� ������ Barrier ��� ��������� ������������
		barrier();

		consts def_tmp;
		def_tmp.locNx = 0;
		def_tmp.locNy = 0;
		def_tmp.locNz = 0;
		def_tmp.Nx = def.Nx;
		def_tmp.Ny = def.Ny;
		def_tmp.Nz = def.Nz;
		def_tmp.size = def.size;
		def_tmp.sizex = def.sizex;
		def_tmp.sizey = def.sizey;
		def_tmp.sizez = def.sizez;

		if ((def.rank)==cpu)
		{
			if(!(f_save=fopen("save/save.dat","rb")))
			{
					printf("\nError: Not open file \"save.dat\"!\n");
					exit(0);
			}
			for (int queue=0;queue<=(def.rank);queue++)
			{
				def_tmp.rank = queue;
				global_to_local_vars(&def_tmp); 
				fread(j, sizeof(int), 1, f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.S_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#endif
				fread(HostArraysPtr.P_w, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.S_n, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.x, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.y, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				//fread(HostArraysPtr.z, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_w_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
				fread(HostArraysPtr.roS_n_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#ifdef THREE_PHASE
				fread(HostArraysPtr.roS_g_old, sizeof(double), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
#endif
				fread(HostArraysPtr.media, sizeof(int), (def_tmp.locNx) * (def_tmp.locNy) * (def_tmp.locNy), f_save);
			}
			fclose(f_save);
		}
	}
}



//------------------------------------------------------------------------------------------

// �������� � �������� ������ ����������� ���������
void read_version(void)
{

	FILE *rev;
	char str[250]="";
	int revision;

	if(!(rev=fopen("../.svn/entries","rt")))
		revision=0;
	else
	{
		for(int i=0;i<4;i++)
			fgets(str,250,rev);
		revision=atoi(str);
	}

	printf("Version %s.%d compiled %s %s.\n\n", VERSION, revision, __DATE__, __TIME__);
}

// ���������� ���������� ������ �� �����
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

#ifdef B_L
		if(!strcmp(attr_name,"Q")) 
		{(*def).Q = atof(attr_value); continue;}
#endif

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
		if(!strcmp(attr_name,"S_W_INIT")) 
			{(*def).S_w_init = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_N_INIT")) 
			{(*def).S_n_init = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_NR_0")) 
			{(*def).S_nr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_NR_1")) 
			{(*def).S_nr[1] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_GR_0")) 
			{(*def).S_gr[0] = atof(attr_value); continue;}
		if(!strcmp(attr_name,"S_GR_1")) 
			{(*def).S_gr[1] = atof(attr_value); continue;}
#endif
		if(!strcmp(attr_name,"S_N_GR")) 
			{(*def).S_n_gr = atof(attr_value); continue;}
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
