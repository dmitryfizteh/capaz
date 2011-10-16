#ifndef DEFINES_H
#define DEFINES_H 

#define TWO_PHASE
// Считать направленными разностями (если не определена NR, то считается без них)
//#define NR 
#define TEST
#define DEFINES_FILE "..//defines.ini"

// Нитей в блоке ускорителя
#define BlockNX 4
#define BlockNY 4
#define BlockNZ 4

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

const double K[2]={6.64e-11,7.15e-12};
const double lambda[2]={2.7,2.0};
const double S_wr[2]={0.09,0.12};
const double m[2]={0.4,0.39};
const double P_d[2]={755,2060};

struct ptr_Arrays_tag 
{
	double *x, *y, *z;
	double *P_w, *P_n, *S_n, *ro_w, *ro_n, *ux_w, *uy_w, *uz_w, *ux_n, *uy_n, *uz_n, *Xi_w, *Xi_n,*roS_w,*roS_w_old,*roS_n,*roS_n_old;
	int *media;
};
typedef struct ptr_Arrays_tag ptr_Arrays;

// Структура параметров сред
struct consts_tag
{
	double K[2];
	double lambda[2];
	double S_wr[2];
	double m[2];
	double P_d[2];
	double hx, hy, hz, dt, tau, l_w, l_n, c, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, S_n_gr;
	int Nx, Ny, Nz;
	int source, timeX, save_plots, print_screen, newton_iterations;
};
typedef struct consts_tag consts;

extern void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t, int localNx, int localNy, int rank,int size, int blocksX, int blocksY, int blocksZ);
extern void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int* j, int* localNx, int* localNy, int* size, int* rank, int* blocksX, int* blocksY, int* blocksZ, int argc, char* argv[], consts def);
extern void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int localNx, int nY, consts def);
extern void host_memory_allocation(ptr_Arrays* ArraysPtr, int nX, int nY, consts def);
extern void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, int nX, consts def);
extern void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr);
extern void host_memory_free(ptr_Arrays HostArraysPtr);
extern void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, int size, int rank, int localNx, consts def);
extern void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def);
extern void communication_initialization(int argc, char* argv[], int* size, int* rank, consts def);
extern void communication_finalization(void);
extern void global_to_local_vars (int* localNx, int* localNy, int size, int rank, consts def);
extern int is_active_point(int i, int localNx, int rank, int size);
extern void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, int localNx, consts def);
extern void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, int localNx, consts def);
extern void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, int localNx, consts def);

extern void device_initialization(int rank, int* blocksX, int* blocksY, int* blocksZ, int localNx, consts def);

// Служебные
extern void print_plots_top (double t, consts def);
extern void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, int size, int localNx, consts def);
extern void barrier(void);
extern void restore (ptr_Arrays HostArraysPtr, int* j, int rank, int size, int localNx, consts def);
extern void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, int rank, int size, int localNx, consts def);
extern void read_defines(int argc, char *argv[], consts* def);

// Unit-тесты
extern void test_correct_Pw_Sn(ptr_Arrays HostArraysPtr, int nX, int rank, consts def);

// Расчеты в каждой точке
extern int i_to_I(int i, int rank, int size, consts def);
extern void assign_Pn_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def);
extern void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def);
extern void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void Border_Sn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def);
extern void Border_Pw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);

extern void ro_Pn_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void Pn_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def);
extern void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void Pw_Sn_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def);
extern void Pw_Sn_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);

extern void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);

#endif