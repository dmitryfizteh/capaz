#ifndef DEFINES_H
#define DEFINES_H 

// Считать направленными разностями (если не определена NR, то считается без них)
//#define NR 

#ifdef TWO_PHASE
//Если 1, то считать двухслойной РС
#define TWO_LAYERS 1
#else
#define TWO_LAYERS 0
#endif

#define DEFINES_FILE "..//defines.ini"

// Нитей в блоке ускорителя
#define BlockNX 4
#define BlockNY 4
#define BlockNZ 4

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

#ifndef THREE_PHASE
	const double K[2]={6.64e-11,7.15e-12};
	const double lambda[2]={2.7,2.0};
	const double S_wr[2]={0.09,0.12};
	const double m[2]={0.4,0.39};
	const double P_d[2]={755,2060};
#else
	// Коэффициенты прямых, продолжающих функции капиллярных давлений на границах интервала изменения насыщенностей [0,1]
	const double aw[2] = {-39640, -26530};
	const double bw[2] = {12560, 27110};
	const double ag[2] = {14170, 35380};
	const double bg[2] = {158, -2969};
	// Ключевые точки интервала изменения насыщенностей для вычисления проницаемостей и капиллярных давлений
	const double S_w_range[3] = {0.001, 0.1, 0.99};
	const double S_g_range[3] = {0.001, 0.005, 0.95};
#endif

struct ptr_Arrays_tag 
{
	double *x, *y, *z;
	double *P_w, *P_n, *ro_w, *ro_n, *ux_w, *uy_w, *uz_w, *ux_n, *uy_n, *uz_n, *Xi_w, *Xi_n,*roS_w,*roS_w_old,*roS_n,*roS_n_old;
	int *media;
#ifdef THREE_PHASE 
	double *P_g, *S_w, *S_g, *ro_g, *ux_g, *uy_g, *uz_g, *Xi_g, *roS_g, *roS_g_old;
#else
	double *S_n;
#endif
};
typedef struct ptr_Arrays_tag ptr_Arrays;

// Структура параметров сред
struct consts_tag
{
	double K[2];
	double lambda[2];
	double S_wr[2];
	double m[2];
	double hx, hy, hz, dt, tau, l, c_w, c_n, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, timeX;
	int Nx, Ny, Nz;
	int source, save_plots, print_screen, newton_iterations;
#ifdef TWO_PHASE
	double P_d[2];
#endif
#ifdef THREE_PHASE
	double S_w_init, S_g_init;
	double S_nr[2];
	double S_gr[2];
	double P_d_nw[2];
	double P_d_gn[2];
	double c_g, beta_g, mu_g, ro0_g, S_w_gr, S_g_gr;
	double aw[2], bw[2], ag[2], bg[2];
	double S_w_range[3];
	double S_g_range[3];
#else
	double S_n_gr;
#endif
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
extern void device__finalization(void);
extern void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def);

// Служебные
extern void print_plots_top (double t, consts def);
extern void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, int size, int localNx, consts def);
extern void barrier(void);
extern void restore (ptr_Arrays HostArraysPtr, int* j, int rank, int size, int localNx, consts def);
extern void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, int rank, int size, int localNx, consts def);
extern void read_defines(int argc, char *argv[], consts* def);

// Unit-тесты
#ifndef THREE_PHASE
extern void test_correct_P_S(ptr_Arrays HostArraysPtr, int nX, int rank, consts def);
#endif

extern void test_nan (double x, char *file, int line);
extern void test_positive (double x, char *file, int line);
extern void test_arrowhead(double big, double small, char *file, int line);
extern void read_defines_test(consts def);

// Расчеты в каждой точке
extern int i_to_I(int i, int rank, int size, consts def);
extern void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def);
extern void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def);
extern void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);
extern void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def);
extern void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);

extern void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def);
extern void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ);
extern void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def);
extern void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);

extern void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);
extern void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def);

#endif