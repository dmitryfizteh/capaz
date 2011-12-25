#ifndef DEFINES_H
#define DEFINES_H 

#define VERSION "0.4"

#ifdef TWO_PHASE
//���� 1, �� ������� ����������� ��
#define TWO_LAYERS 0
// ������� ������������� ���������� (���� �� ���������� NR, �� ��������� ��� ���)
#define NR 
#else
#define TWO_LAYERS 0
#define NR 
#endif

// ���������� ���������������� �� ���� ��������
// ��� �-100 - 3, ��� ���-�������� 1 ��� 2
#define GPU_PER_NODE 3
#define DEFINES_FILE "..//defines.ini"

// ����� � ����� ����������
#define BlockNX 8
#define BlockNY 8
#define BlockNZ 8

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

#ifdef THREE_PHASE	
	// ������������ ������, ������������ ������� ����������� �������� �� �������� ��������� ��������� ������������� [0,1]
	const double aw[2] = {-39640, -26530};
	const double bw[2] = {12560, 27110};
	const double ag[2] = {14170, 35380};
	const double bg[2] = {158, -2969};
	// �������� ����� ��������� ��������� ������������� ��� ���������� �������������� � ����������� ��������
	const double S_w_range[3] = {0.001, 0.1, 0.99};
	const double S_g_range[3] = {0.001, 0.005, 0.95};
#endif

struct ptr_Arrays_tag 
{
	double *P_w, *P_n, *ro_w, *ro_n, *ux_w, *uy_w, *uz_w, *ux_n, *uy_n, *uz_n, *Xi_w, *Xi_n,*roS_w,*roS_w_old,*roS_n,*roS_n_old;
	int *media;
#ifdef THREE_PHASE 
	double *P_g, *S_w, *S_g, *ro_g, *ux_g, *uy_g, *uz_g, *Xi_g, *roS_g, *roS_g_old;
#else
	double *S_n;
#endif
#ifdef B_L
	double *K;
#endif
};
typedef struct ptr_Arrays_tag ptr_Arrays;

struct parts_sizes
{
	int x, y, z;
};
typedef struct parts_sizes parts_sizes;

struct localN
{
	int x, y, z;
};
typedef struct localN localN;

// ��������� ���������� ����
struct consts_tag
{
	double lambda[2];
	double S_wr[2];
	double m[2];
	double hx, hy, hz, dt, tau, l, c_w, c_n, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, timeX;
	int Nx, Ny, Nz;
	int source, save_plots, print_screen, newton_iterations;
#ifndef B_L
	double K[2];
#endif
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

extern void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t, localN locN, int rank,parts_sizes parts, int blocksX, int blocksY, int blocksZ);
extern void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, int* j, localN* locN, int* size, parts_sizes *parts, int* rank, int* blocksX, int* blocksY, int* blocksZ, int argc, char* argv[], consts def);
extern void load_permeability(double* K, localN locN);
extern void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, localN locN, consts def);
extern void host_memory_allocation(ptr_Arrays* ArraysPtr, localN locN, consts def);
extern void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, localN locN, consts def);
extern void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr);
extern void host_memory_free(ptr_Arrays HostArraysPtr);
extern void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, parts_sizes parts, int rank, localN locN, consts def);
extern void data_initialization(ptr_Arrays HostArraysPtr, int* t, localN locN, int rank, parts_sizes parts, consts def);
extern void parts_initialization(int size, parts_sizes *parts);
extern void communication_initialization(int argc, char* argv[], int* size, int* rank, consts def);
extern void communication_finalization(void);
extern void global_to_local_vars (localN* locN, parts_sizes parts, int rank, consts def);
extern int is_active_point(int i, int j, int k, localN locN, int rank, parts_sizes parts);
extern void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def);
extern void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def);
extern void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, localN locN, consts def);

extern void device_initialization(int rank, int* blocksX, int* blocksY, int* blocksZ, localN locN, consts def);
extern void device__finalization(void);
//extern void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def);

// ���������
extern void print_plots_top (double t, consts def);
extern void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, parts_sizes parts, localN locN, consts def);
extern void barrier(void);
extern void restore (ptr_Arrays HostArraysPtr, int* j, int rank, parts_sizes parts, localN locN, consts def);
extern void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int j, int rank, parts_sizes parts, localN locN, consts def);
extern void read_defines(int argc, char *argv[], consts* def);
extern void read_version(void);

// Unit-�����
#ifndef THREE_PHASE
extern void test_correct_P_S(ptr_Arrays HostArraysPtr, localN locN, int rank, consts def);
#endif

extern void test_nan (double x, char *file, int line);
extern void test_positive (double x, char *file, int line);
extern void test_S (double S, char *file, int line);
extern void test_arrowhead(double big, double small, char *file, int line);
extern void test_tau(double S_old, double S_now, double S_new, int media, consts def, char *file, int line);
extern void read_defines_test(consts def);


// ������� � ������ �����
extern int i_to_I(int i, int rank, parts_sizes parts, consts def);
extern double ro_eff_gdy(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);
extern void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);
extern void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);
extern void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);
extern void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, localN locN, consts def);
extern void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, localN locN, consts def);
extern void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);
extern void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, int rank, parts_sizes parts, consts def);
extern void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def);

extern void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ);
extern void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def);
extern void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def);
extern void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def);
extern void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ);
extern void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ);
extern void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def);
extern void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def);

extern void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def);
extern void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def);

#endif