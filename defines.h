#ifndef DEFINES_H
#define DEFINES_H

//#define GTEST

#ifdef GTEST
#include "gtest/gtest.h"
#endif

#define VERSION "0.52"

#ifdef TWO_PHASE
//Если 1, то считать двухслойной РС
#define TWO_LAYERS 0
// Считать направленными разностями (если не определена NR, то считается без них)
#define NR
#else
#define TWO_LAYERS 1
#define NR
#endif

// Количество видеоускорителей на узле кластера
// Для К-100 - 3, для МВС-Экспресс 1 или 2
#define GPU_PER_NODE 3
#define DEFINES_FILE "..//defines.ini"

// Нитей в блоке ускорителя
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

// Псевдо-функция минимума/максимума
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

// Вывод графиков BjnIO
//#include "bjnio.h"

struct ptr_Arrays_tag
{
	double *S_n, *P_w, *P_n, *ro_w, *ro_n, *ux_w, *uy_w, *uz_w, *ux_n, *uy_n, *uz_n, *Xi_w, *Xi_n, *roS_w, *roS_w_old, *roS_n, *roS_n_old;
	double *m, *K;
#ifdef THREE_PHASE
	double *P_g, *S_w, *ro_g, *ux_g, *uy_g, *uz_g, *Xi_g, *roS_g, *roS_g_old;
#endif
};
typedef struct ptr_Arrays_tag ptr_Arrays;

// Структура параметров сред и характерных размеров задачи
struct consts_tag
{
	double upscale_l, upscale_t;
	double lambda[2];
	double S_wr[2];
	double porosity[2];
	double hx, hy, hz, dt, tau, l, c_w, c_n, beta_w, beta_n, P_atm, g_const, mu_w, mu_n, ro0_w, ro0_n, timeX;
	double S_n_gr;
	int Nx, Ny, Nz;
	int source, save_plots, print_screen, newton_iterations;
	double K[2];
	double Q, InjWell_Pw, InjWell_Sn, OutWell_Pw, OutWell_Sn, Background_Sn, Background_Pw; // Дебит скважины
#ifdef TWO_PHASE
	double P_d[2];
#endif
#ifdef THREE_PHASE
	double S_w_init, S_n_init;
	double S_nr[2];
	double S_gr[2];
	double P_d_nw[2];
	double P_d_gn[2];
	double c_g, beta_g, mu_g, ro0_g, S_w_gr;
#endif
	// Локальные размеры
	int locNx, locNy, locNz;
	// Число процессоров и ранг процессора
	int size, rank;
	// Число частей дробления области по измерениям
	int sizex, sizey, sizez;
	// Разложение ранга для трехмерного дробления области между процессорами
	int rankx, ranky, rankz;
	// Количество блоков ускорителя
	int blocksX, blocksY, blocksZ;
};
typedef struct consts_tag consts;

extern void time_step_function(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer, consts def, double t);
extern void initialization(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, long int* time_counter, int argc, char* argv[], consts* def);
extern void load_permeability(double* K, consts def);
extern void finalization(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void memory_allocation(ptr_Arrays* HostArraysPtr, ptr_Arrays* DevArraysPtr, consts def);
extern void host_memory_allocation(ptr_Arrays* ArraysPtr, consts def);
extern void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, consts def);
extern void memory_free(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr);
extern void host_memory_free(ptr_Arrays HostArraysPtr);
extern void device_memory_free(ptr_Arrays DevArraysPtr, double* DevBuffer);
extern void save_data_plots(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def);
extern void data_initialization(ptr_Arrays HostArraysPtr, long int* time_counter, consts def);
extern void sizes_initialization(consts* def);
extern void blocks_initialization(consts* def);
extern void communication_initialization(int argc, char* argv[], consts* def);
extern void communication_finalization(void);
extern void global_to_local_vars(consts* def);
extern int local_to_global(int local_index, char axis, consts def);
extern int is_active_point(int i, int j, int k, consts def);
extern void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, consts def);
extern void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, consts def);
extern void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, consts def);

extern void device_initialization(consts* def);
extern void device_finalization(void);

// Служебные
extern void print_plots_top(double t, consts def);
extern void print_plots(ptr_Arrays HostArraysPtr, double t, consts def);
extern void barrier(void);
extern void restore(ptr_Arrays HostArraysPtr, long int* time_counter, consts def);
extern void save(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, long int time_counter, consts def);
extern void read_defines(int argc, char *argv[], consts* def);
extern void read_version(void);
extern void print_error(char *error, char *file, int line);

// Unit-тесты
#ifndef THREE_PHASE
extern void test_correct_P_S(ptr_Arrays HostArraysPtr, consts def);
#endif

extern void test_nan(double x, char *file, int line);
extern void test_positive(double x, char *file, int line);
extern void test_S(double S, char *file, int line);
extern void test_u(double u, char *file, int line);
extern void test_arrowhead(double big, double small, char *file, int line);
extern void test_tau(double S_old, double S_now, double S_new, int local, consts def, char *file, int line);
extern void read_defines_test(consts def);


// Расчеты в каждой точке
extern double ro_eff_gdy(ptr_Arrays HostArraysPtr, int local, consts def);
extern void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, consts def);
extern void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, consts def);
extern void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def);
extern int set_boundary_basic_coordinate(int i, int j, int k, consts def);

extern void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def);
extern void P_ro_Xi_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def);
extern void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def);
extern void u_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def);
extern void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double t, consts def);
extern void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def);
extern void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def);
extern void P_S_exchange(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, double* HostBuffer, double* DevBuffer, consts def);

extern void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def);
extern void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, consts def);

extern int is_injection_well(int i, int j, int k, consts def);
extern int is_output_well(int i, int j, int k, consts def);
extern void assing_k(double* k_w, double* k_n, double S_w);

#endif

