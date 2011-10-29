#include "defines.h"

// Тестирование

// Функция проверки на выход из допустимого диапазона значений
// во всех точках расчетной области процессора
#ifndef THREE_PHASE
void test_correct_P_S(ptr_Arrays HostArraysPtr, int localNx, int rank, consts def)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
			{
				if (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]<0)
					printf ("\nWarning! S2<0 in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);
				if (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]<=0)
					printf ("\nWarning! P<=0 in point i=%d, j=%d, k=%d, rank=%d\n",i,j,k,rank);

				test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
				test_nan(HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			}
}
#endif

// Тест на NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_nan (double x, char *file, int line)
{
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
}

// Тест на положительное целое и не NaN
// Синтаксис вызова test_nan(x, __FILE__, __LINE__);
void test_positive (int x, char *file, int line)
{
	if ( isnan(x) )
		printf("Error: NaN\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
	if ( x < 0 )
		printf("Error: x<0\nFile:\"%s\"\nLine:\"%d\"\n\n", file, line);
}


void read_defines_test(consts def)
{
	test_positive(def.hx, __FILE__, __LINE__);

	test_positive(def.hy, __FILE__, __LINE__);

	test_positive(def.hz, __FILE__, __LINE__);

	test_positive(def.tau, __FILE__, __LINE__);
	
/*	
	if(!strcmp(attr_name,"DT")) 
	{(*def).dt = atof(attr_value); continue;}
	if(!strcmp(attr_name,"L_W")) 
	{(*def).l_w = atof(attr_value); continue;}
	if(!strcmp(attr_name,"L_N")) 
	{(*def).l_n = atof(attr_value); continue;}
	if(!strcmp(attr_name,"C")) 
	{(*def).c = atof(attr_value); continue;}
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
	if(!strcmp(attr_name,"L_G")) 
	{(*def).l_g = atof(attr_value); continue;}
	if(!strcmp(attr_name,"BETA_G")) 
	{(*def).beta_g = atof(attr_value); continue;}
	if(!strcmp(attr_name,"RO_G")) 
	{(*def).ro0_g = atof(attr_value); continue;}
	if(!strcmp(attr_name,"MU_G")) 
	{(*def).mu_g = atof(attr_value); continue;}
	if(!strcmp(attr_name,"S_W_GR")) 
	{(*def).S_w_gr = atof(attr_value); continue;}
	if(!strcmp(attr_name,"S_G_GR")) 
	{(*def).S_g_gr = atof(attr_value); continue;}
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
	*/
}