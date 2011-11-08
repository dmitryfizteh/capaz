#include "../defines.h"
#include "three-phase.h"

//Функция вычисления значений давлений, плотностей и коэффициентов в законе Дарси в точке (i,j,k) среды media,
//исходя из известных значений основных параметров (Pn,Sw,Sg)
//1. Запоминаем, с какой именно из сред работаем
//2. Вычисление значения насыщенности фазы n из условия равенства насыщенностей в сумме единице
//3. Вычисление эффективных насыщенностей по формулам модели трехфазной фильтрации
//4. Вычисление относительных фазовых проницаемостей в соответствии с приближением Стоуна в модификации Азиза и Сеттари
//5. Вычисление капиллярных давлений в соответствии с приближенной моделью Паркера
//6. Вычисление фазовых давлений c помощью капиллярных
//7. Вычисление коэффициентов закона Дарси
void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{

	int media = HostArraysPtr.media[i + j * localNx + k * localNx * (def.Ny)];                                                                                                                                                                            
	double k_w, k_g, k_n, P_k_nw, P_k_gn;
	double A = def.lambda[media]; 
	double S_w_e = (HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);                                                                                                                            
	double S_g_e = (HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]); 
	double S_n_e = 1. - S_w_e - S_g_e;

	if(S_w_e <= def.S_w_range[0])
	{
		S_w_e = def.S_w_range[0];
		k_w = 0.;
	}
	else
		k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)),(A - 1.) / A), 2.);

	if(S_g_e <= def.S_g_range[0])
	{
		S_g_e = def.S_g_range[0];
		k_g = 0.;
	}
	else
		k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);

	double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);    
	double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A /  (A - 1.)), (A - 1.) / A), 2.);
	k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e); 

	if(S_w_e <= def.S_w_range[1])
		P_k_nw = (def.aw[0]) * S_w_e + (def.bw[0]); 
	else if(S_w_e >= def.S_w_range[2])
		P_k_nw = (def.aw[1]) * S_w_e + (def.bw[1]); 
	else
		P_k_nw = def.P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A); 
	
	if(S_g_e <= def.S_g_range[1])
		P_k_gn = (def.ag[0]) * S_g_e + (def.bg[0]); 
	else if(S_g_e >= def.S_g_range[2])
		P_k_gn = (def.ag[1]) * S_g_e + (def.bg[1]);    
	else
		P_k_gn = def.P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A);                                                           

	HostArraysPtr.P_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw;                                                             
	HostArraysPtr.P_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn;                                                             

	HostArraysPtr.Xi_w[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_w / def.mu_w;                                                                                                                 
	HostArraysPtr.Xi_n[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_n / def.mu_n;                                                                                                                 
	HostArraysPtr.Xi_g[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_g / def.mu_g;                                                                                                                 

	test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
}

//Функция решения системы 3*3 на основные параметры (Pn,Sw,Sg) методом Ньютона в точке (i,j,k) среды media
//1. Вычисление эффективных насыщенностей
//2. Переобозначение степенного коэффициента
//3. Вычисление капиллярных давлений
//4. Вычисление насыщенности фазы n
//5. Нахождение значения трех функций системы
//6. Вычисление частных производных производных капиллярных давлений по насыщенностям
//7. Вычисление матрицы частных производных
//8. Вычисление детерминанта матрицы частных производных
//9. Получение решения системы методом Крамера в явном виде
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i != 0) && (i != localNx - 1) && (j != 0) && (j != (def.Ny) - 1) && (((k != 0) && (k != (def.Nz) - 1)) || ((def.Nz) < 2)))
	{
		int media = HostArraysPtr.media[i + j * localNx + k * localNx * (def.Ny)];
		double S_w_e, S_g_e, P_k_nw, P_k_gn, A, Sn, F1, F2, F3;
		double PkSw, PkSg, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sg, F2Sg, F3Sg, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = (HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);   
			S_g_e = (HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);  

			A = def.lambda[media];                                                                                                                                                                                                                                                                  /*2*/

			if(S_w_e <= def.S_w_range[1])
			{
				P_k_nw = (def.aw[0]) * S_w_e + (def.bw[0]);
				PkSw = (def.aw[0]);
			}
			else if(S_w_e >= def.S_w_range[2])
			{
				P_k_nw = (def.aw[1]) * S_w_e + (def.bw[1]);
				PkSw = (def.aw[1]);
			}
			else
			{
				P_k_nw = def.P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A); 
				PkSw = def.P_d_nw[media] * pow(pow(S_w_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(S_w_e, (A / (1. - A) - 1.)) / (1. - A)                                                                    
					/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);   
			}

			if(S_g_e <= def.S_g_range[1])
			{
				P_k_gn = (def.ag[0]) * S_g_e + (def.bg[0]);
				PkSg = (def.ag[0]);
			}
			else if(S_g_e >= def.S_g_range[2])
			{
				P_k_gn = (def.ag[1]) * S_g_e + (def.bg[1]);
				PkSg = (def.ag[1]);
			}
			else
			{
				P_k_gn = def.P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A); 
				PkSg = (-1) * def.P_d_gn[media] * pow(pow(1. - S_g_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S_g_e, A / (1. - A) - 1.) / (1. - A)           
					/(1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);  
			}   

			Sn = 1. - HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)];                     

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw - def.P_atm))                                                               
				* HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.roS_w[i + j * localNx + k * localNx * (def.Ny)];                  
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - def.P_atm)) * Sn                                                                  
				- HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
			F3 = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn - def.P_atm))                                                              
				* HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.roS_g[i + j * localNx + k * localNx * (def.Ny)];

			F1P = def.ro0_w * def.beta_w * HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)];                                                                                                                             
			F2P = def.ro0_n * def.beta_n * Sn;                                                                                                                                                                                                                                    
			F3P = def.ro0_g * def.beta_g * HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)];                                                                                                                                                
			F1Sw = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw - def.P_atm
				- HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] * PkSw));                                                                                                                                                         
			F2Sw = F2Sg = (-1) * def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - def.P_atm));                                         
			F3Sw = F1Sg = 0;                                                                                                                                                                                                                                                                     
			F3Sg = def.ro0_g * (1 + def.beta_g * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn - def.P_atm                                                                  
				+ HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] * PkSg));

			det = F1P * F2Sw * F3Sg - F1Sw * (F2P * F3Sg - F2Sg * F3P);                                                                                                                                                                                         

			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)]                                                              
			- (1. / det) * (F2Sw * F3Sg * F1 - F1Sw * F3Sg * F2 + F1Sw * F2Sg * F3);
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)]                                                              
			- (1. / det) * ((F2Sg * F3P - F2P * F3Sg) * F1 + F1P * F3Sg * F2 - F1P * F2Sg * F3);
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)]                                                              
			- (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);

			test_nan(HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			test_nan(HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		}  
	}
}

//Задание граничных условий отдельно для (Sw,Sg),Pn
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((j == 0) && ((def.Ny) > 2))
	{
			if((i != 0) && (i != localNx - 1))
			{
				HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + (j + 1) * localNx + k * localNx * (def.Ny)];
				HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + (j + 1) * localNx + k * localNx * (def.Ny)];
			}
			else if(i == 0)
			{
				HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)];
				HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)];
			}
			else
			{
				HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)];
				HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)];
			}
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = def.S_w_gr;
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = def.S_g_gr;
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny) > 2))
	{
		if ((i != 0) && (i != localNx - 1))
		{
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)];
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)];
		}
		else if(i == 0)
		{
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)];
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)];
		}
		else
		{
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)];
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)];
		}

		return;
	}


	if ((i == 0) && ((def.Nx) > 2) && (j != 0) && (j != (def.Ny) - 1))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + 1 + j * localNx + k * localNx * (def.Ny)];
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx) > 2) && (j != 0) && (j != (def.Ny) - 1))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i - 1 + j * localNx + k * localNx * (def.Ny)];
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i - 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((k == 0) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + (k + 1) * localNx * (def.Ny)];
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + (k + 1) * localNx * (def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		return;
	}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((j == 0) && ((def.Ny) > 2))
	{
		//Если верхний слой - это открытая поверхность
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = def.P_atm;
/*
		//Если ставим условие только на градиент давления.
		if((i != 0) && (i != localNx - 1))
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + (j + 1) * localNx + k * localNx *(def.Ny)] 
			- (HostArraysPtr.ro_n[i + (j + 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i + (j + 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + (j + 1) * localNx + k * localNx * (def.Ny)]) 
				+ HostArraysPtr.ro_w[i + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i + (j + 1) * localNx + k * localNx * (def.Ny)]
			+ HostArraysPtr.ro_g[i + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i + (j + 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}
		else if(i == 0)
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + 1 + (j + 1) * localNx + k * localNx *(def.Ny)] 
			- (HostArraysPtr.ro_n[i + 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)]) 
				+ HostArraysPtr.ro_w[i + 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)]
			+ HostArraysPtr.ro_g[i + 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i + 1 + (j + 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}
		else
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i - 1 + (j + 1) * localNx + k * localNx *(def.Ny)] 
			- (HostArraysPtr.ro_n[i - 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)]) 
				+ HostArraysPtr.ro_w[i - 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)]
			+ HostArraysPtr.ro_g[i - 1 + (j + 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i - 1 + (j + 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}
*/
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny) > 2))
	{
		// Если на нижнем слое давление фиксированное
//		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = 1.054e+5;
		//HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = (def.P_atm) * ((def.ro0_w) * (def.S_w_init) + (def.ro0_n) * (1. - (def.S_w_init) - (def.S_g_init)) + (def.ro0_g) * (def.S_g_init)) * (def.g_const) * (def.hy) * (def.Ny);
		// Если на нижних двух слоях давление одинаковое
/*		if ((i != 0) && (i != localNx - 1))
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + (j - 1) * localNx + k * localNx *(def.Ny)];
		else if(i == 0)
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + 1 + (j - 1) * localNx + k * localNx *(def.Ny)]; 
		else
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i - 1 + (j - 1) * localNx + k * localNx *(def.Ny)];
*/
		// Если на нижней границе есть градиент давления
		if ((i != 0) && (i != localNx - 1))
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + (j - 1) * localNx + k * localNx *(def.Ny)] 
		+ (HostArraysPtr.ro_n[i + (j - 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)]) 
			+ HostArraysPtr.ro_w[i + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)]
		+ HostArraysPtr.ro_g[i + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}
		else if(i == 0)
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + 1 + (j - 1) * localNx + k * localNx *(def.Ny)] 
			+ (HostArraysPtr.ro_n[i + 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)]) 
				+ HostArraysPtr.ro_w[i + 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)]
			+ HostArraysPtr.ro_g[i + 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i + 1 + (j - 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}
		else
		{
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i - 1 + (j - 1) * localNx + k * localNx *(def.Ny)] 
			+ (HostArraysPtr.ro_n[i - 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * (1. - HostArraysPtr.S_w[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)]) 
				+ HostArraysPtr.ro_w[i - 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_w[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)]
			+ HostArraysPtr.ro_g[i - 1 + (j - 1) * localNx + k * localNx *(def.Ny)] * HostArraysPtr.S_g[i - 1 + (j - 1) * localNx + k * localNx * (def.Ny)]) * (def.g_const) * (def.hy);
		}

		return;
	}
	
	if ((i == 0) && (j != 0) && (j != (def.Ny) - 1))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + 1 + j * localNx + k * localNx * (def.Ny)]; 
		return;
	}

	if ((i == localNx - 1) && (j != 0) && (j != (def.Ny) - 1))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i - 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((k == 0) && ((def.Nz) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + (k + 1) * localNx * (def.Ny)]; 
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		return;
	}
}