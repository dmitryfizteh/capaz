#include "../defines.h"
#include "three-phase.h"

//Функция вычисления значений давлений, плотностей и коэффициентов в законе Дарси в точке (i,j,k) среды media,
//исходя из известных значений основных параметров (Pw,Sw,Sn)
//1. Запоминаем, с какой именно из сред работаем
//2. Вычисление значения насыщенности фазы n из условия равенства насыщенностей в сумме единице
//3. Вычисление эффективных насыщенностей по формулам модели трехфазной фильтрации
//4. Вычисление относительных фазовых проницаемостей в соответствии с приближением Стоуна в модификации Азиза и Сеттари
//5. Вычисление капиллярных давлений в соответствии с приближенной моделью Паркера
//6. Вычисление фазовых давлений c помощью капиллярных
//7. Вычисление коэффициентов закона Дарси
void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{

	int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];                                                                                                                                                                            
	double k_w, k_g, k_n, P_k_nw, P_k_gn;
	double A = def.lambda[media]; 
	double S_w_e = (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);                                                                                                                            
	double S_n_e = (HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_nr[media]) / (1. - def.S_nr[media] - def.S_nr[media] - def.S_nr[media]); 
	double S_g_e = 1. - S_w_e - S_n_e;

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

	if(S_n_e <= 0.)
	{
		k_n = 0.;
	}
	else
	{
		double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);    
		double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A /  (A - 1.)), (A - 1.) / A), 2.);
		k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e); 
	}
	
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

	HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw;                                                             
	HostArraysPtr.P_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_gn;                                                             

	HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_w / def.mu_w;                                                                                                                 
	HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_n / def.mu_n;                                                                                                                 
	HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_g / def.mu_g;                                                                                                                 

	test_positive(HostArraysPtr.P_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_g[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
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
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double S_w_e, S_g_e, S_n_e, P_k_nw, P_k_gn, A, Sg, F1, F2, F3;
		double PkSw, PkSn, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sn, F2Sn, F3Sn, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);   
			S_n_e = (HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_nr[media]) / (1. - def.S_nr[media] - def.S_nr[media] - def.S_nr[media]); 
			S_g_e = 1. - S_w_e - S_n_e;
			A = def.lambda[media];                                                                                                                                                                                                                                                                  /*2*/

			// По краям интервала [0, 1] функции капиллярных давлений гладко заменяем линейными, производные меняются соответственно. 
			// Описание можно посмотреть в файле mathcad.
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
				PkSn = (-1) * (def.ag[0]);
			}
			else if(S_g_e >= def.S_g_range[2])
			{
				P_k_gn = (def.ag[1]) * S_g_e + (def.bg[1]);
				PkSn = (-1) * (def.ag[1]);
			}
			else
			{
				P_k_gn = def.P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A); 
				PkSn = def.P_d_gn[media] * pow(pow(1. - S_g_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S_g_e, A / (1. - A) - 1.) / (1. - A)           
					/(1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);  
			}   

			Sg = 1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];                     

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm))                                                               
				* HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];                  
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw - def.P_atm))
				* HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.roS_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
			F3 = def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn) / def.P_atm                                                              
				* Sg - HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

			F1P = def.ro0_w * def.beta_w * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];                                                                                                                             
			F2P = def.ro0_n * def.beta_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];                                                                                                                                                                                                                                    
			F3P = def.ro0_g * Sg;   

			F1Sw = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));                                                                                                                                                         
			F2Sn = def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw - def.P_atm));                                         
			F2Sw = F1Sn = 0;                                                                                                                                                                                                                                                                     
			F3Sn = (-1) * def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn - Sg * PkSn) / def.P_atm;
			F3Sw = (-1) * def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn - Sg * (PkSn + PkSw)) / def.P_atm;

			det = F1P * F2Sw * F3Sn - F1Sw * (F2P * F3Sn - F2Sn * F3P);                                                                                                                                                                                         

			HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]                                                              
			- (1. / det) * (F2Sw * F3Sn * F1 - F1Sw * F3Sn* F2 + F1Sw * F2Sn * F3);
			HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]                                                              
			- (1. / det) * ((F2Sn * F3P - F2P * F3Sn) * F1 + F1P * F3Sn * F2 - F1P * F2Sn * F3);
			HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]                                                              
			- (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);
		}  
		
		test_S(HostArraysPtr.S_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	}
}

//Задание граничных условий отдельно для (Sw,Sg),Pn

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (def.locNx) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (def.locNy) - 1)
		j1 --;
	if((k == 0) && ((def.locNz) > 2))
		k1 ++;
	if((k == (def.locNz) - 1) && ((def.locNz) > 2))
		k1 --;

	if((j != 0) || ((def.source) <= 0))
	{
		HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}

	if((j == 0) && ((def.source) > 0))
	{
		HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.S_w_gr;
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.S_n_gr;
	}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (def.locNx) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (def.locNy) - 1)
		j1 --;
	if((k == 0) && ((def.locNz) > 2))
		k1 ++;
	if((k == (def.locNz) - 1) && ((def.locNz) > 2))
		k1 --;

	if((j != 0) && (j != (def.locNy) - 1))
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	else if(j == 0)
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	else
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] + ro_eff_gdy(HostArraysPtr, i1, j1, k1, def);
}

void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for(int i = 0; i < def.locNx; i++)
		for(int j = 0; j < def.locNy; j++)
			for(int k = 0; k < def.locNz; k++)
				if(is_active_point(i, j, k, def))
					{
						// Преобразование локальных координат процессора к глобальным
						int I = local_to_global(i, 'x', def);

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

						///!!!! Не учитываются каппилярные силы! Или надо считать перед этим шагом P_w, P_g
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
						HostArraysPtr.ro_g[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));	

					test_nan(HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
#ifdef THREE_PHASE 
					test_nan(HostArraysPtr.S_w[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
#endif
					}
}
