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

	int media = HostArraysPtr.media[i + j * localNx + k * localNx * (def.Ny)];																						/*1*/

	double Sn = (1. - HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)]);					/*2*/

	double S_w_e = (HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);	/*3*/
	double S_n_e = (Sn - def.S_nr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);																/*3*/
	double S_g_e = (HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);	/*3*/

	double k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, def.lambda[media] / (def.lambda[media] - 1.)),(def.lambda[media] - 1.) / def.lambda[media]), 2.);	/*4*/
	double k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, def.lambda[media] / (def.lambda[media] - 1.)), 2. * (def.lambda[media] - 1.) / def.lambda[media]);		/*4*/
	double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, def.lambda[media] / (def.lambda[media] - 1.)), 2. * (def.lambda[media] - 1.) / def.lambda[media]);	/*4*/
	double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, def.lambda[media] / (def.lambda[media] - 1.)), (def.lambda[media] - 1.) / def.lambda[media]), 2.);/*4*/
	double k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (S_w_e + S_n_e);																								/*4*/
	
	double P_k_nw = def.P_d_nw[media] * pow(pow(S_w_e, def.lambda[media] / (1. - def.lambda[media])) - 1., 1. / def.lambda[media]);									/*5*/
	double P_k_gn = def.P_d_nw[media] * pow(pow(1. - S_g_e, def.lambda[media] / (1. - def.lambda[media])) - 1., 1. / def.lambda[media]);							/*5*/

	HostArraysPtr.P_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw;								/*6*/
	HostArraysPtr.P_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn;								/*6*/
		
	HostArraysPtr.Xi_w[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_w / def.mu_w;															/*7*/
	HostArraysPtr.Xi_n[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_n / def.mu_n;															/*7*/
	HostArraysPtr.Xi_g[i + j * localNx + k * localNx * (def.Ny)] = (-1.) * (def.K[media]) * k_g / def.mu_g;															/*7*/

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
        	S_w_e = (HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);	/*1*/
			S_g_e = (HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);	/*1*/

        	A = def.lambda[media];																																	/*2*/

        	P_k_nw = def.P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A);																				/*3*/
			P_k_gn = def.P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A);																			/*3*/

        	Sn = (1. - HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)]);					/*4*/
        	
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw - def.P_atm))								/*5*/
				* HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.roS_w[i + j * localNx + k * localNx * (def.Ny)];			
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - def.P_atm)) * Sn									/*5*/
				- HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
			F3 = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn - def.P_atm))								/*5*/
				* HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] - HostArraysPtr.roS_g[i + j * localNx + k * localNx * (def.Ny)];
        	
			PkSw = def.P_d_nw[media] * pow(pow(S_w_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(S_w_e, (A / (1. - A) - 1.)) / (1. - A)									/*6*/					
				/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);		
			PkSg = (-1) * def.P_d_gn[media] * pow(pow(1. - S_g_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S_g_e, A / (1. - A) - 1.) / (1. - A)					/*6*/
				/(1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);	
        	
			F1P = def.ro0_w * def.beta_w * HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)];																/*7*/
			F2P = def.ro0_n * def.beta_n * Sn;																														/*7*/
			F3P = def.ro0_g * def.beta_g * HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)];																/*7*/			
			F1Sw = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - P_k_nw - def.P_atm
				- HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] * PkSw));																				/*7*/
			F2Sw = F2Sg = (-1) * def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] - def.P_atm));				     		/*7*/
			F3Sw = F1Sg = 0;																																		/*7*/
			F3Sg = def.ro0_g * (1 + def.beta_g * (HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] + P_k_gn - def.P_atm									/*7*/
				+ HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] * PkSg));
        	
			det = F1P * F2Sw * F3Sg - F1Sw * (F2P * F3Sg - F2Sg * F3P);																								/*8*/
        	
			HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)]								/*9*/
			- (1. / det) * (F2Sw * F3Sg * F1 - F1Sw * F3Sg * F2 + F1Sw * F2Sg * F3);
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)]								/*9*/
			- (1. / det) * ((F2Sg * F3P - F2P * F3Sg) * F1 + F1P * F3Sg * F2 - F1P * F2Sg * F3);
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)]								/*9*/
			- (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);

			test_nan(HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			test_nan(HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		}  
	}
}

//Вызов "персональных" граничных условий
void Border(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	Border_Sw(HostArraysPtr, i, j, k, localNx, rank, size, def);
	Border_Sg(HostArraysPtr, i, j, k, localNx, rank, size, def);
	Border_Pn(HostArraysPtr, i, j, k, localNx, def);
	return;
}

//Задание граничных условий отдельно для Sw,Sg,Pn
void Border_Sw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i - 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + (j - 1) * localNx + k * localNx * (def.Ny)];
		return;
	}

	
	if ((j == 0) && ((def.Ny) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + (j + 1) * localNx + k * localNx * (def.Ny)];
		return;

	/*
		int I = i_to_I(i, rank, size, def);
		if ((I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = 0.1;//def.S_w_gr;
		else
			HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = 0;
	*/
	}
	


	if ((k == 0) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + (k + 1) * localNx * (def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_w[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_w[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		return;
	}
}

void Border_Sg(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i - 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + (j - 1) * localNx + k * localNx * (def.Ny)];
		return;
	}


	if ((j == 0) && ((def.Ny) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + (j + 1) * localNx + k * localNx * (def.Ny)];
		return;
	/*
		int I = i_to_I(i, rank, size, def);
		if ((I >= (def.Nx) / 2 - (def.source)) && (I <= (def.Nx) / 2 + (def.source)) && (k >= (def.Nz) / 2 - (def.source)) && (k <= (def.Nz) / 2 + (def.source)))
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = 0.8; //def.S_g_gr;
		else
			HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = 0;
	*/
	}


	if ((k == 0) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + (k + 1) * localNx * (def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz) > 2))
	{
		HostArraysPtr.S_g[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.S_g[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		return;
	}
}


void Border_Pn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i == 0) && ((def.Nx) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + 1 + j * localNx + k * localNx * (def.Ny)]; 
		return;
	}

	if ((i == localNx - 1) && ((def.Nx) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i - 1 + j * localNx + k * localNx * (def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny) > 2))
	{
                HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + (j - 1) * localNx + k * localNx *(def.Ny)] 
				+ HostArraysPtr.ro_n[i + localNx * 1] * (def.g_const) * (def.hy);
		return;
	}

	if ((j==0) && ((def.Ny) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = def.P_atm;
		return;
	}

	if ((k == 0) && ((def.Nz) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + (k + 1) * localNx*(def.Ny)]; 
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz) > 2))
	{
		HostArraysPtr.P_n[i + j * localNx + k * localNx * (def.Ny)] = HostArraysPtr.P_n[i + j * localNx + (k - 1) * localNx * (def.Ny)];
		return;
	}
}