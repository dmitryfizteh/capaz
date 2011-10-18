#include "../defines.h"

//Вычисление значений давлений, плотностей и коэффициентов в законе Дарси в точке (i,j,k) среды media,
//исходя из известных значений основных параметров (Pn,Sw,Sg)
void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
//Запоминаем, с какой именно из сред работаем
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

//Вычисление значения насыщенности фазы n из условия равенства насыщенностей в сумме единице
	double Sn = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);
//Вычисление эффективных насыщенностей по формулам модели трехфазной фильтрации
	double S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_n_e = (Sn - def.S_nr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
//Вычисление относительных фазовых проницаемостей в соответствии с приближением Стоуна в модификации Азиза и Сеттари
	double k_w = pow(S_w_e,0.5)*pow(1.-pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_g = pow(S_g_e,0.5)*pow(1.-pow(1.-S_g_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_w = pow(1.-S_w_e,0.5)*pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_g = pow(S_n_e,0.5)*pow(1.-pow(1.-pow(S_n_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_n = S_n_e*k_n_w*k_n_g/(1-S_w_e)/(S_w_e+S_n_e);
//Вычисление капиллярных давлений в соответствии с приближенной моделью Паркера
	double P_k_nw = def.P_d_nw[media]*pow(pow(S_w_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);
	double P_k_gn = def.P_d_nw[media]*pow(pow(1.-S_g_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);

//Вычисление фазовых давлений c помощью капиллярных
	HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]-P_k_nw;
	HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]+P_k_gn;
//Вычисление коэффициентов закона Дарси
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;
	HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_g / def.mu_g;
}

//Вычисление (roS) фаз в точке (i,j,k) среды media из гиперболического уравнения непрерывности(модель Морозова)
void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_n = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);

		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * S_n;
		HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)];

		double divgrad1, divgrad2, divgrad3, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, Tx3, Ty3, Tz3, A1=0, A2=0, A3=0;

		if ((def.Nz)<2)
		{
			divgrad1=0;
			divgrad2=0;
			divgrad3=0;
			Tz1=0;
			Tz2=0;
			Tz3=0;
		}
		else
		{
			divgrad1 = (def.m[media] * (def.l_w) * (def.c) / 2.) * (HostArraysPtr.ro_w[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+(k+1)*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+(k-1)*localNx*(def.Ny)]) / ((def.hz) * (def.hz));
			divgrad2 = (def.m[media] * (def.l_n) * (def.c) / 2.) * (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * S_n - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * S_n + HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * S_n) / ((def.hz) * (def.hz));
			divgrad1 = (def.m[media] * (def.l_g) * (def.c) / 2.) * (HostArraysPtr.ro_g[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.S_g[i+j*localNx+(k+1)*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_g[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.S_g[i+j*localNx+(k-1)*localNx*(def.Ny)]) / ((def.hz) * (def.hz));
			Tz1 = (HostArraysPtr.ro_w[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.uz_w[i+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.uz_w[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));
			Tz2 = (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.uz_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.uz_n[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));		
			Tz3 = (HostArraysPtr.ro_g[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.uz_g[i+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_g[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.uz_g[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));		

		}

		divgrad1 += (def.m[media] * (def.l_w) * (def.c) / 2.) *
			((HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+1+j*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+(j+1)*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		divgrad2 += (def.m[media] * (def.l_n) * (def.c) / 2.) *
			((HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * S_n - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * S_n + HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * S_n) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * S_n - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * S_n + HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * S_n) / ((def.hy) * (def.hy)));
		
		divgrad3 += (def.m[media] * (def.l_g) * (def.c) / 2.) *
			((HostArraysPtr.ro_g[i+1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i+1+j*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_g[i-1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_g[i+(j+1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i+(j+1)*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_g[i+(j-1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_g[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		Tx1 = (HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty1 = (HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));
		Tx2 = (HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty2 = (HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));
		Tx3 = (HostArraysPtr.ro_g[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_g[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_g[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_g[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty3 = (HostArraysPtr.ro_g[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_g[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_g[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_g[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));
		
		if (t < 2 * (def.dt))
		{
			A1 = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad2 - Tx2 - Ty2 - Tz2);
			A3 = HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad3 - Tx3 - Ty3 - Tz3);
		}
		else
		{
			A1 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad1 - Tx1 - Ty1 - Tz1 + (2. * HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
			A2 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad2 - Tx2 - Ty2 - Tz2 + (2. * HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
			A3 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad1 - Tx3 - Ty3 - Tz3 + (2. * HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_g_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
		}
		HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_g_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = A1;
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = A2;
		HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)] = A3;
	}
}
//Решение системы 3*3 на основные параметры (Pn,Sw,Sg) методом Ньютона в точке (i,j,k) среды media
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_w_e, S_g_e, P_k_nw, P_k_gn, A, Sn, F1, F2, F3;
		double PkSw, PkSg, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sg, F2Sg, F3Sg, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
        //Вычисление эффективных насыщенностей
			S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
			S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
        //Переобозначение степенного коэффициента
			A = def.lambda[media]/(1.-def.lambda[media]);
        //Вычисление капиллярных давлений
			P_k_nw = def.P_d_nw[media]*pow((pow(S_w_e,A)-1.),1./def.lambda[media]);
			P_k_gn = (-1)*def.P_d_nw[media]*pow(pow((1.-S_g_e),A-1.),1./def.lambda[media]);
        //Вычисление насыщенности фазы n
			Sn = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);
        //Нахождение значения трех функций системы
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - P_k_nw - def.P_atm)) * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm)) * Sn - HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
			F3 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] + P_k_gn - def.P_atm)) * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)];
        //Вычисление частных производных производных капиллярных давлений по насыщенностям
			PkSw = def.P_d_nw[media]*(pow((S_w_e-1.),A)-1.)*pow(S_w_e,(A-1.))/(1.-def.lambda[media])/(1-def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
			PkSg = def.P_d_gn[media]*(pow((1.-S_g_e),A)-1.)*pow((1-S_g_e),(A-1.))/(1.-def.lambda[media])/(1-def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
        //Вычисление матрицы частных производных
			F1P = def.ro0_w * (def.beta_w) * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)];
			F2P = def.ro0_n * (def.beta_n) * Sn;
			F3P = def.ro0_g * (def.beta_g) * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)];
			F1Sw = (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - P_k_nw - def.P_atm - Sn*PkSw));
			F2Sw = F2Sg = (-1)*def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
			F3Sw = F1Sg = 0;
			F3Sg = (def.ro0_g) * (1 + (def.beta_g) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] + P_k_gn - def.P_atm + HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]*PkSg));
        //Вычисление детерминанта матрицы частных производных
			det = F1P * F2Sw * F3Sg - F1Sw * (F2P * F3Sg - F2Sg * F3P);
        //Получение решения системы методом Крамера в явном виде
			HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * (F2Sw * F3Sg * F1 - F1Sw * F3Sg * F2 + F1Sw * F2Sg * F3);
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * ((F2Sg * F3P - F2P * F3Sg) * F1 + F1P * F3Sg * F2 - F1P * F2Sg * F3);
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);
		}  
	}
}

void Border_Sw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+(j-1)*localNx+k*localNx*(def.Ny)];
		return;
	}
/*
	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = def.S_w_gr;
		else
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}
*/

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+j*localNx+(k+1)*localNx*(def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}

//Задание граничных условий отдельно для Sg,Sw,Pn
void Border_Sg(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+(j-1)*localNx+k*localNx*(def.Ny)];
		return;
	}

/*
	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = def.S_g_gr;
		else
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}
*/

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+j*localNx+(k+1)*localNx*(def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}


void Border_Pn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+1+j*localNx+k*localNx*(def.Ny)]; 
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
                HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i+localNx*1] * (def.g_const) * (def.hy);
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
		return;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+(k+1)*localNx*(def.Ny)]; 
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}
//Сохранение результатов, подлежит радикальному исправлению!!!
void print_plots_top (double t, consts def)
{
	char fname_Sn[30],fname_Sw[30],fname_Sg[30];
	FILE *fp_Sn,*fp_Sw,*fp_Sg;

	sprintf(fname_Sn,"plot_Sn/Sn=%012.4f.dat",t);
	sprintf(fname_Sw,"plot_Sw/Sw=%012.4f.dat",t);
	sprintf(fname_Sg,"plot_Sg/Sg=%012.4f.dat",t);

#ifdef _WIN32
	_mkdir("plot_Sn");
	_mkdir("plot_Sw");
	_mkdir("plot_Sg");
#else
	mkdir("plot_Sn",0000777);
	mkdir("plot_Sw",0000777);
	mkdir("plot_Sg",0000777);
#endif
	// Создание (или перезапись файлов) с графиками
	// 1. Для распределения насыщенностей NAPL S2
	// 2. Для распределения давлений воды P1
	// 3. Для распределения скоростей {u_x, u_y}
	// 4. Для среза насыщенностей S2 по оси Y
	// 5. Для среза насыщенностей S2 по оси X
	// 6. Для распределения типов грунтов
	if(!(fp_Sn=fopen(fname_Sn,"wt")) || !(fp_Sw=fopen(fname_Sw,"wt")) || !(fp_Sg=fopen(fname_Sg,"wt")))
		std::cout << "Not open file(s) in function SAVE_DATA_PLOTS! \n";

	//fprintf(fp_Sn,"TITLE =  \"Saturation of DNALP in time=%5.2f\" \n", t); // (1)
	//fprintf(fp_Sn,"VARIABLES = \"X\",\"Y\",\"Sn\" \n");
	//fprintf(fp_Sn,"ZONE T = \"BIG ZONE\", J=%d,I=%d, F = POINT\n", (def.NX), (def.Ny));

	fprintf(fp_Sn,"TITLE =  \"Saturation of DNALP in time=%5.2f\" \n", t); // (1)
	fprintf(fp_Sn,"VARIABLES = \"X\",\"Y\",\"Z\",\"Sn\",\"Sw\",\"Sg\"\n");
	fprintf(fp_Sn,"ZONE T = \"BIG ZONE\", K=%d,J=%d,I=%d, F = POINT\n", (def.Nx), (def.Nz), (def.Ny));

	fprintf(fp_Sw,"TITLE =  \"Saturation of water in time=%5.2f\" \n", t); // (2)
	fprintf(fp_Sw,"VARIABLES = \"X\",\"Y\",\"Z\",\"Sn\",\"Sw\",\"Sg\"\n");
	fprintf(fp_Sw,"ZONE T = \"BIG ZONE\", J=%d,I=%d, F = POINT\n", (def.Nx), (def.Nz), (def.Ny)); 

	fprintf(fp_Sg,"TITLE =  \"Saturation of gas in time=%5.2f\" \n", t); // (3)
	fprintf(fp_Sg,"VARIABLES = \"X\",\"Y\",\"Z\",\"Sn\",\"Sw\",\"Sg\"\n");
	fprintf(fp_Sg,"ZONE T = \"BIG ZONE\", J=%d,I=%d, F = POINT\n", (def.Nx), (def.Nz), (def.Ny));

	fclose(fp_Sn);
	fclose(fp_Sw);
	fclose(fp_Sg);
}

// Функция сохранения данных в файлы графиков (!3D)
void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, int size, int localNx, consts def)
{
	char fname_Sn[30],fname_Sw[30],fname_Sg[30];
	FILE *fp_Sn,*fp_Sw,*fp_Sg;
	int local;

	sprintf(fname_Sn,"plot_Sn/Sn=%012.4f.dat",t);
	sprintf(fname_Sw,"plot_Sw/Sw=%012.4f.dat",t);
	sprintf(fname_Sg,"plot_Sg/Sg=%012.4f.dat",t);
	
	// Открытие на дозапись и сохранение графиков

	if(!(fp_Sn=fopen(fname_Sn,"wt")) || !(fp_Sw=fopen(fname_Sw,"wt")) || !(fp_Sg=fopen(fname_Sg,"wt")))
		std::cout << "Not open file(s) in function SAVE_DATA_PLOTS! \n";

	for(int i=0; i<localNx; i++)
		for(int j=0; j<(def.Ny); j++)
			for(int k=0; k<(def.Nz); k++)
				if(is_active_point(i, localNx, rank, size))
				{
					local=i+j*localNx+k*localNx*(def.Ny);

                    fprintf(fp_Sn,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny)*(def.hy)-HostArraysPtr.y[local], 1.-HostArraysPtr.S_w[local]-HostArraysPtr.S_g[local], HostArraysPtr.P_n[local], HostArraysPtr.ux_n[local], HostArraysPtr.uz_n[local], (-1)*HostArraysPtr.uy_n[local]);
                    fprintf(fp_Sw,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny)*(def.hy)-HostArraysPtr.y[local], HostArraysPtr.S_w[local], HostArraysPtr.P_w[local], HostArraysPtr.ux_w[local], HostArraysPtr.uz_w[local], (-1)*HostArraysPtr.uy_w[local]);
                    fprintf(fp_Sg,"%.2e %.2e %.2e %.3e %.3e %.3e %.3e %.3e\n", HostArraysPtr.x[local], HostArraysPtr.z[local], (def.Ny)*(def.hy)-HostArraysPtr.y[local], HostArraysPtr.S_g[local], HostArraysPtr.P_g[local], HostArraysPtr.ux_g[local], HostArraysPtr.uz_g[local], (-1)*HostArraysPtr.uy_g[local]);
				}

	fclose(fp_Sn);
	fclose(fp_Sw);
	fclose(fp_Sg);
}
