#include "../defines.h"
#include "b-l.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int media = HostArraysPtr.media[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
	double k_w, k_n;
	double S = 1 - HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];

	double S_wc=0.2;
	double S_wi=0.2;
	double S_or=0.2;
	double S_e=(S-S_wc)/(1-S_wc-S_or);

	k_w=S_e*S_e;
	k_n=(1-S_e)*(1-S_e);

	//krw(Sor) = kro(Swc) = 1.0

	if (S<S_wc)
	{
		k_w=0;
		k_n=1;
	}

	if (S>(1-S_or))
	{
		k_w=1;
		k_n=0;
	}

	HostArraysPtr.P_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
	HostArraysPtr.Xi_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = -1. * (HostArraysPtr.K[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = -1. * (HostArraysPtr.K[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)]) * k_n / def.mu_n;

//	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
}

void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i!=0) && (i!=(def.locNx)-1) && (j!=0) && (j!=(def.locNy)-1) && (((k!=0) && (k!=(def.locNz)-1)) || ((def.locNz)<2)))
	{
		double A1 = HostArraysPtr.roS_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
		double A2 = HostArraysPtr.roS_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
		double a = def.beta_w * (def.beta_n);
		double b = def.beta_w + def.beta_n - A2 * (def.beta_w) / (def.ro0_n) - A1 * (def.beta_n) / (def.ro0_w);
		double c = 1 - A2 / def.ro0_n  - A1 / def.ro0_w;
		double D = b*b - 4*a*c;
		double P1 = def.P_atm + (-1 * b + sqrt(D)) / (2 * a);
		double P2 = def.P_atm + (-1 * b - sqrt(D)) / (2 * a);
		
		if (P1 < 0)
			HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = P2;
		else
			HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = P1;

		HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = HostArraysPtr.roS_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] / (def.ro0_n * (1 + def.beta_n * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)]-def.P_atm)));
		/*
		if ( HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] < 0.2)
			HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = 0.2;

		if ( HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] > 0.8)
			HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = 0.8;
		*/	
		/*
		int media = HostArraysPtr.media[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
		double S_e, F1, F2, F1P, F2P, F1S, F2S, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
			S_e = (1 - HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.S_wr[media]) / (1 - def.S_wr[media]);

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)]) - HostArraysPtr.roS_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm)) * HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - HostArraysPtr.roS_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
			
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - (1 / det) * (F1P * F2 - F2P * F1);
		}  
		*/

		test_positive(HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
	}
	
}

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

	HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];

	test_S(HostArraysPtr.S_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
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

	// Без учета силы тяжести
	//HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
/*
	// В левом нижнем углу нагнетающая скважина
	if (((i==0) && (j==def.Ny-2)) || ((i==1) && (j==def.Ny-1)) || ((i==0) && (j==def.Ny-1)))
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 1e5;
*/
	
	if((j != 0) && (j != (def.locNy) - 1))
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	//else if(j == 0)
	//	HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	else
	{
		double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] 
		+ def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[(HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy])]) * (def.g_const) * (def.hy);

		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] + ro_g_dy;//HostArraysPtr.ro_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)] * (def.g_const) * (def.hy);
	}
		

	test_positive(HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)], __FILE__, __LINE__);
}

// Присвоение начальных условий
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


						HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy]=0;
						HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy]=0.5;
						//HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy] =0.3 + 0.3 * j / def.Ny;

						double ro_g_dy = (def.ro0_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] 
						+ def.ro0_w * (1 - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)])) * (def.m[(HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy])]) * (def.g_const) * (def.hy);

						//HostArraysPtr.ro_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
						if(j == 0)
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=def.P_atm;
						else
							HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=HostArraysPtr.P_w[i+(j-1)*def.locNx+k*def.locNx*def.locNy] + ro_g_dy;

						// Не учитывается сила тяжести
						//HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy]=def.P_atm;

						
						///!!!! Не учитываются капиллярные силы! Или надо считать перед этим шагом P_n
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
/*
						// В левом нижнем углу нагнетающая скважина
						if (((i==0) && (j==def.Ny-2)) || ((i==1) && (j==def.Ny-1)) || ((i==0) && (j==def.Ny-1)) || ((i==1) && (j==def.Ny-2)))
						{
							HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 1e6;
							//HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 0;
						}

						HostArraysPtr.ro_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));

						///!!!! Не учитываются капиллярные силы! Или надо считать перед этим шагом P_n
						HostArraysPtr.ro_n[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(def.locNx)+k*(def.locNx)*(def.locNy)] - def.P_atm));
						*/

						test_nan(HostArraysPtr.S_n[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
						test_nan(HostArraysPtr.P_w[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
						test_nan(HostArraysPtr.media[i+j*def.locNx+k*def.locNx*def.locNy], __FILE__, __LINE__);
					}
}