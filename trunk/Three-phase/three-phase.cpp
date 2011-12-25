#include "../defines.h"
#include "three-phase.h"

//������� ���������� �������� ��������, ���������� � ������������� � ������ ����� � ����� (i,j,k) ����� media,
//������ �� ��������� �������� �������� ���������� (Pn,Sw,Sg)
//1. ����������, � ����� ������ �� ���� ��������
//2. ���������� �������� ������������ ���� n �� ������� ��������� ������������� � ����� �������
//3. ���������� ����������� ������������� �� �������� ������ ���������� ����������
//4. ���������� ������������� ������� �������������� � ������������ � ������������ ������ � ����������� ����� � �������
//5. ���������� ����������� �������� � ������������ � ������������ ������� �������
//6. ���������� ������� �������� c ������� �����������
//7. ���������� ������������� ������ �����
void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{

	int media = HostArraysPtr.media[i + j * (locN.x) + k * (locN.x) * (locN.y)];                                                                                                                                                                            
	double k_w, k_g, k_n, P_k_nw, P_k_gn;
	double A = def.lambda[media]; 
	double S_w_e = (HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);                                                                                                                            
	double S_g_e = (HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]); 
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

	HostArraysPtr.P_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] - P_k_nw;                                                             
	HostArraysPtr.P_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] + P_k_gn;                                                             

	HostArraysPtr.Xi_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = (-1.) * (def.K[media]) * k_w / def.mu_w;                                                                                                                 
	HostArraysPtr.Xi_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = (-1.) * (def.K[media]) * k_n / def.mu_n;                                                                                                                 
	HostArraysPtr.Xi_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] = (-1.) * (def.K[media]) * k_g / def.mu_g;                                                                                                                 

	test_positive(HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
}

//������� ������� ������� 3*3 �� �������� ��������� (Pn,Sw,Sg) ������� ������� � ����� (i,j,k) ����� media
//1. ���������� ����������� �������������
//2. ��������������� ���������� ������������
//3. ���������� ����������� ��������
//4. ���������� ������������ ���� n
//5. ���������� �������� ���� ������� �������
//6. ���������� ������� ����������� ����������� ����������� �������� �� �������������
//7. ���������� ������� ������� �����������
//8. ���������� ������������ ������� ������� �����������
//9. ��������� ������� ������� ������� ������� � ����� ����
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	if ((i != 0) && (i != (locN.x) - 1) && (j != 0) && (j != (locN.y) - 1) && (((k != 0) && (k != (locN.z) - 1)) || ((locN.z) < 2)))
	{
		int media = HostArraysPtr.media[i + j * (locN.x) + k * (locN.x) * (locN.y)];
		double S_w_e, S_g_e, P_k_nw, P_k_gn, A, Sn, F1, F2, F3;
		double PkSw, PkSg, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sg, F2Sg, F3Sg, det;

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = (HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);   
			S_g_e = (HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.S_gr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);  

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

			Sn = 1. - HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] - HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)];                     

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] - P_k_nw - def.P_atm))                                                               
				* HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] - HostArraysPtr.roS_w[i + j * (locN.x) + k * (locN.x) * (locN.y)];                  
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.P_atm)) * Sn                                                                  
				- HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
			F3 = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] + P_k_gn - def.P_atm))                                                              
				* HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] - HostArraysPtr.roS_g[i + j * (locN.x) + k * (locN.x) * (locN.y)];

			F1P = def.ro0_w * def.beta_w * HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)];                                                                                                                             
			F2P = def.ro0_n * def.beta_n * Sn;                                                                                                                                                                                                                                    
			F3P = def.ro0_g * def.beta_g * HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)];                                                                                                                                                
			F1Sw = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] - P_k_nw - def.P_atm
				- HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] * PkSw));                                                                                                                                                         
			F2Sw = F2Sg = (-1) * def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] - def.P_atm));                                         
			F3Sw = F1Sg = 0;                                                                                                                                                                                                                                                                     
			F3Sg = def.ro0_g * (1 + def.beta_g * (HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] + P_k_gn - def.P_atm                                                                  
				+ HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] * PkSg));

			det = F1P * F2Sw * F3Sg - F1Sw * (F2P * F3Sg - F2Sg * F3P);                                                                                                                                                                                         

			HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)]                                                              
			- (1. / det) * (F2Sw * F3Sg * F1 - F1Sw * F3Sg * F2 + F1Sw * F2Sg * F3);
			HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)]                                                              
			- (1. / det) * ((F2Sg * F3P - F2P * F3Sg) * F1 + F1P * F3Sg * F2 - F1P * F2Sg * F3);
			HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)]                                                              
			- (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);

			test_S(HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
			test_S(HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
			test_positive(HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		}  
	}
}

//������� ��������� ������� �������� ��� (Sw,Sg),Pn

// ������� ��������� ������� � ������� ������ ��������, �� � ��������� �������������� ����������
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, int rank, parts_sizes parts, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (locN.x) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (locN.y) - 1)
		j1 --;
	if((k == 0) && ((locN.z) > 2))
		k1 ++;
	if((k == (locN.z) - 1) && ((locN.z) > 2))
		k1 --;

	if((j != 0) || ((def.source) <= 0))
	{
		HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.S_w[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)];
		HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.S_g[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)];
	}

	if((j == 0) && ((def.source) > 0))
	{
		HostArraysPtr.S_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = def.S_w_gr;
		HostArraysPtr.S_g[i + j * (locN.x) + k * (locN.x) * (locN.y)] = def.S_g_gr;
	}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (locN.x) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (locN.y) - 1)
		j1 --;
	if((k == 0) && ((locN.z) > 2))
		k1 ++;
	if((k == (locN.z) - 1) && ((locN.z) > 2))
		k1 --;

	if((j != 0) && (j != (locN.y) - 1))
		HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_n[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)];
	else if(j == 0)
		HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = def.P_atm;
	else
		HostArraysPtr.P_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_n[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)] + ro_eff_gdy(HostArraysPtr, i1, j1, k1, locN, def);
}
