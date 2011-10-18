#include "../defines.h"
#include "three-phase.h"

//���������� �������� ��������, ���������� � ������������� � ������ ����� � ����� (i,j,k) ����� media,
//������ �� ��������� �������� �������� ���������� (Pn,Sw,Sg)
void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
//����������, � ����� ������ �� ���� ��������
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

//���������� �������� ������������ ���� n �� ������� ��������� ������������� � ����� �������
	double Sn = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);
//���������� ����������� ������������� �� �������� ������ ���������� ����������
	double S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_n_e = (Sn - def.S_nr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
//���������� ������������� ������� �������������� � ������������ � ������������ ������ � ����������� ����� � �������
	double k_w = pow(S_w_e,0.5)*pow(1.-pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_g = pow(S_g_e,0.5)*pow(1.-pow(1.-S_g_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_w = pow(1.-S_w_e,0.5)*pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_g = pow(S_n_e,0.5)*pow(1.-pow(1.-pow(S_n_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_n = S_n_e*k_n_w*k_n_g/(1-S_w_e)/(S_w_e+S_n_e);
//���������� ����������� �������� � ������������ � ������������ ������� �������
	double P_k_nw = def.P_d_nw[media]*pow(pow(S_w_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);
	double P_k_gn = def.P_d_nw[media]*pow(pow(1.-S_g_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);

//���������� ������� �������� c ������� �����������
	HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]-P_k_nw;
	HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]+P_k_gn;
//���������� ������������� ������ �����
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;
	HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_g / def.mu_g;
}

//������� ������� 3*3 �� �������� ��������� (Pn,Sw,Sg) ������� ������� � ����� (i,j,k) ����� media
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_w_e, S_g_e, P_k_nw, P_k_gn, A, Sn, F1, F2, F3;
		double PkSw, PkSg, F1P, F2P, F3P, F1Sw, F2Sw, F3Sw, F1Sg, F2Sg, F3Sg, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
        //���������� ����������� �������������
			S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
			S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
        //��������������� ���������� ������������
			A = def.lambda[media]/(1.-def.lambda[media]);
        //���������� ����������� ��������
			P_k_nw = def.P_d_nw[media]*pow((pow(S_w_e,A)-1.),1./def.lambda[media]);
			P_k_gn = (-1)*def.P_d_nw[media]*pow(pow((1.-S_g_e),A-1.),1./def.lambda[media]);
        //���������� ������������ ���� n
			Sn = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);
        //���������� �������� ���� ������� �������
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - P_k_nw - def.P_atm)) * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm)) * Sn - HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
			F3 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] + P_k_gn - def.P_atm)) * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)];
        //���������� ������� ����������� ����������� ����������� �������� �� �������������
			PkSw = def.P_d_nw[media]*(pow((S_w_e-1.),A)-1.)*pow(S_w_e,(A-1.))/(1.-def.lambda[media])/(1-def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
			PkSg = def.P_d_gn[media]*(pow((1.-S_g_e),A)-1.)*pow((1-S_g_e),(A-1.))/(1.-def.lambda[media])/(1-def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
        //���������� ������� ������� �����������
			F1P = def.ro0_w * (def.beta_w) * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)];
			F2P = def.ro0_n * (def.beta_n) * Sn;
			F3P = def.ro0_g * (def.beta_g) * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)];
			F1Sw = (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - P_k_nw - def.P_atm - Sn*PkSw));
			F2Sw = F2Sg = (-1)*def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
			F3Sw = F1Sg = 0;
			F3Sg = (def.ro0_g) * (1 + (def.beta_g) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] + P_k_gn - def.P_atm + HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]*PkSg));
        //���������� ������������ ������� ������� �����������
			det = F1P * F2Sw * F3Sg - F1Sw * (F2P * F3Sg - F2Sg * F3P);
        //��������� ������� ������� ������� ������� � ����� ����
			HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * (F2Sw * F3Sg * F1 - F1Sw * F3Sg * F2 + F1Sw * F2Sg * F3);
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * ((F2Sg * F3P - F2P * F3Sg) * F1 + F1P * F3Sg * F2 - F1P * F2Sg * F3);
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - (1. / det) * (F3P * F1Sw * F2 - F3P * F2Sw * F1 + (F1P*F2Sw - F2P * F1Sw) * F3);
		}  
	}
}

void Border(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	Border_Sw(HostArraysPtr, i, j, k, localNx, rank, size, def);
	Border_Sg(HostArraysPtr, i, j, k, localNx, rank, size, def);
	Border_Pn(HostArraysPtr, i, j, k, localNx, def);
	return;
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

//������� ��������� ������� �������� ��� Sg,Sw,Pn
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

// ���������� ��������� ������ �� ���� ������
void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def)
{
	*t=0;
	for(int i=0;i<localNx;i++)
		for(int j=0;j<localNy;j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
					{
						// �������������� ��������� ��������� ���������� � ����������
						int I=i_to_I(i,rank,size,def);
						// ���� ����� �� ������� �������, �� ����� (def.source) ����� �� ������,
						// �� � ��� ��������� ������������. �����, �������
/*
						if ((j==0) && (I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
						{
							HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)]=def.S_w_gr;
							HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]=def.S_g_gr;
						}
						else
						{
							HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)]=0;
							HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]=0;
						}
*/
						HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]=def.P_atm+j * (def.ro0_n) * (def.g_const)*(def.hy);
						HostArraysPtr.x[i+j*localNx+k*localNx*(def.Ny)]=I*(def.hx);
						HostArraysPtr.y[i+j*localNx+k*localNx*(def.Ny)]=j*(def.hy);
						HostArraysPtr.z[i+j*localNx+k*localNx*(def.Ny)]=k*(def.hz);

						HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)]=0;
					}
}

//����� �����������
void print_plots_top (double t, consts def)
{
	return;
}

void print_plots(ptr_Arrays HostArraysPtr, double t, int rank, int size, int localNx, consts def)
{
	return;
}