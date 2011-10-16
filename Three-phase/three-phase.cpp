#include "../defines.h"

void assign_Pn_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

	double Sn = (1. - HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)]);
	double S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_n_e = (Sn - def.S_nr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double k_w = pow(S_w_e,0.5)*pow(1.-pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_g = pow(S_g_e,0.5)*pow(1.-pow(1.-S_g_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_w = pow(1.-S_w_e,0.5)*pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_g = pow(S_n_e,0.5)*pow(1.-pow(1.-pow(S_n_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_n = S_n_e*k_n_w*k_n_g/(1-S_w_e)/(S_w_e+S_n_e);
	double P_k_nw = def.P_d_nw[media]*pow(pow(S_w_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);
	double P_k_gn = def.P_d_nw[media]*pow(pow(1.-S_g_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);


	HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]-P_k_nw;
	HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]+P_k_gn;
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;
	HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_g / def.mu_g;
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

	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = def.S_w_gr;
		else
			HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}

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

	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = def.S_g_gr;
		else
			HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}

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
		HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_n[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+localNx*1] * (def.g_const) * (def.hy); ;
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
