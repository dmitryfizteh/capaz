#include "../defines.h"
#include "b-l.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
	double S_e = (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]);
	double k_w = pow(S_e, (2. + 3. * (def.lambda[media])) / def.lambda[media]);
	double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));
	double P_k = def.P_d[media] * pow((1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]), -1. / def.lambda[media]);

	HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k;
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;
}

void Border(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	Border_Sn(HostArraysPtr, i, j, k, localNx, rank, size, def);
	Border_Pw(HostArraysPtr, i, j, k, localNx, def);
	return;
}

void Border_Sn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = def.S_n_gr;
		else
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}

void Border_Pw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)]; 
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+localNx*1] * (def.g_const) * (def.hy); ;
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
		return;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k+1)*localNx*(def.Ny)]; 
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}